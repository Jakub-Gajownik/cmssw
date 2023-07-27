// Check that ALPAKA_HOST_ONLY is not defined during device compilation:
#ifdef ALPAKA_HOST_ONLY
#error ALPAKA_HOST_ONLY defined in device compilation
#endif

#include <cmath>
#include <limits>
#include <alpaka/alpaka.hpp>

#include "CondFormats/EcalObjects/interface/EcalPulseCovariances.h"
#include "DataFormats/CaloRecHit/interface/MultifitComputations.h"
#include "FWCore/Utilities/interface/CMSUnrollLoop.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "AmplitudeComputationKernels.h"
#include "KernelHelpers.h"
#include "EcalUncalibRecHitMultiFitAlgoPortable.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace ecal {
    namespace multifit {

      using namespace ::ecal::multifit;

      template <typename MatrixType>
      ALPAKA_FN_ACC ALPAKA_FN_INLINE void update_covariance(EcalPulseCovariance const& pulse_covariance,
                                                            MatrixType& inverse_cov,
                                                            SampleVector const& amplitudes) {
        constexpr auto nsamples = SampleVector::RowsAtCompileTime;
        constexpr auto npulses = BXVectorType::RowsAtCompileTime;

        CMS_UNROLL_LOOP
        for (unsigned int ipulse = 0; ipulse < npulses; ++ipulse) {
          auto const amplitude = amplitudes.coeff(ipulse);
          if (amplitude == 0)
            continue;

          // FIXME: ipulse - 5 -> ipulse - firstOffset
          int bx = ipulse - 5;
          int first_sample_t = std::max(0, bx + 3);
          int offset = -3 - bx;

          auto const value_sq = amplitude * amplitude;

          for (int col = first_sample_t; col < nsamples; ++col) {
            for (int row = col; row < nsamples; ++row) {
              inverse_cov(row, col) += value_sq * pulse_covariance.covval[row + offset][col + offset];
            }
          }
        }
      }

      ///
      /// launch ctx parameters are (nchannels / block, blocks)
      /// TODO: trivial impl for now, there must be a way to improve
      ///
      /// Conventions:
      ///   - amplitudes -> solution vector, what we are fitting for
      ///   - samples -> raw detector responses
      ///   - passive constraint - satisfied constraint
      ///   - active constraint - unsatisfied (yet) constraint
      ///
      class kernel_minimize {
        public:
          template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
          ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                        v1::InputProduct::ConstView const & digisDevEB,
                                        v1::InputProduct::ConstView const & digisDevEE,
                                        v1::OutputProduct::View uncalibRecHitsEB,
                                        v1::OutputProduct::View uncalibRecHitsEE,
                                        EcalMultifitConditionsPortableDevice::ConstView conditionsDev,
                                        ::ecal::multifit::SampleMatrix const* noisecov,
                                        ::ecal::multifit::PulseMatrixType const* pulse_matrix,
                                        ::ecal::multifit::BXVectorType* bxs,
                                        ::ecal::multifit::SampleVector const* samples,
                                        bool* hasSwitchToGain6,
                                        bool* hasSwitchToGain1,
                                        bool* isSaturated,
                                        char* acState,
                                        int max_iterations) const {
            // FIXME: ecal has 10 samples and 10 pulses....
            // but this needs to be properly treated and renamed everywhere
            constexpr auto NSAMPLES = SampleMatrix::RowsAtCompileTime;
            constexpr auto NPULSES = SampleMatrix::ColsAtCompileTime;
            static_assert(NSAMPLES == NPULSES);

            using DataType = SampleVector::Scalar;

            DataType* shrmem = alpaka::getDynSharedMem<DataType>(acc);

            auto const blockDim(alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]);
            auto const threadIdx(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
            auto const blockIdx(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);

            auto const nchannelsEB = digisDevEB.size();
            auto const nchannels = nchannelsEB + digisDevEE.size();
            auto const offsetForHashes = nchannelsEB; // FIXME: Check this

            auto const* pulse_covariance = reinterpret_cast<const EcalPulseCovariance*>(conditionsDev.pulseCovariance());

            DataType* shrMatrixLForFnnlsStorage = shrmem + calo::multifit::MapSymM<DataType, NPULSES>::total * threadIdx;
            DataType* shrAtAStorage = shrmem + calo::multifit::MapSymM<DataType, NPULSES>::total * (threadIdx + blockDim);

            // channel
            const auto idx = threadIdx + blockDim * blockIdx;

            auto* amplitudes = reinterpret_cast<SampleVector*>(idx >= nchannelsEB ? uncalibRecHitsEE.outOfTimeAmplitudes()->array.data() : uncalibRecHitsEB.outOfTimeAmplitudes()->array.data());
            auto* energies = idx >= nchannelsEB ? uncalibRecHitsEE.amplitude() : uncalibRecHitsEB.amplitude();
            auto* chi2s = idx >= nchannelsEB ? uncalibRecHitsEE.chi2() : uncalibRecHitsEB.chi2();

            if (idx < nchannels) {
              if (static_cast<MinimizationState>(acState[idx]) == MinimizationState::Precomputed)
                return;

              // get the hash
              int const inputCh = idx >= nchannelsEB ? idx - nchannelsEB : idx;
              auto const* dids = idx >= nchannelsEB ? digisDevEE.id() : digisDevEB.id();
              auto const did = DetId{dids[inputCh]};
              auto const isBarrel = did.subdetId() == EcalBarrel;
              auto const hashedId = isBarrel ? ecal::reconstruction::hashedIndexEB(did.rawId())
                                             : offsetForHashes + ecal::reconstruction::hashedIndexEE(did.rawId());

              // inits
              int iter = 0;
              int npassive = 0;

              calo::multifit::ColumnVector<NPULSES, int> pulseOffsets;
              CMS_UNROLL_LOOP
              for (int i = 0; i < NPULSES; ++i)
                pulseOffsets(i) = i;

              calo::multifit::ColumnVector<NPULSES, DataType> resultAmplitudes;
              CMS_UNROLL_LOOP
              for (int counter = 0; counter < NPULSES; ++counter)
                resultAmplitudes(counter) = 0;

              // inits
              //SampleDecompLLT covariance_decomposition;
              //SampleMatrix inverse_cov;
              //        SampleVector::Scalar chi2 = 0, chi2_now = 0;
              float chi2 = 0, chi2_now = 0;

              // loop until converge
              while (true) {
                if (iter >= max_iterations)
                  break;

                //inverse_cov = noisecov[idx];
                //DataType covMatrixStorage[MapSymM<DataType, NSAMPLES>::total];
                DataType* covMatrixStorage = shrMatrixLForFnnlsStorage;
                calo::multifit::MapSymM<DataType, NSAMPLES> covMatrix{covMatrixStorage};
                int counter = 0;
                CMS_UNROLL_LOOP
                for (int col = 0; col < NSAMPLES; ++col) {
                  CMS_UNROLL_LOOP
                  for (int row = col; row < NSAMPLES; ++row) {
                    covMatrixStorage[counter++] = noisecov[idx].coeffRef(row, col);
                  }
                }
                update_covariance(pulse_covariance[hashedId], covMatrix, resultAmplitudes);

                // compute actual covariance decomposition
                //covariance_decomposition.compute(inverse_cov);
                //auto const& matrixL = covariance_decomposition.matrixL();
                DataType matrixLStorage[calo::multifit::MapSymM<DataType, NSAMPLES>::total];
                calo::multifit::MapSymM<DataType, NSAMPLES> matrixL{matrixLStorage};
                calo::multifit::compute_decomposition_unrolled(matrixL, covMatrix);

                // L * A = P
                calo::multifit::ColMajorMatrix<NSAMPLES, NPULSES> A;
                calo::multifit::solve_forward_subst_matrix(A, pulse_matrix[idx], matrixL);

                // L b = s
                float reg_b[NSAMPLES];
                calo::multifit::solve_forward_subst_vector(reg_b, samples[idx], matrixL);

                // FIXME: shared mem
                //DataType AtAStorage[MapSymM<DataType, NPULSES>::total];
                calo::multifit::MapSymM<DataType, NPULSES> AtA{shrAtAStorage};
                //SampleMatrix AtA;
                SampleVector Atb;
                CMS_UNROLL_LOOP
                for (int icol = 0; icol < NPULSES; ++icol) {
                  float reg_ai[NSAMPLES];

                  // load column icol
                  CMS_UNROLL_LOOP
                  for (int counter = 0; counter < NSAMPLES; ++counter)
                    reg_ai[counter] = A(counter, icol);

                  // compute diagoanl
                  float sum = 0.f;
                  CMS_UNROLL_LOOP
                  for (int counter = 0; counter < NSAMPLES; ++counter)
                    sum += reg_ai[counter] * reg_ai[counter];

                  // store
                  AtA(icol, icol) = sum;

                  // go thru the other columns
                  CMS_UNROLL_LOOP
                  for (int j = icol + 1; j < NPULSES; ++j) {
                    // load column j
                    float reg_aj[NSAMPLES];
                    CMS_UNROLL_LOOP
                    for (int counter = 0; counter < NSAMPLES; ++counter)
                      reg_aj[counter] = A(counter, j);

                    // accum
                    float sum = 0.f;
                    CMS_UNROLL_LOOP
                    for (int counter = 0; counter < NSAMPLES; ++counter)
                      sum += reg_aj[counter] * reg_ai[counter];

                    // store
                    //AtA(icol, j) = sum;
                    AtA(j, icol) = sum;
                  }

                  // Atb accum
                  float sum_atb = 0.f;
                  CMS_UNROLL_LOOP
                  for (int counter = 0; counter < NSAMPLES; ++counter)
                    sum_atb += reg_ai[counter] * reg_b[counter];

                  // store atb
                  Atb(icol) = sum_atb;
                }

                // FIXME: shared mem
                //DataType matrixLForFnnlsStorage[MapSymM<DataType, NPULSES>::total];
                calo::multifit::MapSymM<DataType, NPULSES> matrixLForFnnls{shrMatrixLForFnnlsStorage};

                calo::multifit::fnnls(AtA,
                                      Atb,
                                      //amplitudes[idx],
                                      resultAmplitudes,
                                      npassive,
                                      pulseOffsets,
                                      matrixLForFnnls,
                                      1e-11,
                                      500,
                                      16,
                                      2);

                calo::multifit::calculateChiSq(matrixL, pulse_matrix[idx], resultAmplitudes, samples[idx], chi2_now);

                auto const deltachi2 = chi2_now - chi2;
                chi2 = chi2_now;

                if (std::abs(deltachi2) < 1e-3)
                  break;

                //---- AM: TEST
                //---- it was 3 lines above, now here as in the CPU version
                ++iter;
              }

              // store to global output values
              // FIXME: amplitudes are used in global directly
              chi2s[inputCh] = chi2;
              energies[inputCh] = resultAmplitudes(5);

              CMS_UNROLL_LOOP
              for (int counter = 0; counter < NPULSES; ++counter)
                amplitudes[inputCh](counter) = resultAmplitudes(counter);
            }
          }
      };

      namespace v1 {

        void minimization_procedure(InputProduct const& digisDevEB,
                                    InputProduct const& digisDevEE,
                                    OutputProduct& uncalibRecHitsDevEB,
                                    OutputProduct& uncalibRecHitsDevEE,
                                    EventDataForScratchGPU& scratch,
                                    EcalMultifitConditionsPortableDevice const& conditionsDev,
                                    ConfigurationParameters const& configParams,
                                    Queue& queue) {
          using DataType = SampleVector::Scalar;
          auto const totalChannels = static_cast<uint32_t>(digisDevEB->metadata().size() + digisDevEE->metadata().size());  // FIXME check this
          // TODO: configure from python
          auto threads_min = configParams.kernelMinimizeThreads[0];
          auto blocks_min = threads_min > totalChannels ? 1 : (totalChannels + threads_min - 1) / threads_min;
          //uint32_t const offsetForHashes = conditionsDev.offsetForHashes;  // FIXME check this

          auto workDivMinimize = cms::alpakatools::make_workdiv<Acc1D>(blocks_min, threads_min);
          alpaka::exec<Acc1D>(queue,
                              workDivMinimize,
                              kernel_minimize{},
                              digisDevEB.const_view(),
                              digisDevEE.const_view(),
                              uncalibRecHitsDevEB.view(),
                              uncalibRecHitsDevEE.view(),
                              conditionsDev.const_view(),
                              reinterpret_cast<::ecal::multifit::SampleMatrix*>(scratch.noisecovDevBuf.value().data()),
                              reinterpret_cast<::ecal::multifit::PulseMatrixType*>(scratch.pulse_matrixDevBuf.value().data()),
                              reinterpret_cast<::ecal::multifit::BXVectorType*>(scratch.activeBXsDevBuf.value().data()),
                              reinterpret_cast<::ecal::multifit::SampleVector*>(scratch.samplesDevBuf.value().data()),
                              scratch.hasSwitchToGain6DevBuf.value().data(),
                              scratch.hasSwitchToGain1DevBuf.value().data(),
                              scratch.isSaturatedDevBuf.value().data(),
                              scratch.acStateDevBuf.value().data(),
                              50);

        }

      }  // namespace v1
    }  // namespace multifit
  }  // namespace ecal
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

