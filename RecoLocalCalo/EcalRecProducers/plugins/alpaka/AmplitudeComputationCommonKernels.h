#ifndef RecoLocalCalo_EcalRecProducers_plugins_alpaka_AmplitudeComputationCommonKernels_h
#define RecoLocalCalo_EcalRecProducers_plugins_alpaka_AmplitudeComputationCommonKernels_h

#include <cstdlib>
#include <limits>
#include <alpaka/alpaka.hpp>

#include "CondFormats/EcalObjects/interface/alpaka/EcalMultifitConditionsPortable.h"
#include "DataFormats/EcalDigi/interface/alpaka/EcalDigiDeviceCollection.h"
#include "DataFormats/EcalRecHit/interface/alpaka/EcalUncalibratedRecHitDeviceCollection.h"
//#include "CondFormats/EcalObjects/interface/EcalPulseCovariances.h"
#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"
//#include "CondFormats/EcalObjects/interface/EcalSamplesCorrelation.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "DataFormats/EcalDigi/interface/EcalMGPASample.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
//#include "DataFormats/Math/interface/approx_exp.h"
//#include "DataFormats/Math/interface/approx_log.h"
#include "FWCore/Utilities/interface/CMSUnrollLoop.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"

#include "DeclsForKernels.h"
#include "../EigenMatrixTypes_gpu.h"
#include "KernelHelpers.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace ecal {
    namespace multifit {

      using namespace cms::alpakatools;
      ///
      /// assume kernel launch configuration is
      /// (MAXSAMPLES * nchannels, blocks)
      /// TODO: is there a point to split this kernel further to separate reductions
      ///
      class kernel_prep_1d_and_initialize {
        public:
          template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
          ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                        DigiDeviceCollection::ConstView digisDevEB,
                                        DigiDeviceCollection::ConstView digisDevEE,
                                        UncalibratedRecHitDeviceCollection::View uncalibRecHitsEB,
                                        UncalibratedRecHitDeviceCollection::View uncalibRecHitsEE,
                                        EcalMultifitConditionsPortableDevice::ConstView conditionsDev,
                                        ::ecal::multifit::SampleVector* amplitudes,
                                        ::ecal::multifit::SampleGainVector* gainsNoise,
                                        bool* hasSwitchToGain6,
                                        bool* hasSwitchToGain1,
                                        bool* isSaturated,
                                        char* acState,
                                        ::ecal::multifit::BXVectorType* bxs,
                                        bool const gainSwitchUseMaxSampleEB,
                                        bool const gainSwitchUseMaxSampleEE) const {
            constexpr bool dynamicPedestal = false;  //---- default to false, ok
            constexpr auto nsamples = EcalDataFrame::MAXSAMPLES;
            constexpr int sample_max = 5;
            constexpr int full_pulse_max = 9;
            auto const offsetForHashes = conditionsDev.offsetEE();

            auto const nchannelsEB = digisDevEB.size();
            auto const nchannels = nchannelsEB + digisDevEE.size();
            auto const threadIdx = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u];
            auto const blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];
            auto const blockDim = alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u];
            auto const tx = threadIdx + blockIdx * blockDim;
            auto const nchannels_per_block = blockDim / nsamples;
            auto const ch = tx / nsamples;
            // for accessing input arrays
            int const inputCh = ch >= nchannelsEB ? ch - nchannelsEB : ch;
            int const inputTx = ch >= nchannelsEB ? tx - nchannelsEB * nsamples : tx;
            // eb is first and then ee
            auto const* digis_in = ch >= nchannelsEB ? digisDevEE.data()->array.data() : digisDevEB.data()->array.data(); // is there a way to use digisDevEE.data()->data() instead?
            auto const* dids = ch >= nchannelsEB ? digisDevEE.id() : digisDevEB.id();
            auto const sample = threadIdx % nsamples;

            auto* amplitudesForMinimization = reinterpret_cast<::ecal::multifit::SampleVector*>(ch >= nchannelsEB ? uncalibRecHitsEE.outOfTimeAmplitudes()->array.data() : uncalibRecHitsEB.outOfTimeAmplitudes()->array.data());
            auto* energies = ch >= nchannelsEB ? uncalibRecHitsEE.amplitude() : uncalibRecHitsEB.amplitude();
            auto* chi2 = ch >= nchannelsEB ? uncalibRecHitsEE.chi2() : uncalibRecHitsEB.chi2();
            auto* g_pedestal = ch >= nchannelsEB ? uncalibRecHitsEE.pedestal() : uncalibRecHitsEB.pedestal();
            auto* dids_out = ch >= nchannelsEB ? uncalibRecHitsEE.id() : uncalibRecHitsEB.id();
            auto* flags = ch >= nchannelsEB ? uncalibRecHitsEE.flags() : uncalibRecHitsEB.flags();

            auto const* shapes_in = reinterpret_cast<const EcalPulseShape*>(conditionsDev.pulseShapes()->data());

            if (ch < nchannels) {
              char* shared_mem = alpaka::getDynSharedMem<char>(acc);
              auto* shr_hasSwitchToGain6 = reinterpret_cast<bool*>(shared_mem);
              auto* shr_hasSwitchToGain1 = shr_hasSwitchToGain6 + nchannels_per_block * nsamples;
              auto* shr_hasSwitchToGain0 = shr_hasSwitchToGain1 + nchannels_per_block * nsamples;
              auto* shr_isSaturated = shr_hasSwitchToGain0 + nchannels_per_block * nsamples;
              auto* shr_hasSwitchToGain0_tmp = shr_isSaturated + nchannels_per_block * nsamples;
              auto* shr_counts = reinterpret_cast<char*>(shr_hasSwitchToGain0_tmp) + nchannels_per_block * nsamples;

              //
              // indices
              //
              auto const did = DetId{dids[inputCh]};
              auto const isBarrel = did.subdetId() == EcalBarrel;
              // TODO offset for ee, 0 for eb
              auto const hashedId = isBarrel ? reconstruction::hashedIndexEB(did.rawId())
                                             : offsetForHashes + reconstruction::hashedIndexEE(did.rawId());

              //
              // pulse shape template

              // will be used in the future for setting state
              auto const rmsForChecking = conditionsDev.pedestals_rms_x12()[hashedId];

              //
              // amplitudes
              //
              auto const adc = ecalMGPA::adc(digis_in[inputTx]);
              auto const gainId = ecalMGPA::gainId(digis_in[inputTx]);
              ::ecal::multifit::SampleVector::Scalar amplitude = 0.;
              ::ecal::multifit::SampleVector::Scalar pedestal = 0.;
              ::ecal::multifit::SampleVector::Scalar gainratio = 0.;

              // store into shared mem for initialization
              shr_hasSwitchToGain6[threadIdx] = gainId == EcalMgpaBitwiseGain6;
              shr_hasSwitchToGain1[threadIdx] = gainId == EcalMgpaBitwiseGain1;
              shr_hasSwitchToGain0_tmp[threadIdx] = gainId == EcalMgpaBitwiseGain0;
              shr_hasSwitchToGain0[threadIdx] = shr_hasSwitchToGain0_tmp[threadIdx];
              shr_counts[threadIdx] = 0;

              alpaka::syncBlockThreads(acc);

              // non-divergent branch (except for the last 4 threads)
              if (threadIdx <= blockDim - 5) {
                CMS_UNROLL_LOOP
                for (int i = 0; i < 5; ++i)
                  shr_counts[threadIdx] += shr_hasSwitchToGain0[threadIdx + i];
              }
              shr_isSaturated[threadIdx] = shr_counts[threadIdx] == 5;

              //
              // unrolled reductions
              //
              if (sample < 5) {
                shr_hasSwitchToGain6[threadIdx] =
                    shr_hasSwitchToGain6[threadIdx] || shr_hasSwitchToGain6[threadIdx + 5];
                shr_hasSwitchToGain1[threadIdx] =
                    shr_hasSwitchToGain1[threadIdx] || shr_hasSwitchToGain1[threadIdx + 5];

                // duplication of hasSwitchToGain0 in order not to
                // introduce another syncthreads
                shr_hasSwitchToGain0_tmp[threadIdx] =
                    shr_hasSwitchToGain0_tmp[threadIdx] || shr_hasSwitchToGain0_tmp[threadIdx + 5];
              }

              alpaka::syncBlockThreads(acc);

              if (sample < 2) {
                // note, both threads per channel take value [3] twice to avoid another if
                shr_hasSwitchToGain6[threadIdx] = shr_hasSwitchToGain6[threadIdx] ||
                                                  shr_hasSwitchToGain6[threadIdx + 2] ||
                                                  shr_hasSwitchToGain6[threadIdx + 3];
                shr_hasSwitchToGain1[threadIdx] = shr_hasSwitchToGain1[threadIdx] ||
                                                  shr_hasSwitchToGain1[threadIdx + 2] ||
                                                  shr_hasSwitchToGain1[threadIdx + 3];

                shr_hasSwitchToGain0_tmp[threadIdx] = shr_hasSwitchToGain0_tmp[threadIdx] ||
                                                      shr_hasSwitchToGain0_tmp[threadIdx + 2] ||
                                                      shr_hasSwitchToGain0_tmp[threadIdx + 3];

                // sample < 2 -> first 2 threads of each channel will be used here
                // => 0 -> will compare 3 and 4 and put into 0
                // => 1 -> will compare 4 and 5 and put into 1
                shr_isSaturated[threadIdx] = shr_isSaturated[threadIdx + 3] || shr_isSaturated[threadIdx + 4];
              }

              alpaka::syncBlockThreads(acc);

              bool check_hasSwitchToGain0 = false;

              if (sample == 0) {
                shr_hasSwitchToGain6[threadIdx] =
                    shr_hasSwitchToGain6[threadIdx] || shr_hasSwitchToGain6[threadIdx + 1];
                shr_hasSwitchToGain1[threadIdx] =
                    shr_hasSwitchToGain1[threadIdx] || shr_hasSwitchToGain1[threadIdx + 1];
                shr_hasSwitchToGain0_tmp[threadIdx] =
                    shr_hasSwitchToGain0_tmp[threadIdx] || shr_hasSwitchToGain0_tmp[threadIdx + 1];

                hasSwitchToGain6[ch] = shr_hasSwitchToGain6[threadIdx];
                hasSwitchToGain1[ch] = shr_hasSwitchToGain1[threadIdx];

                // set only for the threadIdx corresponding to sample==0
                check_hasSwitchToGain0 = shr_hasSwitchToGain0_tmp[threadIdx];

                shr_isSaturated[threadIdx + 3] = shr_isSaturated[threadIdx] || shr_isSaturated[threadIdx + 1];
                isSaturated[ch] = shr_isSaturated[threadIdx + 3];
              }

              // TODO: w/o this sync, there is a race
              // if (threadIdx == sample_max) below uses max sample thread, not for 0 sample
              // check if we can remove it

              alpaka::syncBlockThreads(acc);

              // TODO: divergent branch
              if (gainId == 0 || gainId == 3) {
                pedestal = conditionsDev.pedestals_mean_x1()[hashedId];
                gainratio = conditionsDev.gain6Over1()[hashedId] * conditionsDev.gain12Over6()[hashedId];
                gainsNoise[ch](sample) = 2;
              } else if (gainId == 1) {
                pedestal = conditionsDev.pedestals_mean_x12()[hashedId];
                gainratio = 1.;
                gainsNoise[ch](sample) = 0;
              } else if (gainId == 2) {
                pedestal = conditionsDev.pedestals_mean_x6()[hashedId];
                gainratio = conditionsDev.gain12Over6()[hashedId];
                gainsNoise[ch](sample) = 1;
              }

              // TODO: compile time constant -> branch should be non-divergent
              if (dynamicPedestal)
                amplitude = static_cast<::ecal::multifit::SampleVector::Scalar>(adc) * gainratio;
              else
                amplitude = (static_cast<::ecal::multifit::SampleVector::Scalar>(adc) - pedestal) * gainratio;
              amplitudes[ch][sample] = amplitude;

#ifdef ECAL_RECO_ALPAKA_DEBUG
              printf("%d %d %d %d %f %f %f\n", tx, ch, sample, adc, amplitude, pedestal, gainratio);
              if (adc == 0)
                printf("adc is zero\n");
#endif

              //
              // initialization
              //
              amplitudesForMinimization[inputCh](sample) = 0;
              bxs[ch](sample) = sample - 5;

              // select the thread for the max sample
              //---> hardcoded above to be 5th sample, ok
              if (sample == sample_max) {
                //
                // initialization
                //
                acState[ch] = static_cast<char>(MinimizationState::NotFinished);
                energies[inputCh] = 0;
                chi2[inputCh] = 0;
                g_pedestal[inputCh] = 0;
                uint32_t flag = 0;
                dids_out[inputCh] = did.rawId();

                // start of this channel in shared mem
                auto const chStart = threadIdx - sample_max;
                // thread for the max sample in shared mem
                auto const threadMax = threadIdx;
                auto const gainSwitchUseMaxSample = isBarrel ? gainSwitchUseMaxSampleEB : gainSwitchUseMaxSampleEE;

                // this flag setting is applied to all of the cases
                if (shr_hasSwitchToGain6[chStart])
                  flag |= 0x1 << EcalUncalibratedRecHit::kHasSwitchToGain6;
                if (shr_hasSwitchToGain1[chStart])
                  flag |= 0x1 << EcalUncalibratedRecHit::kHasSwitchToGain1;

                // this corresponds to cpu branching on lastSampleBeforeSaturation
                // likely false
                if (check_hasSwitchToGain0) {
                  // assign for the case some sample having gainId == 0
                  //energies[inputCh] = amplitudes[ch][sample_max];
                  energies[inputCh] = amplitude;

                  // check if samples before sample_max have true
                  bool saturated_before_max = false;
                  CMS_UNROLL_LOOP
                  for (char ii = 0; ii < 5; ++ii)
                    saturated_before_max = saturated_before_max || shr_hasSwitchToGain0[chStart + ii];

                  // if saturation is in the max sample and not in the first 5
                  if (!saturated_before_max && shr_hasSwitchToGain0[threadMax])
                    energies[inputCh] = 49140;  // 4095 * 12 (maximum ADC range * MultiGainPreAmplifier (MGPA) gain)
                                                // This is the actual maximum range that is set when we saturate.
                                                //---- AM FIXME : no pedestal subtraction???
                                                //It should be "(4095. - pedestal) * gainratio"

                  // set state flag to terminate further processing of this channel
                  acState[ch] = static_cast<char>(MinimizationState::Precomputed);
                  flag |= 0x1 << EcalUncalibratedRecHit::kSaturated;
                  flags[inputCh] = flag;
                  return;
                }

                // according to cpu version
                //            auto max_amplitude = amplitudes[ch][sample_max];
                auto const max_amplitude = amplitude;
                // according to cpu version
                auto shape_value = shapes_in[hashedId].pdfval[full_pulse_max - 7];
                // note, no syncing as the same thread will be accessing here
                bool hasGainSwitch =
                    shr_hasSwitchToGain6[chStart] || shr_hasSwitchToGain1[chStart] || shr_isSaturated[chStart + 3];

                // pedestal is final unconditionally
                g_pedestal[inputCh] = pedestal;
                if (hasGainSwitch && gainSwitchUseMaxSample) {
                  // thread for sample=0 will access the right guys
                  energies[inputCh] = max_amplitude / shape_value;
                  acState[ch] = static_cast<char>(MinimizationState::Precomputed);
                  flags[inputCh] = flag;
                  return;
                }

                // this happens cause sometimes rms_x12 is 0...
                // needs to be checkec why this is the case
                // general case here is that noisecov is a Zero matrix
                if (rmsForChecking == 0) {
                  acState[ch] = static_cast<char>(MinimizationState::Precomputed);
                  flags[inputCh] = flag;
                  return;
                }

                // for the case when no shortcuts were taken
                flags[inputCh] = flag;
              }
            }
          }
      };

      ///
      /// assume kernel launch configuration is
      /// ([MAXSAMPLES, MAXSAMPLES], nchannels)
      ///
      class kernel_prep_2d {
        public:
          template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
          ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                        DigiDeviceCollection::ConstView digisDevEB,
                                        DigiDeviceCollection::ConstView digisDevEE,
                                        EcalMultifitConditionsPortableDevice::ConstView conditionsDev,
                                        ::ecal::multifit::SampleGainVector const* gainsNoise,
                                        ::ecal::multifit::SampleMatrix* noisecov,
                                        ::ecal::multifit::PulseMatrixType* pulse_matrix,
                                        bool const* hasSwitchToGain6,
                                        bool const* hasSwitchToGain1,
                                        bool const* isSaturated) const {
            auto const offsetForHashes = conditionsDev.offsetEE();
            auto const nchannelsEB = digisDevEB.size();
            auto const ch = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];
            auto const tx = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u];
            auto const ty = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[1u];
            constexpr float addPedestalUncertainty = 0.f;
            constexpr bool dynamicPedestal = false;
            constexpr bool simplifiedNoiseModelForGainSwitch = true;  //---- default is true

            // to access input arrays (ids and digis only)
            int const inputCh = ch >= nchannelsEB ? ch - nchannelsEB : ch;
            auto const* dids = ch >= nchannelsEB ? digisDevEE.id() : digisDevEB.id();

            auto const did = DetId{dids[inputCh]};
            auto const isBarrel = did.subdetId() == EcalBarrel;
            auto const hashedId = isBarrel ? ecal::reconstruction::hashedIndexEB(did.rawId())
                                           : offsetForHashes + ecal::reconstruction::hashedIndexEE(did.rawId());
            auto const* G12SamplesCorrelation = isBarrel ? conditionsDev.sampleCorrelation_EB_G12().data() : conditionsDev.sampleCorrelation_EE_G12().data();
            auto const* G6SamplesCorrelation = isBarrel ? conditionsDev.sampleCorrelation_EB_G6().data() : conditionsDev.sampleCorrelation_EE_G6().data();
            auto const* G1SamplesCorrelation = isBarrel ? conditionsDev.sampleCorrelation_EB_G1().data() : conditionsDev.sampleCorrelation_EE_G1().data();
            auto const hasGainSwitch = hasSwitchToGain6[ch] || hasSwitchToGain1[ch] || isSaturated[ch];
            auto const vidx = std::abs(static_cast<int>(ty) - static_cast<int>(tx));

            // non-divergent branch for all threads per block
            if (hasGainSwitch) {
              // TODO: did not include simplified noise model
              float noise_value = 0;

              // non-divergent branch - all threads per block
              // TODO: all of these constants indicate that
              // that these parts could be splitted into completely different
              // kernels and run one of them only depending on the config
              if (simplifiedNoiseModelForGainSwitch) {
                constexpr int isample_max = 5;  // according to cpu defs
                auto const gainidx = gainsNoise[ch][isample_max];

                // non-divergent branches
                if (gainidx == 0) {
                  auto const rms_x12 = conditionsDev.pedestals_rms_x12()[hashedId];
                  noise_value = rms_x12 * rms_x12 * G12SamplesCorrelation[vidx];
                } else if (gainidx == 1) {
                  auto const gain12Over6 = conditionsDev.gain12Over6()[hashedId];
                  auto const rms_x6 = conditionsDev.pedestals_rms_x6()[hashedId];
                  noise_value = gain12Over6 * gain12Over6 * rms_x6 * rms_x6 * G6SamplesCorrelation[vidx];
                } else if (gainidx == 2) {
                  auto const gain12Over6 = conditionsDev.gain12Over6()[hashedId];
                  auto const gain6Over1 = conditionsDev.gain6Over1()[hashedId];
                  auto const gain12Over1 = gain12Over6 * gain6Over1;
                  auto const rms_x1 = conditionsDev.pedestals_rms_x1()[hashedId];
                  noise_value = gain12Over1 * gain12Over1 * rms_x1 * rms_x1 * G1SamplesCorrelation[vidx];
                }
                if (!dynamicPedestal && addPedestalUncertainty > 0.f)
                  noise_value += addPedestalUncertainty * addPedestalUncertainty;
              } else {
                int gainidx = 0;
                char mask = gainidx;
                int pedestal = gainsNoise[ch][ty] == mask ? 1 : 0;
                //            NB: gainratio is 1, that is why it does not appear in the formula
                auto const rms_x12 = conditionsDev.pedestals_rms_x12()[hashedId];
                noise_value += rms_x12 * rms_x12 * pedestal * G12SamplesCorrelation[vidx];
                // non-divergent branch
                if (!dynamicPedestal && addPedestalUncertainty > 0.f) {
                  noise_value += addPedestalUncertainty * addPedestalUncertainty * pedestal;  // gainratio is 1
                }

                //
                gainidx = 1;
                mask = gainidx;
                pedestal = gainsNoise[ch][ty] == mask ? 1 : 0;
                auto const gain12Over6 = conditionsDev.gain12Over6()[hashedId];
                auto const rms_x6 = conditionsDev.pedestals_rms_x6()[hashedId];
                noise_value += gain12Over6 * gain12Over6 * rms_x6 * rms_x6 * pedestal * G6SamplesCorrelation[vidx];
                // non-divergent branch
                if (!dynamicPedestal && addPedestalUncertainty > 0.f) {
                  noise_value += gain12Over6 * gain12Over6 * addPedestalUncertainty * addPedestalUncertainty * pedestal;
                }

                //
                gainidx = 2;
                mask = gainidx;
                pedestal = gainsNoise[ch][ty] == mask ? 1 : 0;
                auto const gain6Over1 = conditionsDev.gain6Over1()[hashedId];
                auto const gain12Over1 = gain12Over6 * gain6Over1;
                auto const rms_x1 = conditionsDev.pedestals_rms_x1()[hashedId];
                noise_value += gain12Over1 * gain12Over1 * rms_x1 * rms_x1 * pedestal * G1SamplesCorrelation[vidx];
                // non-divergent branch
                if (!dynamicPedestal && addPedestalUncertainty > 0.f) {
                  noise_value += gain12Over1 * gain12Over1 * addPedestalUncertainty * addPedestalUncertainty * pedestal;
                }
              }

              noisecov[ch](ty, tx) = noise_value;
            } else {
              auto const rms = conditionsDev.pedestals_rms_x12()[hashedId];
              float noise_value = rms * rms * G12SamplesCorrelation[vidx];
              if (!dynamicPedestal && addPedestalUncertainty > 0.f) {
                //----  add fully correlated component to noise covariance to inflate pedestal uncertainty
                noise_value += addPedestalUncertainty * addPedestalUncertainty;
              }
              noisecov[ch](ty, tx) = noise_value;
            }

            // pulse matrix
            auto const* pulse_shapes = reinterpret_cast<const EcalPulseShape*>(conditionsDev.pulseShapes()->data());
            auto const posToAccess = 9 - static_cast<int>(tx + ty);  // see cpu for reference
            float const value = posToAccess >= 7 ? pulse_shapes[hashedId].pdfval[posToAccess - 7] : 0;
            pulse_matrix[ch](ty, tx) = value;
          }
      };

    }  // namespace multifit
  }  // namespace ecal
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

namespace alpaka::trait
{
  using namespace ALPAKA_ACCELERATOR_NAMESPACE::ecal::multifit;

  //! The trait for getting the size of the block shared dynamic memory for kernel_prep_1d_and_initialize.
  template<typename TAcc>
  struct BlockSharedMemDynSizeBytes<kernel_prep_1d_and_initialize, TAcc>
  {
    //! \return The size of the shared memory allocated for a block.
    template<typename TVec, typename... TArgs>
    ALPAKA_FN_HOST_ACC static auto getBlockSharedMemDynSizeBytes(
      kernel_prep_1d_and_initialize const&,
      TVec const& threadsPerBlock,
      TVec const&,
      TArgs const&...) -> std::size_t
    {
      // return the amount of dynamic shared memory needed
      std::size_t bytes = threadsPerBlock[0u] * EcalDataFrame::MAXSAMPLES * (5 * sizeof(bool) + sizeof(char));
      return bytes;
    }
  };
} // namespace alpaka::trait

#endif  // RecoLocalCalo_EcalRecProducers_plugins_AmplitudeComputationCommonKernels_h
