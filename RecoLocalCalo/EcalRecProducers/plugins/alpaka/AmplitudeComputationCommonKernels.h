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

// this flag setting is applied to all of the cases
struct EcalPulseCovariance;
  
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
            constexpr int nsamples = EcalDataFrame::MAXSAMPLES;
            constexpr int sample_max = 5;
            constexpr int full_pulse_max = 9;

            auto const nchannelsEB = digisDevEB.size();
            auto const nchannels = nchannelsEB + digisDevEE.size();
            auto const threadIdx = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u];
            auto const blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];
            auto const blockDim = alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u];
            auto const tx = threadIdx + blockIdx * blockDim;
            //auto const nchannels_per_block = blockDim / nsamples;
            uint32_t constexpr nchannels_per_block = 32;
            auto const ch = tx / nsamples;
            // for accessing input arrays
            int const inputCh = ch >= nchannelsEB ? ch - nchannelsEB : ch;
            int const inputTx = ch >= nchannelsEB ? tx - nchannelsEB * nsamples : tx;
            // eb is first and then ee
            auto const* digis_in = ch >= nchannelsEB ? digisDevEE.data()->array.data() : digisDevEB.data()->array.data(); // is there a way to use digisDevEE.data()->data() instead?
            auto const* dids = ch >= nchannelsEB ? digisDevEE.id() : digisDevEB.id();
            auto const sample = threadIdx % nsamples;

            auto* amplitudesForMinimization = (::ecal::multifit::SampleVector*)(ch >= nchannelsEB ? uncalibRecHitsEE.outOfTimeAmplitudes()->array.data() : uncalibRecHitsEB.outOfTimeAmplitudes()->array.data());
            auto* energies = ch >= nchannelsEB ? uncalibRecHitsEE.amplitude() : uncalibRecHitsEB.amplitude();
            auto* chi2 = ch >= nchannelsEB ? uncalibRecHitsEE.chi2() : uncalibRecHitsEB.chi2();
            auto* g_pedestal = ch >= nchannelsEB ? uncalibRecHitsEE.pedestal() : uncalibRecHitsEB.pedestal();
            auto* dids_out = ch >= nchannelsEB ? uncalibRecHitsEE.id() : uncalibRecHitsEB.id();
            auto* flags = ch >= nchannelsEB ? uncalibRecHitsEE.flags() : uncalibRecHitsEB.flags();

            auto const* shapes_in = (EcalPulseShape*)conditionsDev.pulseShapes()->data();

            if (ch < nchannels) {
              auto& shr_hasSwitchToGain6 = alpaka::declareSharedVar<bool[nchannels_per_block * nsamples], __COUNTER__>(acc);
              auto& shr_hasSwitchToGain1 = alpaka::declareSharedVar<bool[nchannels_per_block * nsamples], __COUNTER__>(acc);
              auto& shr_hasSwitchToGain0 = alpaka::declareSharedVar<bool[nchannels_per_block * nsamples], __COUNTER__>(acc);
              auto& shr_isSaturated = alpaka::declareSharedVar<bool[nchannels_per_block * nsamples], __COUNTER__>(acc);
              auto& shr_hasSwitchToGain0_tmp = alpaka::declareSharedVar<bool[nchannels_per_block * nsamples], __COUNTER__>(acc);
              auto& shr_counts = alpaka::declareSharedVar<uint32_t[nchannels_per_block * nsamples], __COUNTER__>(acc);

              //
              // indices
              //
              auto const did = DetId{dids[inputCh]};
              auto const isBarrel = did.subdetId() == EcalBarrel;
              // TODO offset for ee, 0 for eb
              auto const hashedId = isBarrel ? reconstruction::hashedIndexEB(did.rawId())
                                             : nchannelsEB + reconstruction::hashedIndexEE(did.rawId());

              //
              // pulse shape template

              // will be used in the future for setting state
              auto const rmsForChecking = conditionsDev.pedestals_rms_x12()[hashedId];

              //
              // amplitudes
              //
              int const adc = ecalMGPA::adc(digis_in[inputTx]);
              int const gainId = ecalMGPA::gainId(digis_in[inputTx]);
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
                for (int i = 0; i < 5; i++)
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
                int const chStart = threadIdx - sample_max;
                // thread for the max sample in shared mem
                int const threadMax = threadIdx;
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
                  for (char ii = 0; ii < 5; ii++)
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
                                        DigiDeviceCollection::ConstView digisDevEE) const {

      //__global__ void kernel_prep_2d(SampleGainVector const* gainNoise,
      //                               uint32_t const* dids_eb,
      //                               uint32_t const* dids_ee,
      //                               float const* rms_x12,
      //                               float const* rms_x6,
      //                               float const* rms_x1,
      //                               float const* gain12Over6,
      //                               float const* gain6Over1,
      //                               double const* G12SamplesCorrelationEB,
      //                               double const* G6SamplesCorrelationEB,
      //                               double const* G1SamplesCorrelationEB,
      //                               double const* G12SamplesCorrelationEE,
      //                               double const* G6SamplesCorrelationEE,
      //                               double const* G1SamplesCorrelationEE,
      //                               SampleMatrix* noisecov,
      //                               PulseMatrixType* pulse_matrix,
      //                               EcalPulseShape const* pulse_shape,
      //                               bool const* hasSwitchToGain6,
      //                               bool const* hasSwitchToGain1,
      //                               bool const* isSaturated,
      //                               uint32_t const offsetForHashes,
      //                               uint32_t const offsetForInputs);
      //      int const ch = blockIdx.x;
      //      int const tx = threadIdx.x;
      //      int const ty = threadIdx.y;
      //      constexpr float addPedestalUncertainty = 0.f;
      //      constexpr bool dynamicPedestal = false;
      //      constexpr bool simplifiedNoiseModelForGainSwitch = true;  //---- default is true

      //      // to access input arrays (ids and digis only)
      //      int const inputCh = ch >= offsetForInputs ? ch - offsetForInputs : ch;
      //      auto const* dids = ch >= offsetForInputs ? dids_ee : dids_eb;

      //      bool tmp0 = hasSwitchToGain6[ch];
      //      bool tmp1 = hasSwitchToGain1[ch];
      //      auto const did = DetId{dids[inputCh]};
      //      auto const isBarrel = did.subdetId() == EcalBarrel;
      //      auto const hashedId = isBarrel ? ecal::reconstruction::hashedIndexEB(did.rawId())
      //                                     : offsetForHashes + ecal::reconstruction::hashedIndexEE(did.rawId());
      //      auto const G12SamplesCorrelation = isBarrel ? G12SamplesCorrelationEB : G12SamplesCorrelationEE;
      //      auto const* G6SamplesCorrelation = isBarrel ? G6SamplesCorrelationEB : G6SamplesCorrelationEE;
      //      auto const* G1SamplesCorrelation = isBarrel ? G1SamplesCorrelationEB : G1SamplesCorrelationEE;
      //      bool tmp2 = isSaturated[ch];
      //      bool hasGainSwitch = tmp0 || tmp1 || tmp2;
      //      auto const vidx = std::abs(ty - tx);

      //      // non-divergent branch for all threads per block
      //      if (hasGainSwitch) {
      //        // TODO: did not include simplified noise model
      //        float noise_value = 0;

      //        // non-divergent branch - all threads per block
      //        // TODO: all of these constants indicate that
      //        // that these parts could be splitted into completely different
      //        // kernels and run one of them only depending on the config
      //        if (simplifiedNoiseModelForGainSwitch) {
      //          int isample_max = 5;  // according to cpu defs
      //          int gainidx = gainNoise[ch][isample_max];

      //          // non-divergent branches
      //          if (gainidx == 0)
      //            noise_value = rms_x12[hashedId] * rms_x12[hashedId] * G12SamplesCorrelation[vidx];
      //          if (gainidx == 1)
      //            noise_value = gain12Over6[hashedId] * gain12Over6[hashedId] * rms_x6[hashedId] * rms_x6[hashedId] *
      //                          G6SamplesCorrelation[vidx];
      //          if (gainidx == 2)
      //            noise_value = gain12Over6[hashedId] * gain12Over6[hashedId] * gain6Over1[hashedId] * gain6Over1[hashedId] *
      //                          rms_x1[hashedId] * rms_x1[hashedId] * G1SamplesCorrelation[vidx];
      //          if (!dynamicPedestal && addPedestalUncertainty > 0.f)
      //            noise_value += addPedestalUncertainty * addPedestalUncertainty;
      //        } else {
      //          int gainidx = 0;
      //          char mask = gainidx;
      //          int pedestal = gainNoise[ch][ty] == mask ? 1 : 0;
      //          //            NB: gainratio is 1, that is why it does not appear in the formula
      //          noise_value += rms_x12[hashedId] * rms_x12[hashedId] * pedestal * G12SamplesCorrelation[vidx];
      //          // non-divergent branch
      //          if (!dynamicPedestal && addPedestalUncertainty > 0.f) {
      //            noise_value += addPedestalUncertainty * addPedestalUncertainty * pedestal;  // gainratio is 1
      //          }

      //          //
      //          gainidx = 1;
      //          mask = gainidx;
      //          pedestal = gainNoise[ch][ty] == mask ? 1 : 0;
      //          noise_value += gain12Over6[hashedId] * gain12Over6[hashedId] * rms_x6[hashedId] * rms_x6[hashedId] *
      //                         pedestal * G6SamplesCorrelation[vidx];
      //          // non-divergent branch
      //          if (!dynamicPedestal && addPedestalUncertainty > 0.f) {
      //            noise_value += gain12Over6[hashedId] * gain12Over6[hashedId] * addPedestalUncertainty *
      //                           addPedestalUncertainty * pedestal;
      //          }

      //          //
      //          gainidx = 2;
      //          mask = gainidx;
      //          pedestal = gainNoise[ch][ty] == mask ? 1 : 0;
      //          float tmp = gain6Over1[hashedId] * gain12Over6[hashedId];
      //          noise_value += tmp * tmp * rms_x1[hashedId] * rms_x1[hashedId] * pedestal * G1SamplesCorrelation[vidx];
      //          // non-divergent branch
      //          if (!dynamicPedestal && addPedestalUncertainty > 0.f) {
      //            noise_value += tmp * tmp * addPedestalUncertainty * addPedestalUncertainty * pedestal;
      //          }
      //        }

      //        noisecov[ch](ty, tx) = noise_value;
      //      } else {
      //        auto rms = rms_x12[hashedId];
      //        float noise_value = rms * rms * G12SamplesCorrelation[vidx];
      //        if (!dynamicPedestal && addPedestalUncertainty > 0.f) {
      //          //----  add fully correlated component to noise covariance to inflate pedestal uncertainty
      //          noise_value += addPedestalUncertainty * addPedestalUncertainty;
      //        }
      //        noisecov[ch](ty, tx) = noise_value;
      //      }

      //      // pulse matrix
      //      int const posToAccess = 9 - tx + ty;  // see cpu for reference
      //      float const value = posToAccess >= 7 ? pulse_shape[hashedId].pdfval[posToAccess - 7] : 0;
      //      pulse_matrix[ch](ty, tx) = value;
          }
      };
 
      class kernel_permute_results {
        public:
          template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
          ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                        DigiDeviceCollection::ConstView digisDevEB,
                                        DigiDeviceCollection::ConstView digisDevEE,
                                        int const nchannels) const {

      //__global__ void kernel_permute_results(SampleVector* amplitudes,
      //                                       BXVectorType const* activeBXs,
      //                                       ::ecal::reco::StorageScalarType* energies,
      //                                       char const* acState,
      //                                       int const nchannels);
      //      // constants
      //      constexpr int nsamples = EcalDataFrame::MAXSAMPLES;

      //      // indices
      //      int const tx = threadIdx.x + blockIdx.x * blockDim.x;
      //      int const ch = tx / nsamples;
      //      int const sampleidx = tx % nsamples;  // this is to address activeBXs

      //      if (ch >= nchannels)
      //        return;

      //      // channels that have amplitude precomputed do not need results to be permuted
      //      auto const state = static_cast<MinimizationState>(acState[ch]);
      //      if (state == MinimizationState::Precomputed)
      //        return;

      //      // configure shared memory and cp into it
      //      extern __shared__ char smem[];
      //      SampleVector::Scalar* values = reinterpret_cast<SampleVector::Scalar*>(smem);
      //      values[threadIdx.x] = amplitudes[ch](sampleidx);
      //      __syncthreads();

      //      // get the sample for this bx
      //      auto const sample = static_cast<int>(activeBXs[ch](sampleidx)) + 5;

      //      // store back to global
      //      amplitudes[ch](sample) = values[threadIdx.x];

      //      // store sample 5 separately
      //      // only for the case when minimization was performed
      //      // not for cases with precomputed amplitudes
      //      if (sample == 5)
      //        energies[ch] = values[threadIdx.x];
          }
      };  
 
  ///
  /// Build an Ecal RecHit.
  /// TODO: Use SoA data structures on the host directly
  /// the reason for removing this from minimize kernel is to isolate the minimize +
  /// again, building an aos rec hit involves strides... -> bad memory access pattern
  ///
//  #ifdef RUN_BUILD_AOS_RECHIT
//      __global__ void kernel_build_rechit(
//          float const* energies, float const* chi2s, uint32_t* dids, EcalUncalibratedRecHit* rechits, int nchannels);
//  #endif  // RUN_BUILD_AOS_RECHIT
  
    }  // namespace multifit
  }  // namespace ecal
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoLocalCalo_EcalRecProducers_plugins_AmplitudeComputationCommonKernels_h
