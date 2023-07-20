// Check that ALPAKA_HOST_ONLY is not defined during device compilation:
#ifdef ALPAKA_HOST_ONLY
#error ALPAKA_HOST_ONLY defined in device compilation
#endif

#include <limits>
#include <alpaka/alpaka.hpp>

//#include "CondFormats/EcalObjects/interface/EcalMGPAGainRatio.h"
//#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
//#include "CondFormats/EcalObjects/interface/EcalPulseCovariances.h"
//#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"
//#include "CondFormats/EcalObjects/interface/EcalSampleMask.h"
//#include "CondFormats/EcalObjects/interface/EcalSamplesCorrelation.h"
//#include "CondFormats/EcalObjects/interface/EcalXtalGroupId.h"
//#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
//#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "AmplitudeComputationCommonKernels.h"
#include "AmplitudeComputationKernels.h"
#include "EcalUncalibRecHitMultiFitAlgoPortable.h"
#include "TimeComputationKernels.h"

//#define DEBUG
//#define ECAL_RECO_CUDA_DEBUG

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace ecal { 
    namespace multifit {
  
      using namespace cms::alpakatools;

      void entryPoint(InputProduct const& digisDevEB,
                      InputProduct const& digisDevEE,
                      OutputProduct& uncalibRecHitsDevEB,
                      OutputProduct& uncalibRecHitsDevEE,
                      EcalMultifitConditionsPortableDevice const& conditionsDev,
                      ConfigurationParameters const& configParams,
                      Queue& queue) {
  
      //void entryPoint(EventInputDataGPU const& eventInputGPU,
      //                EventOutputDataGPU& eventOutputGPU,
      //                EventDataForScratchGPU& scratch,
      //                ConditionsProducts const& conditions,
      //                ConfigurationParameters const& configParameters,
      //                cudaStream_t cudaStream) {
        using digis_type = std::vector<uint16_t>;
        using dids_type = std::vector<uint32_t>;
        // accodring to the cpu setup  //----> hardcoded
        bool const gainSwitchUseMaxSampleEB = true;
        // accodring to the cpu setup  //----> hardcoded
        bool const gainSwitchUseMaxSampleEE = false;
  
        auto const totalChannels = static_cast<uint32_t>(digisDevEB->metadata().size() + digisDevEE->metadata().size());

        //
        // 1d preparation kernel
        //
        uint32_t const nchannels_per_block = 32;
        auto const threads_1d = EcalDataFrame::MAXSAMPLES * nchannels_per_block;
        auto const blocks_1d = threads_1d > EcalDataFrame::MAXSAMPLES * totalChannels ? 1u : (totalChannels * EcalDataFrame::MAXSAMPLES + threads_1d - 1) / threads_1d;
      //  int shared_bytes = nchannels_per_block * EcalDataFrame::MAXSAMPLES *
      //                     (sizeof(bool) + sizeof(bool) + sizeof(bool) + sizeof(bool) + sizeof(char) + sizeof(bool));
        auto workDivPrep1D = cms::alpakatools::make_workdiv<Acc1D>(blocks_1d, threads_1d);
        alpaka::exec<Acc1D>(queue,
                            workDivPrep1D,
                            kernel_prep_1d_and_initialize{},
                            digisDevEB.const_view(),
                            digisDevEE.const_view(),
                            uncalibRecHitsDevEB.view(),
                            uncalibRecHitsDevEE.view(),
                            conditionsDev.const_view());
  
        //
        // 2d preparation kernel
        //
        Vec2D const blocks_2d{totalChannels, 1u};
        Vec2D const threads_2d{EcalDataFrame::MAXSAMPLES, EcalDataFrame::MAXSAMPLES};
        auto workDivPrep2D = cms::alpakatools::make_workdiv<Acc2D>(blocks_2d, threads_2d);
        alpaka::exec<Acc2D>(queue,
                            workDivPrep2D,
                            kernel_prep_2d{},
                            digisDevEB.const_view(),
                            digisDevEE.const_view());
      //  kernel_prep_2d<<<blocks_2d, threads_2d, 0, cudaStream>>>((SampleGainVector*)scratch.gainsNoise.get(),
      //                                                           eventInputGPU.ebDigis.ids.get(),
      //                                                           eventInputGPU.eeDigis.ids.get(),
      //                                                           conditions.pedestals.rms_x12,
      //                                                           conditions.pedestals.rms_x6,
      //                                                           conditions.pedestals.rms_x1,
      //                                                           conditions.gainRatios.gain12Over6,
      //                                                           conditions.gainRatios.gain6Over1,
      //                                                           conditions.samplesCorrelation.EBG12SamplesCorrelation,
      //                                                           conditions.samplesCorrelation.EBG6SamplesCorrelation,
      //                                                           conditions.samplesCorrelation.EBG1SamplesCorrelation,
      //                                                           conditions.samplesCorrelation.EEG12SamplesCorrelation,
      //                                                           conditions.samplesCorrelation.EEG6SamplesCorrelation,
      //                                                           conditions.samplesCorrelation.EEG1SamplesCorrelation,
      //                                                           (SampleMatrix*)scratch.noisecov.get(),
      //                                                           (PulseMatrixType*)scratch.pulse_matrix.get(),
      //                                                           conditions.pulseShapes.values,
      //                                                           scratch.hasSwitchToGain6.get(),
      //                                                           scratch.hasSwitchToGain1.get(),
      //                                                           scratch.isSaturated.get(),
      //                                                           offsetForHashes,
      //                                                           offsetForInputs);
      //  cudaCheck(cudaGetLastError());
  
      //  // run minimization kernels
      //  v1::minimization_procedure(eventInputGPU, eventOutputGPU, scratch, conditions, configParameters, cudaStream);
  
        if (configParams.shouldRunTimingComputation) {
          //
          // TODO: this guy can run concurrently with other kernels,
          // there is no dependence on the order of execution
          //
          auto const blocks_time_init = blocks_1d;
          auto const threads_time_init = threads_1d;
          auto workDivTimeCompInit1D = cms::alpakatools::make_workdiv<Acc1D>(blocks_time_init, threads_time_init);
          alpaka::exec<Acc1D>(queue,
                              workDivTimeCompInit1D,
                              kernel_time_computation_init{},
                              digisDevEB.const_view(),
                              digisDevEE.const_view());
           //int sharedBytesInit = 2 * threads_time_init * sizeof(SampleVector::Scalar);
      //    kernel_time_computation_init<<<blocks_time_init, threads_time_init, sharedBytesInit, cudaStream>>>(
      //        eventInputGPU.ebDigis.data.get(),
      //        eventInputGPU.ebDigis.ids.get(),
      //        eventInputGPU.eeDigis.data.get(),
      //        eventInputGPU.eeDigis.ids.get(),
      //        conditions.pedestals.rms_x12,
      //        conditions.pedestals.rms_x6,
      //        conditions.pedestals.rms_x1,
      //        conditions.pedestals.mean_x12,
      //        conditions.pedestals.mean_x6,
      //        conditions.pedestals.mean_x1,
      //        conditions.gainRatios.gain12Over6,
      //        conditions.gainRatios.gain6Over1,
      //        scratch.sample_values.get(),
      //        scratch.sample_value_errors.get(),
      //        scratch.ampMaxError.get(),
      //        scratch.useless_sample_values.get(),
      //        scratch.pedestal_nums.get(),
      //        offsetForHashes,
      //        offsetForInputs,
      //        conditions.sampleMask.getEcalSampleMaskRecordEB(),
      //        conditions.sampleMask.getEcalSampleMaskRecordEE(),
      //        totalChannels);
      //    cudaCheck(cudaGetLastError());
  
          //
          // TODO: small kernel only for EB. It needs to be checked if
          /// fusing such small kernels is beneficial in here
          //
          // we are running only over EB digis
          // therefore we need to create threads/blocks only for that
          auto const threadsFixMGPA = threads_1d;
          auto const blocksFixMGPA =
              threadsFixMGPA > EcalDataFrame::MAXSAMPLES * static_cast<unsigned int>(digisDevEB->metadata().size()) // TODO check if metadata.size is OK here
                  ? 1
                  : (EcalDataFrame::MAXSAMPLES * static_cast<unsigned int>(digisDevEB->metadata().size()) + threadsFixMGPA - 1) / threadsFixMGPA;
          auto workDivTimeFixMGPAslew1D = cms::alpakatools::make_workdiv<Acc1D>(blocksFixMGPA, threadsFixMGPA);
          alpaka::exec<Acc1D>(queue,
                              workDivTimeFixMGPAslew1D,
                              kernel_time_compute_fixMGPAslew{},
                              digisDevEB.const_view(),
                              digisDevEE.const_view());
      //    kernel_time_compute_fixMGPAslew<<<blocksFixMGPA, threadsFixMGPA, 0, cudaStream>>>(
      //        eventInputGPU.ebDigis.data.get(),
      //        eventInputGPU.eeDigis.data.get(),
      //        scratch.sample_values.get(),
      //        scratch.sample_value_errors.get(),
      //        scratch.useless_sample_values.get(),
      //        conditions.sampleMask.getEcalSampleMaskRecordEB(),
      //        totalChannels,
      //        offsetForInputs);
      //    cudaCheck(cudaGetLastError());
  
      //    int sharedBytes = EcalDataFrame::MAXSAMPLES * nchannels_per_block * 4 * sizeof(SampleVector::Scalar);
          auto const threads_nullhypot = threads_1d;
          auto const blocks_nullhypot = blocks_1d;
          auto workDivTimeNullhypot1D = cms::alpakatools::make_workdiv<Acc1D>(blocks_nullhypot, threads_nullhypot);
          alpaka::exec<Acc1D>(queue,
                              workDivTimeNullhypot1D,
                              kernel_time_compute_nullhypot{});
      //    kernel_time_compute_nullhypot<<<blocks_nullhypot, threads_nullhypot, sharedBytes, cudaStream>>>(
      //        scratch.sample_values.get(),
      //        scratch.sample_value_errors.get(),
      //        scratch.useless_sample_values.get(),
      //        scratch.chi2sNullHypot.get(),
      //        scratch.sum0sNullHypot.get(),
      //        scratch.sumAAsNullHypot.get(),
      //        totalChannels);
      //    cudaCheck(cudaGetLastError());
  
          unsigned int const nchannels_per_block_makeratio = 10;
          auto const threads_makeratio = 45 * nchannels_per_block_makeratio;
          unsigned int const blocks_makeratio = threads_makeratio > 45 * totalChannels
                                              ? 1
                                              : (totalChannels * 45 + threads_makeratio - 1) / threads_makeratio;
      //    int sharedBytesMakeRatio = 5 * threads_makeratio * sizeof(SampleVector::Scalar);
          auto workDivTimeMakeRatio1D = cms::alpakatools::make_workdiv<Acc1D>(blocks_makeratio, threads_makeratio);
          alpaka::exec<Acc1D>(queue,
                              workDivTimeMakeRatio1D,
                              kernel_time_compute_makeratio{},
                              digisDevEB.const_view(),
                              digisDevEE.const_view());
      //    kernel_time_compute_makeratio<<<blocks_makeratio, threads_makeratio, sharedBytesMakeRatio, cudaStream>>>(
      //        scratch.sample_values.get(),
      //        scratch.sample_value_errors.get(),
      //        eventInputGPU.ebDigis.ids.get(),
      //        eventInputGPU.eeDigis.ids.get(),
      //        scratch.useless_sample_values.get(),
      //        scratch.pedestal_nums.get(),
      //        configParameters.amplitudeFitParametersEB,
      //        configParameters.amplitudeFitParametersEE,
      //        configParameters.timeFitParametersEB,
      //        configParameters.timeFitParametersEE,
      //        scratch.sumAAsNullHypot.get(),
      //        scratch.sum0sNullHypot.get(),
      //        scratch.tMaxAlphaBetas.get(),
      //        scratch.tMaxErrorAlphaBetas.get(),
      //        scratch.accTimeMax.get(),
      //        scratch.accTimeWgt.get(),
      //        scratch.tcState.get(),
      //        configParameters.timeFitParametersSizeEB,
      //        configParameters.timeFitParametersSizeEE,
      //        configParameters.timeFitLimitsFirstEB,
      //        configParameters.timeFitLimitsFirstEE,
      //        configParameters.timeFitLimitsSecondEB,
      //        configParameters.timeFitLimitsSecondEE,
      //        totalChannels,
      //        offsetForInputs);
      //    cudaCheck(cudaGetLastError());
  
          auto const threads_findamplchi2 = threads_1d;
          auto const blocks_findamplchi2 = blocks_1d;
      //    int const sharedBytesFindAmplChi2 = 2 * threads_findamplchi2 * sizeof(SampleVector::Scalar);
          auto workDivTimeFindAmplChi21D = cms::alpakatools::make_workdiv<Acc1D>(blocks_findamplchi2, threads_findamplchi2);
          alpaka::exec<Acc1D>(queue,
                              workDivTimeFindAmplChi21D,
                              kernel_time_compute_findamplchi2_and_finish{},
                              digisDevEB.const_view(),
                              digisDevEE.const_view());
      //    kernel_time_compute_findamplchi2_and_finish<<<blocks_findamplchi2,
      //                                                  threads_findamplchi2,
      //                                                  sharedBytesFindAmplChi2,
      //                                                  cudaStream>>>(scratch.sample_values.get(),
      //                                                                scratch.sample_value_errors.get(),
      //                                                                eventInputGPU.ebDigis.ids.get(),
      //                                                                eventInputGPU.eeDigis.ids.get(),
      //                                                                scratch.useless_sample_values.get(),
      //                                                                scratch.tMaxAlphaBetas.get(),
      //                                                                scratch.tMaxErrorAlphaBetas.get(),
      //                                                                scratch.accTimeMax.get(),
      //                                                                scratch.accTimeWgt.get(),
      //                                                                configParameters.amplitudeFitParametersEB,
      //                                                                configParameters.amplitudeFitParametersEE,
      //                                                                scratch.sumAAsNullHypot.get(),
      //                                                                scratch.sum0sNullHypot.get(),
      //                                                                scratch.chi2sNullHypot.get(),
      //                                                                scratch.tcState.get(),
      //                                                                scratch.ampMaxAlphaBeta.get(),
      //                                                                scratch.ampMaxError.get(),
      //                                                                scratch.timeMax.get(),
      //                                                                scratch.timeError.get(),
      //                                                                totalChannels,
      //                                                                offsetForInputs);
      //    cudaCheck(cudaGetLastError());
  
          auto const threads_timecorr = 32;
          auto const blocks_timecorr =
              threads_timecorr > totalChannels ? 1 : (totalChannels + threads_timecorr - 1) / threads_timecorr;
          auto workDivCorrFinal1D = cms::alpakatools::make_workdiv<Acc1D>(blocks_timecorr, threads_timecorr);
          alpaka::exec<Acc1D>(queue,
                              workDivCorrFinal1D,
                              kernel_time_correction_and_finalize{},
                              digisDevEB.const_view(),
                              digisDevEE.const_view(),
                              uncalibRecHitsDevEB.view(),
                              uncalibRecHitsDevEE.view());
      //    kernel_time_correction_and_finalize<<<blocks_timecorr, threads_timecorr, 0, cudaStream>>>(
      //        eventOutputGPU.recHitsEB.amplitude.get(),
      //        eventOutputGPU.recHitsEE.amplitude.get(),
      //        eventInputGPU.ebDigis.data.get(),
      //        eventInputGPU.ebDigis.ids.get(),
      //        eventInputGPU.eeDigis.data.get(),
      //        eventInputGPU.eeDigis.ids.get(),
      //        conditions.timeBiasCorrections.ebTimeCorrAmplitudeBins,
      //        conditions.timeBiasCorrections.eeTimeCorrAmplitudeBins,
      //        conditions.timeBiasCorrections.ebTimeCorrShiftBins,
      //        conditions.timeBiasCorrections.eeTimeCorrShiftBins,
      //        scratch.timeMax.get(),
      //        scratch.timeError.get(),
      //        conditions.pedestals.rms_x12,
      //        conditions.timeCalibConstants.values,
      //        eventOutputGPU.recHitsEB.jitter.get(),
      //        eventOutputGPU.recHitsEE.jitter.get(),
      //        eventOutputGPU.recHitsEB.jitterError.get(),
      //        eventOutputGPU.recHitsEE.jitterError.get(),
      //        eventOutputGPU.recHitsEB.flags.get(),
      //        eventOutputGPU.recHitsEE.flags.get(),
      //        conditions.timeBiasCorrections.ebTimeCorrAmplitudeBinsSize,
      //        conditions.timeBiasCorrections.eeTimeCorrAmplitudeBinsSize,
      //        configParameters.timeConstantTermEB,
      //        configParameters.timeConstantTermEE,
      //        conditions.timeOffsetConstant.getEBValue(),
      //        conditions.timeOffsetConstant.getEEValue(),
      //        configParameters.timeNconstEB,
      //        configParameters.timeNconstEE,
      //        configParameters.amplitudeThreshEB,
      //        configParameters.amplitudeThreshEE,
      //        configParameters.outOfTimeThreshG12pEB,
      //        configParameters.outOfTimeThreshG12pEE,
      //        configParameters.outOfTimeThreshG12mEB,
      //        configParameters.outOfTimeThreshG12mEE,
      //        configParameters.outOfTimeThreshG61pEB,
      //        configParameters.outOfTimeThreshG61pEE,
      //        configParameters.outOfTimeThreshG61mEB,
      //        configParameters.outOfTimeThreshG61mEE,
      //        offsetForHashes,
      //        offsetForInputs,
      //        totalChannels);
      //    cudaCheck(cudaGetLastError());
        }
      }
  
    }  // namespace multifit
  }  // namespace ecal
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
