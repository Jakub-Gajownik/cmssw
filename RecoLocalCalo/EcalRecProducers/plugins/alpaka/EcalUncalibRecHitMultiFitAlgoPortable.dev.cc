// Check that ALPAKA_HOST_ONLY is not defined during device compilation:
#ifdef ALPAKA_HOST_ONLY
#error ALPAKA_HOST_ONLY defined in device compilation
#endif

#include <limits>
#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "AmplitudeComputationCommonKernels.h"
#include "AmplitudeComputationKernels.h"
#include "EcalUncalibRecHitMultiFitAlgoPortable.h"
#include "TimeComputationKernels.h"

//#define DEBUG
//#define ECAL_RECO_ALPAKA_DEBUG

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
        using digis_type = std::vector<uint16_t>;
        using dids_type = std::vector<uint32_t>;
        // according to the cpu setup  //----> hardcoded
        bool constexpr gainSwitchUseMaxSampleEB = true;
        // according to the cpu setup  //----> hardcoded
        bool constexpr gainSwitchUseMaxSampleEE = false;

        auto const totalChannels = static_cast<uint32_t>(digisDevEB->metadata().size() + digisDevEE->metadata().size()); // FIXME: actual size

        EventDataForScratchGPU scratch;
        scratch.allocate(configParams, totalChannels, queue);

        //
        // 1d preparation kernel
        //
        uint32_t constexpr nchannels_per_block = 32;
        auto constexpr threads_1d = EcalDataFrame::MAXSAMPLES * nchannels_per_block;
        auto const blocks_1d = threads_1d > EcalDataFrame::MAXSAMPLES * totalChannels ? 1u : (totalChannels * EcalDataFrame::MAXSAMPLES + threads_1d - 1) / threads_1d;
        auto workDivPrep1D = cms::alpakatools::make_workdiv<Acc1D>(blocks_1d, threads_1d);
        alpaka::exec<Acc1D>(queue,
                            workDivPrep1D,
                            kernel_prep_1d_and_initialize{},
                            digisDevEB.const_view(),
                            digisDevEE.const_view(),
                            uncalibRecHitsDevEB.view(),
                            uncalibRecHitsDevEE.view(),
                            conditionsDev.const_view(),
                            reinterpret_cast<::ecal::multifit::SampleVector*>(scratch.samplesDevBuf.value().data()),
                            reinterpret_cast<::ecal::multifit::SampleGainVector*>(scratch.gainsNoiseDevBuf.value().data()),
                            scratch.hasSwitchToGain6DevBuf.value().data(),
                            scratch.hasSwitchToGain1DevBuf.value().data(),
                            scratch.isSaturatedDevBuf.value().data(),
                            scratch.acStateDevBuf.value().data(),
                            reinterpret_cast<::ecal::multifit::BXVectorType*>(scratch.activeBXsDevBuf.value().data()),
                            gainSwitchUseMaxSampleEB,
                            gainSwitchUseMaxSampleEE
                           );

        //
        // 2d preparation kernel
        //
        Vec2D const blocks_2d{1u, totalChannels};
        Vec2D const threads_2d{EcalDataFrame::MAXSAMPLES, EcalDataFrame::MAXSAMPLES};
        auto workDivPrep2D = cms::alpakatools::make_workdiv<Acc2D>(blocks_2d, threads_2d);
        alpaka::exec<Acc2D>(queue,
                            workDivPrep2D,
                            kernel_prep_2d{},
                            digisDevEB.const_view(),
                            digisDevEE.const_view(),
                            conditionsDev.const_view(),
                            reinterpret_cast<::ecal::multifit::SampleGainVector*>(scratch.gainsNoiseDevBuf.value().data()),
                            reinterpret_cast<::ecal::multifit::SampleMatrix*>(scratch.noisecovDevBuf.value().data()),
                            reinterpret_cast<::ecal::multifit::PulseMatrixType*>(scratch.pulse_matrixDevBuf.value().data()),
                            scratch.hasSwitchToGain6DevBuf.value().data(),
                            scratch.hasSwitchToGain1DevBuf.value().data(),
                            scratch.isSaturatedDevBuf.value().data());

        // run minimization kernels
        v1::minimization_procedure(digisDevEB, digisDevEE, uncalibRecHitsDevEB, uncalibRecHitsDevEE, scratch, conditionsDev, configParams, queue);

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
                              digisDevEE.const_view(),
                              conditionsDev.const_view(),
                              scratch.sample_valuesDevBuf.value().data(),
                              scratch.sample_value_errorsDevBuf.value().data(),
                              scratch.ampMaxErrorDevBuf.value().data(),
                              scratch.useless_sample_valuesDevBuf.value().data(),
                              scratch.pedestal_numsDevBuf.value().data());

          //
          // TODO: small kernel only for EB. It needs to be checked if
          /// fusing such small kernels is beneficial in here
          //
          // we are running only over EB digis
          // therefore we need to create threads/blocks only for that
          auto const threadsFixMGPA = threads_1d;
          auto const blocksFixMGPA =
              threadsFixMGPA > EcalDataFrame::MAXSAMPLES * static_cast<unsigned int>(digisDevEB->metadata().size()) // FIXME check if metadata.size is OK here
                  ? 1
                  : (EcalDataFrame::MAXSAMPLES * static_cast<unsigned int>(digisDevEB->metadata().size()) + threadsFixMGPA - 1) / threadsFixMGPA;
          auto workDivTimeFixMGPAslew1D = cms::alpakatools::make_workdiv<Acc1D>(blocksFixMGPA, threadsFixMGPA);
          alpaka::exec<Acc1D>(queue,
                              workDivTimeFixMGPAslew1D,
                              kernel_time_compute_fixMGPAslew{},
                              digisDevEB.const_view(),
                              digisDevEE.const_view(),
                              conditionsDev.const_view(),
                              scratch.sample_valuesDevBuf.value().data(),
                              scratch.sample_value_errorsDevBuf.value().data(),
                              scratch.useless_sample_valuesDevBuf.value().data());
                  
          auto const threads_nullhypot = threads_1d;
          auto const blocks_nullhypot = blocks_1d;
          auto workDivTimeNullhypot1D = cms::alpakatools::make_workdiv<Acc1D>(blocks_nullhypot, threads_nullhypot);
          alpaka::exec<Acc1D>(queue,
                              workDivTimeNullhypot1D,
                              kernel_time_compute_nullhypot{},
                              digisDevEB.const_view(),
                              digisDevEE.const_view(),
                              scratch.sample_valuesDevBuf.value().data(),
                              scratch.sample_value_errorsDevBuf.value().data(),
                              scratch.useless_sample_valuesDevBuf.value().data(),
                              scratch.chi2sNullHypotDevBuf.value().data(),
                              scratch.sum0sNullHypotDevBuf.value().data(),
                              scratch.sumAAsNullHypotDevBuf.value().data());

          unsigned int const nchannels_per_block_makeratio = 10;
          auto const threads_makeratio = 45 * nchannels_per_block_makeratio;
          unsigned int const blocks_makeratio = threads_makeratio > 45 * totalChannels
                                              ? 1
                                              : (totalChannels * 45 + threads_makeratio - 1) / threads_makeratio;
          auto workDivTimeMakeRatio1D = cms::alpakatools::make_workdiv<Acc1D>(blocks_makeratio, threads_makeratio);
          alpaka::exec<Acc1D>(queue,
                              workDivTimeMakeRatio1D,
                              kernel_time_compute_makeratio{},
                              digisDevEB.const_view(),
                              digisDevEE.const_view(),
                              scratch.sample_valuesDevBuf.value().data(),
                              scratch.sample_value_errorsDevBuf.value().data(),
                              scratch.useless_sample_valuesDevBuf.value().data(),
                              scratch.pedestal_numsDevBuf.value().data(),
                              scratch.sumAAsNullHypotDevBuf.value().data(),
                              scratch.sum0sNullHypotDevBuf.value().data(),
                              scratch.tMaxAlphaBetasDevBuf.value().data(),
                              scratch.tMaxErrorAlphaBetasDevBuf.value().data(),
                              scratch.accTimeMaxDevBuf.value().data(),
                              scratch.accTimeWgtDevBuf.value().data(),
                              scratch.tcStateDevBuf.value().data(),
                              configParams.amplitudeFitParametersEB,
                              configParams.amplitudeFitParametersEE,
                              configParams.timeFitParametersEB,
                              configParams.timeFitParametersEE,
                              configParams.timeFitParametersSizeEB,
                              configParams.timeFitParametersSizeEE,
                              configParams.timeFitLimitsFirstEB,
                              configParams.timeFitLimitsFirstEE,
                              configParams.timeFitLimitsSecondEB,
                              configParams.timeFitLimitsSecondEE);

          auto const threads_findamplchi2 = threads_1d;
          auto const blocks_findamplchi2 = blocks_1d;
          auto workDivTimeFindAmplChi21D = cms::alpakatools::make_workdiv<Acc1D>(blocks_findamplchi2, threads_findamplchi2);
          alpaka::exec<Acc1D>(queue,
                              workDivTimeFindAmplChi21D,
                              kernel_time_compute_findamplchi2_and_finish{},
                              digisDevEB.const_view(),
                              digisDevEE.const_view(),
                              scratch.sample_valuesDevBuf.value().data(),
                              scratch.sample_value_errorsDevBuf.value().data(),
                              scratch.useless_sample_valuesDevBuf.value().data(),
                              scratch.tMaxAlphaBetasDevBuf.value().data(),
                              scratch.tMaxErrorAlphaBetasDevBuf.value().data(),
                              scratch.accTimeMaxDevBuf.value().data(),
                              scratch.accTimeWgtDevBuf.value().data(),
                              scratch.sumAAsNullHypotDevBuf.value().data(),
                              scratch.sum0sNullHypotDevBuf.value().data(),
                              scratch.chi2sNullHypotDevBuf.value().data(),
                              scratch.tcStateDevBuf.value().data(),
                              scratch.ampMaxAlphaBetaDevBuf.value().data(),
                              scratch.ampMaxErrorDevBuf.value().data(),
                              scratch.timeMaxDevBuf.value().data(),
                              scratch.timeErrorDevBuf.value().data(),
                              configParams.amplitudeFitParametersEB,
                              configParams.amplitudeFitParametersEE);

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
                              uncalibRecHitsDevEE.view(),
                              conditionsDev.const_view(),
                              scratch.timeMaxDevBuf.value().data(),
                              scratch.timeErrorDevBuf.value().data(),
                              configParams.timeConstantTermEB,
                              configParams.timeConstantTermEE,
                              configParams.timeNconstEB,
                              configParams.timeNconstEE,
                              configParams.amplitudeThreshEB,
                              configParams.amplitudeThreshEE,
                              configParams.outOfTimeThreshG12pEB,
                              configParams.outOfTimeThreshG12pEE,
                              configParams.outOfTimeThreshG12mEB,
                              configParams.outOfTimeThreshG12mEE,
                              configParams.outOfTimeThreshG61pEB,
                              configParams.outOfTimeThreshG61pEE,
                              configParams.outOfTimeThreshG61mEB,
                              configParams.outOfTimeThreshG61mEE);
        }
      }

    }  // namespace multifit
  }  // namespace ecal
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
