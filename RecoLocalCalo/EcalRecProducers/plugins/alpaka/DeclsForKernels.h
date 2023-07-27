#ifndef RecoLocalCalo_EcalRecProducers_plugins_alpaka_DeclsForKernels_h
#define RecoLocalCalo_EcalRecProducers_plugins_alpaka_DeclsForKernels_h

#include <vector>

#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatusCode.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
//#include "CondFormats/EcalObjects/interface/EcalGainRatiosGPU.h"
//#include "CondFormats/EcalObjects/interface/EcalIntercalibConstantsGPU.h"
//#include "CondFormats/EcalObjects/interface/EcalLaserAPDPNRatiosGPU.h"
//#include "CondFormats/EcalObjects/interface/EcalLaserAPDPNRatiosRefGPU.h"
//#include "CondFormats/EcalObjects/interface/EcalLaserAlphasGPU.h"
//#include "CondFormats/EcalObjects/interface/EcalLinearCorrectionsGPU.h"
//#include "CondFormats/EcalObjects/interface/EcalMultifitParametersGPU.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
//#include "CondFormats/EcalObjects/interface/EcalPedestalsGPU.h"
//#include "CondFormats/EcalObjects/interface/EcalPulseCovariancesGPU.h"
//#include "CondFormats/EcalObjects/interface/EcalPulseShapesGPU.h"
//#include "CondFormats/EcalObjects/interface/EcalRechitADCToGeVConstantGPU.h"
//#include "CondFormats/EcalObjects/interface/EcalRechitChannelStatusGPU.h"
//#include "CondFormats/EcalObjects/interface/EcalSamplesCorrelationGPU.h"
#include "CondFormats/EcalObjects/interface/EcalTimeBiasCorrections.h"
//#include "CondFormats/EcalObjects/interface/EcalTimeBiasCorrectionsGPU.h"
//#include "CondFormats/EcalObjects/interface/EcalTimeCalibConstantsGPU.h"
#include "CondFormats/EcalObjects/interface/EcalTimeOffsetConstant.h"
#include "CondFormats/EcalObjects/interface/EcalWeightSet.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"

#include "../EigenMatrixTypes_gpu.h"

struct EcalPulseShape;
class EcalSampleMask;
class EcalTimeBiasCorrections;
struct EcalPulseCovariance;
class EcalDigiCollection;
class EcalXtalGroupId;
class EcalSamplesCorrelation;

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace ecal {
    namespace multifit {

      enum class TimeComputationState : char { NotFinished = 0, Finished = 1 };
      enum class MinimizationState : char {
        NotFinished = 0,
        Finished = 1,
        Precomputed = 2,
      };

      // parameters have a fixed type
      // Can we go by with single precision
      struct ConfigurationParameters {
        using type = double;
        // device ptrs
        const type *amplitudeFitParametersEB = nullptr, *amplitudeFitParametersEE = nullptr;

        uint32_t timeFitParametersSizeEB, timeFitParametersSizeEE;
        // device ptrs
        const type *timeFitParametersEB = nullptr, *timeFitParametersEE = nullptr;

        type timeFitLimitsFirstEB, timeFitLimitsFirstEE;
        type timeFitLimitsSecondEB, timeFitLimitsSecondEE;

        type timeConstantTermEB, timeConstantTermEE;

        type timeNconstEB, timeNconstEE;

        type amplitudeThreshEE, amplitudeThreshEB;

        type outOfTimeThreshG12pEB, outOfTimeThreshG12mEB;
        type outOfTimeThreshG12pEE, outOfTimeThreshG12mEE;
        type outOfTimeThreshG61pEE, outOfTimeThreshG61mEE;
        type outOfTimeThreshG61pEB, outOfTimeThreshG61mEB;

        std::array<uint32_t, 3> kernelMinimizeThreads;

        bool shouldRunTimingComputation;
      };

      template <typename EigenM>
      constexpr auto getLength() -> uint32_t {
        return EigenM::RowsAtCompileTime * EigenM::ColsAtCompileTime;
      }

      struct EventDataForScratchGPU {
        using SVT = ::ecal::multifit::SampleVector::Scalar;
        using SGVT = ::ecal::multifit::SampleGainVector::Scalar;
        using SMT = ::ecal::multifit::SampleMatrix::Scalar;
        using PMT = ::ecal::multifit::PulseMatrixType::Scalar;
        using BXVT = ::ecal::multifit::BXVectorType::Scalar;

        EventDataForScratchGPU() = default;

        std::optional<cms::alpakatools::device_buffer<Device, SVT[]>> samplesDevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, SGVT[]>> gainsNoiseDevBuf;

        std::optional<cms::alpakatools::device_buffer<Device, SMT[]>> noisecovDevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, PMT[]>> pulse_matrixDevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, BXVT[]>> activeBXsDevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, char[]>> acStateDevBuf;

        std::optional<cms::alpakatools::device_buffer<Device, bool[]>> hasSwitchToGain6DevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, bool[]>> hasSwitchToGain1DevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, bool[]>> isSaturatedDevBuf;

        std::optional<cms::alpakatools::device_buffer<Device, SVT[]>> sample_valuesDevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, SVT[]>> sample_value_errorsDevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, bool[]>> useless_sample_valuesDevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, SVT[]>> chi2sNullHypotDevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, SVT[]>> sum0sNullHypotDevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, SVT[]>> sumAAsNullHypotDevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, char[]>> pedestal_numsDevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, SVT[]>> tMaxAlphaBetasDevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, SVT[]>> tMaxErrorAlphaBetasDevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, SVT[]>> accTimeMaxDevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, SVT[]>> accTimeWgtDevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, SVT[]>> ampMaxAlphaBetaDevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, SVT[]>> ampMaxErrorDevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, SVT[]>> timeMaxDevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, SVT[]>> timeErrorDevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, TimeComputationState[]>> tcStateDevBuf;

        void allocate(ConfigurationParameters const& configParameters,
                      uint32_t size,
                      Queue& queue) {
          constexpr auto svlength = getLength<::ecal::multifit::SampleVector>();
          constexpr auto sgvlength = getLength<::ecal::multifit::SampleGainVector>();
          constexpr auto smlength = getLength<::ecal::multifit::SampleMatrix>();
          constexpr auto pmlength = getLength<::ecal::multifit::PulseMatrixType>();
          constexpr auto bxvlength = getLength<::ecal::multifit::BXVectorType>();

          samplesDevBuf = cms::alpakatools::make_device_buffer<SVT[]>(queue, size * svlength);
          gainsNoiseDevBuf = cms::alpakatools::make_device_buffer<SGVT[]>(queue, size * sgvlength);

          noisecovDevBuf = cms::alpakatools::make_device_buffer<SMT[]>(queue, size * smlength);
          pulse_matrixDevBuf = cms::alpakatools::make_device_buffer<PMT[]>(queue, size * pmlength);
          activeBXsDevBuf = cms::alpakatools::make_device_buffer<BXVT[]>(queue, size * bxvlength);
          acStateDevBuf = cms::alpakatools::make_device_buffer<char[]>(queue, size);

          hasSwitchToGain6DevBuf = cms::alpakatools::make_device_buffer<bool[]>(queue, size);
          hasSwitchToGain1DevBuf = cms::alpakatools::make_device_buffer<bool[]>(queue, size);
          isSaturatedDevBuf = cms::alpakatools::make_device_buffer<bool[]>(queue, size);

          if (configParameters.shouldRunTimingComputation) {
            sample_valuesDevBuf = cms::alpakatools::make_device_buffer<SVT[]>(queue, size * svlength);
            sample_value_errorsDevBuf = cms::alpakatools::make_device_buffer<SVT[]>(queue, size * svlength);
            useless_sample_valuesDevBuf = cms::alpakatools::make_device_buffer<bool[]>(queue, size * EcalDataFrame::MAXSAMPLES);
            chi2sNullHypotDevBuf = cms::alpakatools::make_device_buffer<SVT[]>(queue, size);
            sum0sNullHypotDevBuf = cms::alpakatools::make_device_buffer<SVT[]>(queue, size);
            sumAAsNullHypotDevBuf = cms::alpakatools::make_device_buffer<SVT[]>(queue, size);
            pedestal_numsDevBuf = cms::alpakatools::make_device_buffer<char[]>(queue, size);

            tMaxAlphaBetasDevBuf = cms::alpakatools::make_device_buffer<SVT[]>(queue, size);
            tMaxErrorAlphaBetasDevBuf = cms::alpakatools::make_device_buffer<SVT[]>(queue, size);
            accTimeMaxDevBuf = cms::alpakatools::make_device_buffer<SVT[]>(queue, size);
            accTimeWgtDevBuf = cms::alpakatools::make_device_buffer<SVT[]>(queue, size);
            ampMaxAlphaBetaDevBuf = cms::alpakatools::make_device_buffer<SVT[]>(queue, size);
            ampMaxErrorDevBuf = cms::alpakatools::make_device_buffer<SVT[]>(queue, size);
            timeMaxDevBuf = cms::alpakatools::make_device_buffer<SVT[]>(queue, size);
            timeErrorDevBuf = cms::alpakatools::make_device_buffer<SVT[]>(queue, size);
            tcStateDevBuf = cms::alpakatools::make_device_buffer<TimeComputationState[]>(queue, size);
          }
        }
      };

      // const refs products to conditions
      struct ConditionsProducts {
        //EcalPedestalsGPU::Product const& pedestals;
        //EcalGainRatiosGPU::Product const& gainRatios;
        //EcalPulseShapesGPU::Product const& pulseShapes;
        //EcalPulseCovariancesGPU::Product const& pulseCovariances;
        //EcalSamplesCorrelationGPU::Product const& samplesCorrelation;
        //EcalTimeBiasCorrectionsGPU::Product const& timeBiasCorrections;
        //EcalTimeCalibConstantsGPU::Product const& timeCalibConstants;
        EcalSampleMask const& sampleMask;
        EcalTimeOffsetConstant const& timeOffsetConstant;
        uint32_t offsetForHashes;
        //EcalMultifitParametersGPU::Product const& multifitParameters;
      };

      struct xyz {
        int x, y, z;
      };

      struct conf_data {
        xyz threads;
        bool runV1;
        //cudaStream_t cuStream;
      };

    }  // namespace multifit
  }  // namespace ecal

  //
  // ECAL Rechit producer
  //

  namespace ecal {
    namespace rechit {

      // parameters that are read in the configuration file for rechit producer
      struct ConfigurationParameters {
        // device ptrs
        const int* ChannelStatusToBeExcluded = nullptr;
        uint32_t ChannelStatusToBeExcludedSize;

        bool killDeadChannels;

        bool recoverEBIsolatedChannels;
        bool recoverEEIsolatedChannels;
        bool recoverEBVFE;
        bool recoverEEVFE;
        bool recoverEBFE;
        bool recoverEEFE;

        float EBLaserMIN;
        float EELaserMIN;
        float EBLaserMAX;
        float EELaserMAX;

        const int* expanded_v_DB_reco_flags;
        const uint32_t* expanded_Sizes_v_DB_reco_flags;
        const uint32_t* expanded_flagbit_v_DB_reco_flags;
        uint32_t expanded_v_DB_reco_flagsSize;

        uint32_t flagmask;
      };

      // const refs products to conditions
      struct ConditionsProducts {
        //EcalRechitADCToGeVConstantGPU::Product const& ADCToGeV;
        //EcalIntercalibConstantsGPU::Product const& Intercalib;
        //EcalRechitChannelStatusGPU::Product const& ChannelStatus;

        //EcalLaserAPDPNRatiosGPU::Product const& LaserAPDPNRatios;
        //EcalLaserAPDPNRatiosRefGPU::Product const& LaserAPDPNRatiosRef;
        //EcalLaserAlphasGPU::Product const& LaserAlphas;
        //EcalLinearCorrectionsGPU::Product const& LinearCorrections;

        uint32_t offsetForHashes;
      };

    }  // namespace rechit
  }  // namespace ecal
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoLocalCalo_EcalRecProducers_plugins_DeclsForKernels_h
