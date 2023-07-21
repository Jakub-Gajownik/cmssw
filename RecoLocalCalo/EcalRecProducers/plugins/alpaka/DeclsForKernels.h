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
 
      //  cms::cuda::device::unique_ptr<SMT[]> noisecov;
      //  cms::cuda::device::unique_ptr<PMT[]> pulse_matrix;
        std::optional<cms::alpakatools::device_buffer<Device, BXVT[]>> activeBXsDevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, char[]>> acStateDevBuf;
 
        std::optional<cms::alpakatools::device_buffer<Device, bool[]>> hasSwitchToGain6DevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, bool[]>> hasSwitchToGain1DevBuf;
        std::optional<cms::alpakatools::device_buffer<Device, bool[]>> isSaturatedDevBuf;
 
      //  cms::cuda::device::unique_ptr<SVT[]> sample_values, sample_value_errors;
      //  cms::cuda::device::unique_ptr<bool[]> useless_sample_values;
      //  cms::cuda::device::unique_ptr<SVT[]> chi2sNullHypot;
      //  cms::cuda::device::unique_ptr<SVT[]> sum0sNullHypot;
      //  cms::cuda::device::unique_ptr<SVT[]> sumAAsNullHypot;
      //  cms::cuda::device::unique_ptr<char[]> pedestal_nums;
      //  cms::cuda::device::unique_ptr<SVT[]> tMaxAlphaBetas, tMaxErrorAlphaBetas;
      //  cms::cuda::device::unique_ptr<SVT[]> accTimeMax, accTimeWgt;
      //  cms::cuda::device::unique_ptr<SVT[]> ampMaxAlphaBeta, ampMaxError;
      //  cms::cuda::device::unique_ptr<SVT[]> timeMax, timeError;
      //  cms::cuda::device::unique_ptr<TimeComputationState[]> tcState;
 
        void allocate(ConfigurationParameters const& configParameters,
                      uint32_t size,
                      Queue& queue) {
          constexpr auto svlength = getLength<::ecal::multifit::SampleVector>();
          constexpr auto sgvlength = getLength<::ecal::multifit::SampleGainVector>();
          //constexpr auto smlength = getLength<::ecal::multifit::SampleMatrix>();
          //constexpr auto pmlength = getLength<::ecal::multifit::PulseMatrixType>();
          constexpr auto bxvlength = getLength<::ecal::multifit::BXVectorType>();
 
          //auto alloc = [queue](auto& var, uint32_t size) {
          //  using element_type = typename std::remove_reference<decltype(var)>::type;
          //  var = cms::alpakatools::make_device_buffer<element_type[]>(queue, size);
          //};
 
          samplesDevBuf = cms::alpakatools::make_device_buffer<SVT[]>(queue, size * svlength);
          gainsNoiseDevBuf = cms::alpakatools::make_device_buffer<SGVT[]>(queue, size * sgvlength);
 
      //    alloc(noisecov, size * smlength);
      //    alloc(pulse_matrix, size * pmlength);
          activeBXsDevBuf = cms::alpakatools::make_device_buffer<BXVT[]>(queue, size * bxvlength);
          acStateDevBuf = cms::alpakatools::make_device_buffer<char[]>(queue, size);
 
          hasSwitchToGain6DevBuf = cms::alpakatools::make_device_buffer<bool[]>(queue, size);
          hasSwitchToGain1DevBuf = cms::alpakatools::make_device_buffer<bool[]>(queue, size);
          isSaturatedDevBuf = cms::alpakatools::make_device_buffer<bool[]>(queue, size);
 
          if (configParameters.shouldRunTimingComputation) {
      //      alloc(sample_values, size * svlength);
      //      alloc(sample_value_errors, size * svlength);
      //      alloc(useless_sample_values, size * EcalDataFrame::MAXSAMPLES);
      //      alloc(chi2sNullHypot, size);
      //      alloc(sum0sNullHypot, size);
      //      alloc(sumAAsNullHypot, size);
      //      alloc(pedestal_nums, size);
 
      //      alloc(tMaxAlphaBetas, size);
      //      alloc(tMaxErrorAlphaBetas, size);
      //      alloc(accTimeMax, size);
      //      alloc(accTimeWgt, size);
      //      alloc(ampMaxAlphaBeta, size);
      //      alloc(ampMaxError, size);
      //      alloc(timeMax, size);
      //      alloc(timeError, size);
      //      alloc(tcState, size);
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
