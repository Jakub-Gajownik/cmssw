#ifndef RecoLocalCalo_EcalRecProducers_plugins_AmplitudeComputationKernels_h
#define RecoLocalCalo_EcalRecProducers_plugins_AmplitudeComputationKernels_h

#include "DeclsForKernels.h"
#include "EigenMatrixTypes_gpu.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class EcalPulseShape;
  class EcalPulseCovariance;
  class EcalUncalibratedRecHit;
  
  namespace ecal {
    namespace multifit {
  
      namespace v1 {
  
        void minimization_procedure(EventInputDataGPU const& eventInputGPU,
                                    EventOutputDataGPU& eventOutputGPU,
                                    EventDataForScratchGPU& scratch,
                                    ConditionsProducts const& conditions,
                                    ConfigurationParameters const& configParameters,
                                    cudaStream_t cudaStream);
  
      }
  
    }  // namespace multifit
  }  // namespace ecal
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoLocalCalo_EcalRecProducers_plugins_AmplitudeComputationKernels_h
