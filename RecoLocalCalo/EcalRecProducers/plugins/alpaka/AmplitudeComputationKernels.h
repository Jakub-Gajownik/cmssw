#ifndef RecoLocalCalo_EcalRecProducers_plugins_alpaka_AmplitudeComputationKernels_h
#define RecoLocalCalo_EcalRecProducers_plugins_alpaka_AmplitudeComputationKernels_h

#include "CondFormats/EcalObjects/interface/alpaka/EcalMultifitConditionsPortable.h"
#include "DataFormats/EcalDigi/interface/alpaka/EcalDigiDeviceCollection.h"
#include "DataFormats/EcalRecHit/interface/alpaka/EcalUncalibratedRecHitDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "DeclsForKernels.h"
#include "../EigenMatrixTypes_gpu.h"

class EcalPulseShape;
class EcalPulseCovariance;
class EcalUncalibratedRecHit;
  
namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace ecal {
    namespace multifit {
  
      namespace v1 {
  
        void minimization_procedure(DigiDeviceCollection const&,
                                    DigiDeviceCollection const&,
                                    UncalibratedRecHitDeviceCollection&,
                                    UncalibratedRecHitDeviceCollection&,
                                    //EventDataForScratchGPU& scratch,
                                    //EcalMultifitConditionsPortableDevice const&,
                                    //ConditionsProducts const& conditions,
                                    ConfigurationParameters const&,
                                    Queue&);
  
      }
  
    }  // namespace multifit
  }  // namespace ecal
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoLocalCalo_EcalRecProducers_plugins_AmplitudeComputationKernels_h
