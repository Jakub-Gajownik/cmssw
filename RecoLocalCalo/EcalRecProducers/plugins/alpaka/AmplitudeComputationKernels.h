#ifndef RecoLocalCalo_EcalRecProducers_plugins_alpaka_AmplitudeComputationKernels_h
#define RecoLocalCalo_EcalRecProducers_plugins_alpaka_AmplitudeComputationKernels_h

#include "CondFormats/EcalObjects/interface/alpaka/EcalMultifitConditionsPortable.h"
#include "DataFormats/EcalDigi/interface/alpaka/EcalDigiDeviceCollection.h"
#include "DataFormats/EcalRecHit/interface/alpaka/EcalUncalibratedRecHitDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "DeclsForKernels.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace ecal {
    namespace multifit {

      namespace v1 {

        using InputProduct = EcalDigiDeviceCollection;
        using OutputProduct = EcalUncalibratedRecHitDeviceCollection;

        void minimization_procedure(InputProduct const&,
                                    InputProduct const&,
                                    OutputProduct&,
                                    OutputProduct&,
                                    EventDataForScratchGPU&,
                                    EcalMultifitConditionsPortableDevice const&,
                                    ConfigurationParameters const&,
                                    uint32_t const,
                                    Queue&);

      }

    }  // namespace multifit
  }  // namespace ecal
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoLocalCalo_EcalRecProducers_plugins_AmplitudeComputationKernels_h
