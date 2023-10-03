#ifndef RecoLocalCalo_EcalRecProducers_plugins_alpaka_EcalUncalibRecHitMultiFitAlgoPortable_h
#define RecoLocalCalo_EcalRecProducers_plugins_alpaka_EcalUncalibRecHitMultiFitAlgoPortable_h

#include <vector>

#include "CondFormats/EcalObjects/interface/alpaka/EcalMultifitConditionsPortable.h"
#include "DataFormats/EcalDigi/interface/alpaka/EcalDigiDeviceCollection.h"
#include "DataFormats/EcalRecHit/interface/alpaka/EcalUncalibratedRecHitDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "DeclsForKernels.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace ecal {
    namespace multifit {

      using InputProduct = EcalDigiDeviceCollection;
      using OutputProduct = EcalUncalibratedRecHitDeviceCollection;
  
      void entryPoint(InputProduct const&,
                      InputProduct const&,
                      OutputProduct&,
                      OutputProduct&,
                      EcalMultifitConditionsPortableDevice const&,
                      ConfigurationParameters const&,
                      Queue&);
  
    }  // namespace multifit
  }  // namespace ecal
} // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoLocalCalo_EcalRecProducers_plugins_alpaka_EcalUncalibRecHitMultiFitAlgoGPU_h
