#ifndef DataFormats_EcalRecHit_alpaka_EcalUncalibratedRecHitDeviceCollection_h
#define DataFormats_EcalRecHit_alpaka_EcalUncalibratedRecHitDeviceCollection_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHitSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace ecal {

    // make the names from the top-level ecal namespace visible for unqualified lookup
    // inside the ALPAKA_ACCELERATOR_NAMESPACE::portabletest namespace
    using namespace ::ecal;

    // EcalUncalibratedRecHitSoA in device global memory
    using UncalibratedRecHitDeviceCollection = PortableCollection<EcalUncalibratedRecHitSoA>;
  }  // namespace ecal

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif