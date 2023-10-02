#ifndef DataFormats_EcalRecHit_EcalUncalibratedRecHitHostCollection_h
#define DataFormats_EcalRecHit_EcalUncalibratedRecHitHostCollection_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHitSoA.h"

namespace ecal {

  // EcalUncalibratedRecHitSoA in host memory
  using UncalibratedRecHitHostCollection = PortableHostCollection<EcalUncalibratedRecHitSoA>;
}  // namespace ecal

#endif