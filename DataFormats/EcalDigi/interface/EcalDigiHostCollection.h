#ifndef DataFormats_EcalDigi_EcalDigiHostCollection_h
#define DataFormats_EcalDigi_EcalDigiHostCollection_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/EcalDigi/interface/EcalDigiSoA.h"

namespace ecal {

  // EcalDigiSoA in host memory
  using DigiHostCollection = PortableHostCollection<EcalDigiSoA>;
}  // namespace ecal

#endif
