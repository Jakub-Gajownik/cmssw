#ifndef DataFormats_EcalDigi_interface_alpaka_EcalDigiDeviceCollection_h
#define DataFormats_EcalDigi_interface_alpaka_EcalDigiDeviceCollection_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/EcalDigi/interface/EcalDigiSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace ecal {

    // make the names from the top-level ecal namespace visible for unqualified lookup
    // inside the ALPAKA_ACCELERATOR_NAMESPACE::portabletest namespace
    using namespace ::ecal;

    // EcalDigiSoA in device global memory
    using DigiDeviceCollection = PortableCollection<EcalDigiSoA>;
  }  // namespace ecal

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif
