#ifndef DataFormats_EcalDigi_interface_alpaka_EcalDigiPhase2DeviceCollection_h
#define DataFormats_EcalDigi_interface_alpaka_EcalDigiPhase2DeviceCollection_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/EcalDigi/interface/EcalDigiPhase2SoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace ecal {

    // make the names from the top-level ecal namespace visible for unqualified lookup
    // inside the ALPAKA_ACCELERATOR_NAMESPACE::portabletest namespace
    using namespace ::ecal;

    // EcalDigiPhase2SoA in device global memory
    using DigiPhase2DeviceCollection = PortableCollection<EcalDigiPhase2SoA>;
  }  // namespace ecal

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif