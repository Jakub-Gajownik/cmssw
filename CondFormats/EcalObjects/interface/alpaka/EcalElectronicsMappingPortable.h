#ifndef CondFormats_EcalObjects_interface_alpaka_EcalElectronicsMappingPortable_h
#define CondFormats_EcalObjects_interface_alpaka_EcalElectronicsMappingPortable_h

#include "CondFormats/EcalObjects/interface/EcalElectronicsMappingPortable.h"
#include "CondFormats/EcalObjects/interface/EcalElectronicsMappingSoA.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using EcalElectronicsMappingPortableHost = EcalElectronicsMappingPortableHost;
  using EcalElectronicsMappingPortableDevice = PortableCollection<EcalElectronicsMappingSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif
