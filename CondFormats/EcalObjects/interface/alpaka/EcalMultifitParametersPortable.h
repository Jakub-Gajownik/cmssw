#ifndef CondFormats_EcalObjects_interface_alpaka_EcalMultifitParametersPortable_h
#define CondFormats_EcalObjects_interface_alpaka_EcalMultifitParametersPortable_h

#include "CondFormats/EcalObjects/interface/EcalMultifitParametersPortable.h"
#include "CondFormats/EcalObjects/interface/EcalMultifitParametersSoA.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using ::EcalMultifitParametersPortableHost;
  using EcalMultifitParametersPortableDevice = PortableCollection<EcalMultifitParametersSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif
