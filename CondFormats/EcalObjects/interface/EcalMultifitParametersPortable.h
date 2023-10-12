#ifndef CondFormats_EcalObjects_interface_EcalMultifitParametersPortable_h
#define CondFormats_EcalObjects_interface_EcalMultifitParametersPortable_h

#include "CondFormats/EcalObjects/interface/EcalMultifitParametersSoA.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

using EcalMultifitParametersPortableHost = PortableHostCollection<EcalMultifitParametersSoA>;

#endif
