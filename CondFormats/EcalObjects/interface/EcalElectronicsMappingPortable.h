#ifndef CondFormats_EcalObjects_interface_EcalElectronicsMappingPortable_h
#define CondFormats_EcalObjects_interface_EcalElectronicsMappingPortable_h

#include "CondFormats/EcalObjects/interface/EcalElectronicsMappingSoA.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

using EcalElectronicsMappingPortableHost = PortableHostCollection<EcalElectronicsMappingSoA>;

#endif
