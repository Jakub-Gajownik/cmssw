#ifndef EventFilter_EcalRawToDigi_plugins_alpaka_UnpackPortable_h
#define EventFilter_EcalRawToDigi_plugins_alpaka_UnpackPortable_h

#include "CondFormats/EcalObjects/interface/alpaka/EcalElectronicsMappingPortable.h"
#include "DataFormats/EcalDigi/interface/alpaka/EcalDigiDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "DeclsForKernels.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace ecal {
    namespace raw {

      void unpackRaw(InputDataHost const&,
                     InputDataDevice&,
                     EcalDigiDeviceCollection&,
                     EcalDigiDeviceCollection&,
                     EcalElectronicsMappingPortableDevice const&,
                     Queue&,
                     uint32_t const,
                     uint32_t const);

    }  // namespace raw
  }  // namespace ecal
} // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // EventFilter_EcalRawToDigi_plugins_alpaka_UnpackPortable_h
