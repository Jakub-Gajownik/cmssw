#ifndef EventFilter_EcalRawToDigi_plugins_alpaka_DeclsForKernels_h
#define EventFilter_EcalRawToDigi_plugins_alpaka_DeclsForKernels_h

#include <optional>

#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace ecal {
    namespace raw {
  
      struct InputDataHost {
        InputDataHost() = default;

        void initialize(size_t size, size_t nFeds) {
          data = cms::alpakatools::make_host_buffer<unsigned char[]>(size);
          offsets = cms::alpakatools::make_host_buffer<uint32_t[]>(nFeds);
          feds = cms::alpakatools::make_host_buffer<int[]>(nFeds);
        }

        std::optional<cms::alpakatools::host_buffer<unsigned char[]>> data;
        std::optional<cms::alpakatools::host_buffer<uint32_t[]>> offsets;
        std::optional<cms::alpakatools::host_buffer<int[]>> feds;
      };
  
      struct ConfigurationParameters {
        uint32_t maxChannelsEE;
        uint32_t maxChannelsEB;
      };
  
      struct InputDataDevice {
        InputDataDevice() = default;

        void initialize(const Queue& queue, size_t size, size_t nFeds) {
          data = cms::alpakatools::make_device_buffer<unsigned char[]>(queue, size);
          offsets = cms::alpakatools::make_device_buffer<uint32_t[]>(queue, nFeds);
          feds = cms::alpakatools::make_device_buffer<int[]>(queue, nFeds);
        }

        std::optional<cms::alpakatools::device_buffer<Device, unsigned char[]>> data;
        std::optional<cms::alpakatools::device_buffer<Device, uint32_t[]>> offsets;
        std::optional<cms::alpakatools::device_buffer<Device, int[]>> feds;
      };
  
    }  // namespace raw
  }  // namespace ecal
}  // ALPAKA_ACCELERATOR_NAMESPACE

#endif  // EventFilter_EcalRawToDigi_plugins_alpaka_DeclsForKernels_h
