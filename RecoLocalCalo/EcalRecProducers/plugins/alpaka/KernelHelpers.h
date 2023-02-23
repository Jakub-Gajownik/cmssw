#ifndef RecoLocalCalo_EcalRecProducers_plugins_alpaka_KernelHelpers_h
#define RecoLocalCalo_EcalRecProducers_plugins_alpaka_KernelHelpers_h

#include "DataFormats/CaloRecHit/interface/MultifitComputations.h"

#include <cmath>
#include <limits>
#include <type_traits>

#include <alpaka/alpaka.hpp>
#include <Eigen/Dense>

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace ecal {
    namespace reconstruction {
  
      ALPAKA_FN_ACC uint32_t hashedIndexEB(uint32_t id);
  
      ALPAKA_FN_ACC uint32_t hashedIndexEE(uint32_t id);
  
      ALPAKA_FN_ACC int laser_monitoring_region_EB(uint32_t id);
  
      ALPAKA_FN_ACC int laser_monitoring_region_EE(uint32_t id);
  
    }  // namespace reconstruction
  }  // namespace ecal
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoLocalCalo_EcalRecProducers_plugins_KernelHelpers_h
