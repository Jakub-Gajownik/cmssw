#ifndef DataFormats_EcalDigi_EcalDigiPhase2SoA_h
#define DataFormats_EcalDigi_EcalDigiPhase2SoA_h
 
#include "DataFormats/EcalDigi/interface/EcalConstants.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
 
#include <array> 
 
namespace ecal {
 
  using DataArrayPhase2 = std::array<uint16_t, ecalPh2::sampleSize>;
 
  GENERATE_SOA_LAYOUT(EcalDigiPhase2SoALayout,
    SOA_COLUMN(uint32_t, id),
    SOA_COLUMN(DataArrayPhase2, data),
    SOA_SCALAR(uint32_t, size)
  )

  using EcalDigiPhase2SoA = EcalDigiPhase2SoALayout<>;
}

#endif
