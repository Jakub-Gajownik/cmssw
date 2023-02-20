#ifndef DataFormats_EcalDigi_EcalDigiPhase2SoA_h
#define DataFormats_EcalDigi_EcalDigiPhase2SoA_h
                                                                                                                  
#include "DataFormats/EcalDigi/interface/EcalConstants.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"                                                          
                                                                                                                  
#include <array>                                                                                                  
                                                                                                                  
namespace ecal {
  
  using DataArray = std::array<uint16_t, ecalPh2::sampleSize>;                                                    
                                                                                                                  
  GENERATE_SOA_LAYOUT(EcalDigiSoALayout,                                                                          
    SOA_COLUMN(uint32_t, id),
    SOA_COLUMN(DataArray, data),
    SOA_SCALAR(uint32_t, size)
  )

  using EcalDigiPhase2SoA = EcalDigiSoALayout<>;
}

#endif