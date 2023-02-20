#ifndef DataFormats_EcalRecHit_EcalUncalibratedRecHitSoA_h
#define DataFormats_EcalRecHit_EcalUncalibratedRecHitSoA_h
                                                                                                                  
#include "DataFormats/SoATemplate/interface/SoALayout.h"                                                          
                                                                                                              
namespace ecal {
                                                                                                                                                                     
  GENERATE_SOA_LAYOUT(EcalUncalibratedRecHitSoALayout,                                                                          
    SOA_COLUMN(uint32_t, id),
    SOA_SCALAR(uint32_t, size),
    SOA_COLUMN(float, amplitude),
    SOA_COLUMN(float, amplitudeError),
    SOA_COLUMN(float, pedestal),
    SOA_COLUMN(float, jitter),
    SOA_COLUMN(float, chi2),
    SOA_COLUMN(uint32_t, flags),
    SOA_COLUMN(uint32_t, aux)
  )

  using EcalUncalibratedRecHitSoA = EcalUncalibratedRecHitSoALayout<>;
}

#endif