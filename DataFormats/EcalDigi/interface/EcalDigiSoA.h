#ifndef DataFormats_EcalDigi_EcalDigiSoA_h
#define DataFormats_EcalDigi_EcalDigiSoA_h

#include "DataFormats/EcalDigi/interface/EcalConstants.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

// due to a ROOT limitation the std::array needs to be wrapped in a struct
// https://github.com/root-project/root/issues/12007
using EcalDataArrayStruct = StdArrayStruct<uint16_t, ecalPh1::sampleSize>;

GENERATE_SOA_LAYOUT(EcalDigiSoALayout,
  SOA_COLUMN(uint32_t, id),
  SOA_COLUMN(EcalDataArrayStruct, data),
  SOA_SCALAR(uint32_t, size)
)

using EcalDigiSoA = EcalDigiSoALayout<>;

#endif
