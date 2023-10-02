#ifndef DataFormats_EcalDigi_EcalDigiSoA_h
#define DataFormats_EcalDigi_EcalDigiSoA_h

#include "DataFormats/EcalDigi/interface/EcalConstants.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace ecal {

  // due to a ROOT limitation the std::array needs to be wrapped in a struct
  // https://github.com/root-project/root/issues/12007
  using DataArrayStruct = StdArrayStruct<uint16_t, ecalPh1::sampleSize>;

  GENERATE_SOA_LAYOUT(EcalDigiSoALayout,
                      SOA_COLUMN(uint32_t, id),
                      SOA_COLUMN(DataArrayStruct, data),
                      SOA_SCALAR(uint32_t, size))

  using EcalDigiSoA = EcalDigiSoALayout<>;
}  // namespace ecal

#endif
