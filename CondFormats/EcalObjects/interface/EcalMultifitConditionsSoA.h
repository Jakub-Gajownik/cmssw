#ifndef CondFormats_EcalObjects_EcalMultifitConditionsSoA_h
#define CondFormats_EcalObjects_EcalMultifitConditionsSoA_h

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/SoATemplate/interface/SoAView.h"


GENERATE_SOA_LAYOUT(EcalMultifitConditionsSoALayout,
                    SOA_COLUMN(uint32_t, rawid),
                    SOA_COLUMN(float, gainRatios_gain12Over6),
                    SOA_COLUMN(float, gainRatios_gain6Over1)
  )

using EcalMultifitConditionsSoA = EcalMultifitConditionsSoALayout<>;


#endif
