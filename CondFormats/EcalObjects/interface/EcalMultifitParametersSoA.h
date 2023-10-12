#ifndef CondFormats_EcalObjects_EcalMultifitParametersSoA_h
#define CondFormats_EcalObjects_EcalMultifitParametersSoA_h

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/SoATemplate/interface/SoAView.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"

constexpr size_t kNTimeFitParams = 8;
constexpr size_t kNAmplitudeFitParams = 2;
using TimeFitParamsArrayStruct = StdArrayStruct<float, kNTimeFitParams>;
using AmplitudeFitParamsArrayStruct = StdArrayStruct<float, kNAmplitudeFitParams>;

GENERATE_SOA_LAYOUT(EcalMultifitParametersSoALayout,
                    SOA_SCALAR(TimeFitParamsArrayStruct, timeFitParamsEB),
                    SOA_SCALAR(TimeFitParamsArrayStruct, timeFitParamsEE),
                    SOA_SCALAR(AmplitudeFitParamsArrayStruct, amplitudeFitParamsEB),
                    SOA_SCALAR(AmplitudeFitParamsArrayStruct, amplitudeFitParamsEE))

using EcalMultifitParametersSoA = EcalMultifitParametersSoALayout<>;

#endif
