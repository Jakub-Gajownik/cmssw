#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"
#include "CondFormats/DataRecord/interface/EcalPulseCovariancesRcd.h"
#include "CondFormats/DataRecord/interface/EcalPulseShapesRcd.h"
#include "CondFormats/DataRecord/interface/EcalSampleMaskRcd.h"
#include "CondFormats/DataRecord/interface/EcalSamplesCorrelationRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeBiasCorrectionsRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeCalibConstantsRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeOffsetConstantRcd.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalPulseCovariances.h"
#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"
#include "CondFormats/EcalObjects/interface/EcalSamplesCorrelation.h"
#include "CondFormats/EcalObjects/interface/EcalSampleMask.h"
#include "CondFormats/EcalObjects/interface/EcalTimeBiasCorrections.h"
#include "CondFormats/EcalObjects/interface/EcalTimeCalibConstants.h"
#include "CondFormats/EcalObjects/interface/EcalTimeOffsetConstant.h"

#include "CondFormats/EcalObjects/interface/alpaka/EcalMultifitConditionsPortable.h"
#include "CondFormats/DataRecord/interface/EcalMultifitConditionsRcd.h"

#include "DataFormats/EcalDetId/interface/EcalElectronicsId.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ModuleFactory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  class EcalMultifitConditionsPortableESProducer : public ESProducer {
  public:
    EcalMultifitConditionsPortableESProducer(edm::ParameterSet const& iConfig) {
      auto cc = setWhatProduced(this);
      pedestalsToken_ = cc.consumes();
      gainRatiosToken_ = cc.consumes();
      pulseShapesToken_ = cc.consumes();
      pulseCovariancesToken_ = cc.consumes();
      samplesCorrelationToken_ = cc.consumes();
      timeBiasCorrectionsToken_ = cc.consumes();
      timeCalibConstantsToken_ = cc.consumes();
      sampleMaskToken_ = cc.consumes();
      timeOffsetConstantToken_ = cc.consumes();
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      descriptions.addWithDefaultLabel(desc);
    }

    std::unique_ptr<EcalMultifitConditionsPortableHost> produce(EcalMultifitConditionsRcd const& iRecord) {
      //auto const& mapping = iRecord.get(token_);
      auto product = std::make_unique<EcalMultifitConditionsPortableHost>(10, cms::alpakatools::host());
      return product;
    }

  private:
    edm::ESGetToken<EcalPedestals, EcalPedestalsRcd> pedestalsToken_;
    edm::ESGetToken<EcalGainRatios, EcalGainRatiosRcd> gainRatiosToken_;
    edm::ESGetToken<EcalPulseShapes, EcalPulseShapesRcd> pulseShapesToken_;
    edm::ESGetToken<EcalPulseCovariances, EcalPulseCovariancesRcd> pulseCovariancesToken_;
    edm::ESGetToken<EcalSamplesCorrelation, EcalSamplesCorrelationRcd> samplesCorrelationToken_;
    edm::ESGetToken<EcalTimeBiasCorrections, EcalTimeBiasCorrectionsRcd> timeBiasCorrectionsToken_;
    edm::ESGetToken<EcalTimeCalibConstants, EcalTimeCalibConstantsRcd> timeCalibConstantsToken_;
    edm::ESGetToken<EcalSampleMask, EcalSampleMaskRcd> sampleMaskToken_;
    edm::ESGetToken<EcalTimeOffsetConstant, EcalTimeOffsetConstantRcd> timeOffsetConstantToken_;

  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_EVENTSETUP_ALPAKA_MODULE(EcalMultifitConditionsPortableESProducer);
