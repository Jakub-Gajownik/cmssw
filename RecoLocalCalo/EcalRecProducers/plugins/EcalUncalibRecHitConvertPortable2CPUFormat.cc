#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/EmptyGroupDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHitHostCollection.h"


class EcalUncalibRecHitConvertPortable2CPUFormat : public edm::stream::EDProducer<> {
public:
  explicit EcalUncalibRecHitConvertPortable2CPUFormat(edm::ParameterSet const &ps);
  ~EcalUncalibRecHitConvertPortable2CPUFormat() override;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  using InputProduct = EcalUncalibratedRecHitHostCollection;
  void produce(edm::Event&, edm::EventSetup const&) override;

private:
  const bool isPhase2_;
  const edm::EDGetTokenT<InputProduct> recHitsGPUEB_;
  const edm::EDGetTokenT<InputProduct> recHitsGPUEE_;
  const edm::EDPutTokenT<EBUncalibratedRecHitCollection> uncalibRecHitsCPUEBToken_;
  const edm::EDPutTokenT<EEUncalibratedRecHitCollection> uncalibRecHitsCPUEEToken_;
};

void EcalUncalibRecHitConvertPortable2CPUFormat::fillDescriptions(edm::ConfigurationDescriptions &confDesc) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("recHitsLabelGPUEB", edm::InputTag("ecalUncalibRecHitProducerGPU", "EcalUncalibRecHitsEB"));
  desc.add<std::string>("recHitsLabelCPUEB", "EcalUncalibRecHitsEB");
  desc.ifValue(
      edm::ParameterDescription<bool>("isPhase2", false, true),
      false >>
              (edm::ParameterDescription<edm::InputTag>(
                   "recHitsLabelGPUEE", edm::InputTag("ecalUncalibRecHitProducerGPU", "EcalUncalibRecHitsEE"), true) and
               edm::ParameterDescription<std::string>("recHitsLabelCPUEE", "EcalUncalibRecHitsEE", true)) or
          true >> edm::EmptyGroupDescription());
  confDesc.add("ecalUncalibRecHitConvertPortable2CPUFormat", desc);
}

EcalUncalibRecHitConvertPortable2CPUFormat::EcalUncalibRecHitConvertPortable2CPUFormat(edm::ParameterSet const &ps)
    : isPhase2_{ps.getParameter<bool>("isPhase2")},
      recHitsGPUEB_{consumes<InputProduct>(ps.getParameter<edm::InputTag>("recHitsLabelGPUEB"))},
      recHitsGPUEE_{isPhase2_ ? edm::EDGetTokenT<InputProduct>{}
                              : consumes<InputProduct>(ps.getParameter<edm::InputTag>("recHitsLabelGPUEE"))},
      uncalibRecHitsCPUEBToken_{produces<EBUncalibratedRecHitCollection>(ps.getParameter<std::string>("recHitsLabelCPUEB"))},
      uncalibRecHitsCPUEEToken_{isPhase2_ ? edm::EDPutTokenT<EEUncalibratedRecHitCollection>{}
                              : produces<EEUncalibratedRecHitCollection>(ps.getParameter<std::string>("recHitsLabelCPUEE"))} {}

EcalUncalibRecHitConvertPortable2CPUFormat::~EcalUncalibRecHitConvertPortable2CPUFormat() {}

void EcalUncalibRecHitConvertPortable2CPUFormat::produce(edm::Event &event, edm::EventSetup const &setup) {
  auto const& uncalRecHitsEBColl = event.get(recHitsGPUEB_);
  auto const& uncalRecHitsEBCollView = uncalRecHitsEBColl.const_view();
  auto recHitsCPUEB = std::make_unique<EBUncalibratedRecHitCollection>();
  recHitsCPUEB->reserve(uncalRecHitsEBCollView.size());

  for (uint32_t i = 0; i < uncalRecHitsEBCollView.size(); ++i) {
    recHitsCPUEB->emplace_back(DetId{uncalRecHitsEBCollView.id()[i]},
                               uncalRecHitsEBCollView.amplitude()[i],
                               uncalRecHitsEBCollView.pedestal()[i],
                               uncalRecHitsEBCollView.jitter()[i],
                               uncalRecHitsEBCollView.chi2()[i],
                               uncalRecHitsEBCollView.flags()[i]);
    if (isPhase2_)
      (*recHitsCPUEB)[i].setAmplitudeError(uncalRecHitsEBCollView.amplitudeError()[i]);
    (*recHitsCPUEB)[i].setJitterError(uncalRecHitsEBCollView.jitterError()[i]);
    for (uint32_t sample = 0; sample < EcalDataFrame::MAXSAMPLES; ++sample)
      (*recHitsCPUEB)[i].setOutOfTimeAmplitude(sample, uncalRecHitsEBCollView.outOfTimeAmplitudes()[i][sample]);
  }

  if (!isPhase2_) {
    auto const& uncalRecHitsEEColl = event.get(recHitsGPUEE_);
    auto const& uncalRecHitsEECollView = uncalRecHitsEEColl.const_view();
    auto recHitsCPUEE = std::make_unique<EEUncalibratedRecHitCollection>();
    recHitsCPUEE->reserve(uncalRecHitsEECollView.size());

    for (uint32_t i = 0; i < uncalRecHitsEECollView.size(); ++i) {
      recHitsCPUEE->emplace_back(DetId{uncalRecHitsEECollView.id()[i]},
                                 uncalRecHitsEECollView.amplitude()[i],
                                 uncalRecHitsEECollView.pedestal()[i],
                                 uncalRecHitsEECollView.jitter()[i],
                                 uncalRecHitsEECollView.chi2()[i],
                                 uncalRecHitsEECollView.flags()[i]);
      (*recHitsCPUEE)[i].setJitterError(uncalRecHitsEECollView.jitterError()[i]);
      for (uint32_t sample = 0; sample < EcalDataFrame::MAXSAMPLES; ++sample) {
        (*recHitsCPUEE)[i].setOutOfTimeAmplitude(sample, uncalRecHitsEECollView.outOfTimeAmplitudes()[i][sample]);
      }
    }
    event.put(uncalibRecHitsCPUEEToken_, std::move(recHitsCPUEE));
  }
  event.put(uncalibRecHitsCPUEBToken_, std::move(recHitsCPUEB));
}

DEFINE_FWK_MODULE(EcalUncalibRecHitConvertPortable2CPUFormat);
