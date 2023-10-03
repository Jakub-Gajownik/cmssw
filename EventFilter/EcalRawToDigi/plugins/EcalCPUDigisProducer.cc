#include <utility>

#include "DataFormats/EcalDetId/interface/EcalDetIdCollections.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiHostCollection.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

class EcalCPUDigisProducer : public edm::stream::EDProducer<> {
public:
  explicit EcalCPUDigisProducer(edm::ParameterSet const& ps);
  ~EcalCPUDigisProducer() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, edm::EventSetup const&) override;

  template <typename ProductType, typename... ARGS>
  edm::EDPutTokenT<ProductType> dummyProduces(ARGS&&... args) {
    return (produceDummyIntegrityCollections_) ? produces<ProductType>(std::forward<ARGS>(args)...)
                                               : edm::EDPutTokenT<ProductType>{};
  }

private:
  // input digi collections on host in SoA format
  using InputProduct = EcalDigiHostCollection;
  edm::EDGetTokenT<InputProduct> digisInEBToken_;
  edm::EDGetTokenT<InputProduct> digisInEEToken_;

  // output digi collections in legacy format
  edm::EDPutTokenT<EBDigiCollection> digisOutEBToken_;
  edm::EDPutTokenT<EEDigiCollection> digisOutEEToken_;

  // whether to produce dummy integrity collections
  bool produceDummyIntegrityCollections_;

  // dummy producer collections
  edm::EDPutTokenT<EBSrFlagCollection> ebSrFlagToken_;
  edm::EDPutTokenT<EESrFlagCollection> eeSrFlagToken_;

  // dummy integrity for xtal data
  edm::EDPutTokenT<EBDetIdCollection> ebIntegrityGainErrorsToken_;
  edm::EDPutTokenT<EBDetIdCollection> ebIntegrityGainSwitchErrorsToken_;
  edm::EDPutTokenT<EBDetIdCollection> ebIntegrityChIdErrorsToken_;

  // dummy integrity for xtal data - EE specific (to be rivisited towards EB+EE common collection)
  edm::EDPutTokenT<EEDetIdCollection> eeIntegrityGainErrorsToken_;
  edm::EDPutTokenT<EEDetIdCollection> eeIntegrityGainSwitchErrorsToken_;
  edm::EDPutTokenT<EEDetIdCollection> eeIntegrityChIdErrorsToken_;

  // dummy integrity errors
  edm::EDPutTokenT<EcalElectronicsIdCollection> integrityTTIdErrorsToken_;
  edm::EDPutTokenT<EcalElectronicsIdCollection> integrityZSXtalIdErrorsToken_;
  edm::EDPutTokenT<EcalElectronicsIdCollection> integrityBlockSizeErrorsToken_;

  edm::EDPutTokenT<EcalPnDiodeDigiCollection> pnDiodeDigisToken_;

  // dummy TCC collections
  edm::EDPutTokenT<EcalTrigPrimDigiCollection> ecalTriggerPrimitivesToken_;
  edm::EDPutTokenT<EcalPSInputDigiCollection> ecalPseudoStripInputsToken_;
};

void EcalCPUDigisProducer::fillDescriptions(edm::ConfigurationDescriptions& confDesc) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("digisInLabelEB", edm::InputTag{"ecalRawToDigiPortable", "ebDigis"});
  desc.add<edm::InputTag>("digisInLabelEE", edm::InputTag{"ecalRawToDigiPortable", "eeDigis"});
  desc.add<std::string>("digisOutLabelEB", "ebDigis");
  desc.add<std::string>("digisOutLabelEE", "eeDigis");

  desc.add<bool>("produceDummyIntegrityCollections", false);

  std::string label = "ecalCPUDigisProducer";
  confDesc.add(label, desc);
}

EcalCPUDigisProducer::EcalCPUDigisProducer(const edm::ParameterSet& ps)
    : // input digi collections on host in SoA format
      digisInEBToken_{consumes(ps.getParameter<edm::InputTag>("digisInLabelEB"))},
      digisInEEToken_{consumes(ps.getParameter<edm::InputTag>("digisInLabelEE"))},

      // output digi collections in legacy format
      digisOutEBToken_{produces<EBDigiCollection>(ps.getParameter<std::string>("digisOutLabelEB"))},
      digisOutEEToken_{produces<EEDigiCollection>(ps.getParameter<std::string>("digisOutLabelEE"))},

      // whether to produce dummy integrity collections
      produceDummyIntegrityCollections_{ps.getParameter<bool>("produceDummyIntegrityCollections")},

      // dummy collections
      ebSrFlagToken_{dummyProduces<EBSrFlagCollection>()},
      eeSrFlagToken_{dummyProduces<EESrFlagCollection>()},

      // dummy integrity for xtal data
      ebIntegrityGainErrorsToken_{dummyProduces<EBDetIdCollection>("EcalIntegrityGainErrors")},
      ebIntegrityGainSwitchErrorsToken_{dummyProduces<EBDetIdCollection>("EcalIntegrityGainSwitchErrors")},
      ebIntegrityChIdErrorsToken_{dummyProduces<EBDetIdCollection>("EcalIntegrityChIdErrors")},

      // dummy integrity for xtal data - EE specific (to be rivisited towards EB+EE common collection)
      eeIntegrityGainErrorsToken_{dummyProduces<EEDetIdCollection>("EcalIntegrityGainErrors")},
      eeIntegrityGainSwitchErrorsToken_{dummyProduces<EEDetIdCollection>("EcalIntegrityGainSwitchErrors")},
      eeIntegrityChIdErrorsToken_{dummyProduces<EEDetIdCollection>("EcalIntegrityChIdErrors")},

      // dummy integrity errors
      integrityTTIdErrorsToken_{dummyProduces<EcalElectronicsIdCollection>("EcalIntegrityTTIdErrors")},
      integrityZSXtalIdErrorsToken_{dummyProduces<EcalElectronicsIdCollection>("EcalIntegrityZSXtalIdErrors")},
      integrityBlockSizeErrorsToken_{dummyProduces<EcalElectronicsIdCollection>("EcalIntegrityBlockSizeErrors")},

      //
      pnDiodeDigisToken_{dummyProduces<EcalPnDiodeDigiCollection>()},

      // dummy TCC collections
      ecalTriggerPrimitivesToken_{dummyProduces<EcalTrigPrimDigiCollection>("EcalTriggerPrimitives")},
      ecalPseudoStripInputsToken_{dummyProduces<EcalPSInputDigiCollection>("EcalPseudoStripInputs")}
{}

void EcalCPUDigisProducer::produce(edm::Event& event, edm::EventSetup const& setup) {
  // output collections
  auto digisEB = std::make_unique<EBDigiCollection>();
  auto digisEE = std::make_unique<EEDigiCollection>();

  auto const& digisEBSoAHostColl = event.get(digisInEBToken_);
  auto const& digisEESoAHostColl = event.get(digisInEEToken_);
  auto& digisEBSoAView = digisEBSoAHostColl.view();
  auto& digisEESoAView = digisEESoAHostColl.view();

  digisEB->resize(digisEBSoAView.size());
  digisEE->resize(digisEESoAView.size());

  // cast constness away
  auto* dataEB = const_cast<uint16_t*>(digisEB->data().data());
  auto* dataEE = const_cast<uint16_t*>(digisEE->data().data());
  auto* idsEB = const_cast<uint32_t*>(digisEB->ids().data());
  auto* idsEE = const_cast<uint32_t*>(digisEE->ids().data());

  // copy data
  std::memcpy(dataEB, digisEBSoAView.data()->data(), digisEBSoAView.size() * ecalPh1::sampleSize * sizeof(uint16_t));
  std::memcpy(dataEE, digisEESoAView.data()->data(), digisEESoAView.size() * ecalPh1::sampleSize * sizeof(uint16_t));
  std::memcpy(idsEB, digisEBSoAView.id(), digisEBSoAView.size() * sizeof(uint32_t));
  std::memcpy(idsEE, digisEESoAView.id(), digisEESoAView.size() * sizeof(uint32_t));

  digisEB->sort();
  digisEE->sort();

  event.put(digisOutEBToken_, std::move(digisEB));
  event.put(digisOutEEToken_, std::move(digisEE));

  if (produceDummyIntegrityCollections_) {
    // dummy collections
    event.emplace(ebSrFlagToken_);
    event.emplace(eeSrFlagToken_);
    // dummy integrity for xtal data
    event.emplace(ebIntegrityGainErrorsToken_);
    event.emplace(ebIntegrityGainSwitchErrorsToken_);
    event.emplace(ebIntegrityChIdErrorsToken_);
    // dummy integrity for xtal data - EE specific (to be rivisited towards EB+EE common collection)
    event.emplace(eeIntegrityGainErrorsToken_);
    event.emplace(eeIntegrityGainSwitchErrorsToken_);
    event.emplace(eeIntegrityChIdErrorsToken_);
    // dummy integrity errors
    event.emplace(integrityTTIdErrorsToken_);
    event.emplace(integrityZSXtalIdErrorsToken_);
    event.emplace(integrityBlockSizeErrorsToken_);
    //
    event.emplace(pnDiodeDigisToken_);
    // dummy TCC collections
    event.emplace(ecalTriggerPrimitivesToken_);
    event.emplace(ecalPseudoStripInputsToken_);
  }
}

DEFINE_FWK_MODULE(EcalCPUDigisProducer);
