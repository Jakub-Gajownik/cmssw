#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h" 
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "DataFormats/EcalDigi/interface/alpaka/EcalDigiPhase2DeviceCollection.h" 
#include "DataFormats/EcalDigi/interface/EcalDigiPhase2HostCollection.h"




namespace ALPAKA_ACCELERATOR_NAMESPACE{
  class EcalPhase2DigiToGPUProducer : public stream::EDProducer<> {
  public:
    explicit EcalPhase2DigiToGPUProducer(edm::ParameterSet const &ps);
    ~EcalPhase2DigiToGPUProducer() override = default;
    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

    void produce(device::Event &event, device::EventSetup const &setup) override;

  private:
    const device::EDGetToken<EBDigiCollectionPh2> inputDigiToken_;
    const device::EDPutToken<ecal::DigiPhase2DeviceCollection> outputDigiDevToken_;
  };

  void EcalPhase2DigiToGPUProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
    edm::ParameterSetDescription desc;

    desc.add<edm::InputTag>("BarrelDigis", edm::InputTag("simEcalUnsuppressedDigis", ""));
    desc.add<std::string>("digisLabelEB", "ebDigis");

    descriptions.addWithDefaultLabel(desc);
}

  EcalPhase2DigiToGPUProducer::EcalPhase2DigiToGPUProducer(edm::ParameterSet const &ps)
      : inputDigiToken_(consumes(ps.getParameter<edm::InputTag>("BarrelDigis"))),
        outputDigiDevToken_(produces(ps.getParameter<std::string>("digisLabelEB"))) {}

  void EcalPhase2DigiToGPUProducer::produce(device::Event &event, device::EventSetup const &setup) {
    
  //input data from event
    const auto &inputDigis = event.get(inputDigiToken_);

    const uint32_t size = inputDigis.size();

    //create host and device collections of desired size
    ecal::DigiPhase2DeviceCollection DigisDevColl{static_cast<int32_t>(size), event.queue()};
    ecal::DigiPhase2HostCollection DigisHostColl{static_cast<int32_t>(size), event.queue()};  

  //iterate over digis
    uint32_t i = 0;
    for (const auto& inputDigi : inputDigis) {
      const int nSamples = inputDigi.size();
    //assign id to host collection
      DigisHostColl.view().id()[i] = inputDigi.id();
    //iterate over sample in digi
      for (int sample = 0; sample < nSamples; ++sample) {
      //get samples from input digi
        EcalLiteDTUSample thisSample = inputDigi[sample];
      //assign adc data to host collection
        DigisHostColl.view().data()[i * nSamples + sample][i] = thisSample.raw();
      }
      ++i;
    }

  //copy collection from host to device
    alpaka::memcpy(event.queue(), DigisDevColl.buffer(), DigisHostColl.buffer());

  //emplace device collection in the event
    event.emplace(outputDigiDevToken_, std::move(DigisDevColl));
  }
}
DEFINE_FWK_ALPAKA_MODULE(EcalPhase2DigiToGPUProducer);
