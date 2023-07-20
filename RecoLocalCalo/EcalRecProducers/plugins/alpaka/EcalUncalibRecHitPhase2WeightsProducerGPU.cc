#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h" 
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"


#include "DataFormats/EcalDigi/interface/EcalDataFrame_Ph2.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"
#include "DataFormats/EcalDigi/interface/alpaka/EcalDigiPhase2DeviceCollection.h" 
#include "DataFormats/EcalDigi/interface/EcalDigiPhase2HostCollection.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHitHostCollection.h"
#include "DataFormats/EcalRecHit/interface/alpaka/EcalUncalibratedRecHitDeviceCollection.h"
#include "DataFormats/Portable/interface/Product.h" 

#include "EcalUncalibRecHitPhase2WeightsAlgoGPU.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE{
  class EcalUncalibRecHitPhase2WeightsProducerGPU : public stream::EDProducer<> {
  public:
    explicit EcalUncalibRecHitPhase2WeightsProducerGPU(edm::ParameterSet const &ps);
    ~EcalUncalibRecHitPhase2WeightsProducerGPU() override = default;
    static void fillDescriptions(edm::ConfigurationDescriptions &);

    void produce(device::Event &, device::EventSetup const &) override;

  private:
    cms::alpakatools::host_buffer<double[]> weights_;   
    
    using InputProduct = ecal::DigiPhase2DeviceCollection;
    const device::EDGetToken<InputProduct> digisToken_;           //assumed both will be stored on device side    
    using OutputProduct = ecal::UncalibratedRecHitDeviceCollection; 
    const device::EDPutToken<OutputProduct> recHitsToken_;
  };

  // constructor with initialisation of elements
  EcalUncalibRecHitPhase2WeightsProducerGPU::EcalUncalibRecHitPhase2WeightsProducerGPU(const edm::ParameterSet &ps) 
      : weights_{cms::alpakatools::make_host_buffer<double[]>(ecalPh2::sampleSize)},
        digisToken_{consumes(ps.getParameter<edm::InputTag>("digisLabelEB"))},
        recHitsToken_{produces(ps.getParameter<std::string>("recHitsLabelEB"))} {
    //extracting the weights, for-loop to save them to the buffer even if the size is different than standard
    const auto weights = ps.getParameter<std::vector<double>>("weights");
    for (unsigned int i = 0; i < ecalPh2::sampleSize; ++i) {
      if (i < weights.size()){
        weights_[i] = weights[i];
      } else {
        weights_[i]= 0;
      }
    }
  }
  
  void EcalUncalibRecHitPhase2WeightsProducerGPU::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
    edm::ParameterSetDescription desc;
  
    desc.add<std::string>("recHitsLabelEB", "EcalUncalibRecHitsEB");
    //The below weights values should be kept up to date with those on the CPU version of this module
    desc.add<std::vector<double>>("weights",                 
                                  {-0.121016,
                                   -0.119899,
                                   -0.120923,
                                   -0.0848959,
                                   0.261041,
                                   0.509881,
                                   0.373591,
                                   0.134899,
                                   -0.0233605,
                                   -0.0913195,
                                   -0.112452,
                                   -0.118596,
                                   -0.121737,
                                   -0.121737,
                                   -0.121737,
                                   -0.121737});
  
    desc.add<edm::InputTag>("digisLabelEB", edm::InputTag("simEcalUnsuppressedDigis", ""));
  
    descriptions.addWithDefaultLabel(desc);
  }

  void EcalUncalibRecHitPhase2WeightsProducerGPU::produce(device::Event &event, const device::EventSetup &setup) {

    // device collection of digis products
    auto const &digis = event.get(digisToken_);

    //get size of digis
    const uint32_t size = digis->metadata().size();

    //allocate product of queue on the device
    OutputProduct recHits{static_cast<int32_t>(size), event.queue()};
  
    // do not run the algo if there are no digis
    if (size > 0) {
    //launch the asynchronous work
    ecal::weights::phase2Weights(digis, recHits, weights_, event.queue());
    }
    // put into the event
    event.emplace(recHitsToken_, std::move(recHits));
  }

} //namespace ALPAKA_ACCELERATOR_NAMESPACE
DEFINE_FWK_ALPAKA_MODULE(EcalUncalibRecHitPhase2WeightsProducerGPU);
