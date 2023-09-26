#include "CondFormats/DataRecord/interface/EcalMultifitConditionsRcd.h"
#include "CondFormats/EcalObjects/interface/alpaka/EcalMultifitConditionsPortable.h"
#include "DataFormats/EcalDigi/interface/alpaka/EcalDigiDeviceCollection.h"
#include "DataFormats/EcalRecHit/interface/alpaka/EcalUncalibratedRecHitDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/global/EDProducer.h"

#include "DeclsForKernels.h"
#include "EcalUncalibRecHitMultiFitAlgoPortable.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class EcalUncalibRecHitProducerPortable : public global::EDProducer<> {
  public:
    explicit EcalUncalibRecHitProducerPortable(edm::ParameterSet const& ps);
    ~EcalUncalibRecHitProducerPortable() override = default;
    static void fillDescriptions(edm::ConfigurationDescriptions&);
  
    void produce(edm::StreamID, device::Event&, device::EventSetup const&) const override;
  
  private:
    using InputProduct = ecal::DigiDeviceCollection;
    const device::EDGetToken<InputProduct> digisTokenEB_;
    const device::EDGetToken<InputProduct> digisTokenEE_;
    using OutputProduct = ecal::UncalibratedRecHitDeviceCollection;
    const device::EDPutToken<OutputProduct> uncalibRecHitsTokenEB_;
    const device::EDPutToken<OutputProduct> uncalibRecHitsTokenEE_;
  
    // conditions tokens
    const device::ESGetToken<EcalMultifitConditionsPortableDevice, EcalMultifitConditionsRcd> multifitConditionsToken_;
    //const edm::ESGetToken<EcalMultifitParametersGPU, JobConfigurationGPURecord> multifitParametersToken_;
  
    // configuration parameters
    ecal::multifit::ConfigurationParameters configParameters_;
  };
  
  void EcalUncalibRecHitProducerPortable::fillDescriptions(edm::ConfigurationDescriptions& confDesc) {
    edm::ParameterSetDescription desc;
  
    desc.add<edm::InputTag>("digisLabelEB", edm::InputTag("ecalRawToDigiPortable", "ebDigis"));
    desc.add<edm::InputTag>("digisLabelEE", edm::InputTag("ecalRawToDigiPortable", "eeDigis"));
  
    desc.add<std::string>("recHitsLabelEB", "EcalUncalibRecHitsEB");
    desc.add<std::string>("recHitsLabelEE", "EcalUncalibRecHitsEE");
  
    desc.add<double>("EBtimeFitLimits_Lower", 0.2);
    desc.add<double>("EBtimeFitLimits_Upper", 1.4);
    desc.add<double>("EEtimeFitLimits_Lower", 0.2);
    desc.add<double>("EEtimeFitLimits_Upper", 1.4);
    desc.add<double>("EBtimeConstantTerm", .6);
    desc.add<double>("EEtimeConstantTerm", 1.0);
    desc.add<double>("EBtimeNconst", 28.5);
    desc.add<double>("EEtimeNconst", 31.8);
    desc.add<double>("outOfTimeThresholdGain12pEB", 5);
    desc.add<double>("outOfTimeThresholdGain12mEB", 5);
    desc.add<double>("outOfTimeThresholdGain61pEB", 5);
    desc.add<double>("outOfTimeThresholdGain61mEB", 5);
    desc.add<double>("outOfTimeThresholdGain12pEE", 1000);
    desc.add<double>("outOfTimeThresholdGain12mEE", 1000);
    desc.add<double>("outOfTimeThresholdGain61pEE", 1000);
    desc.add<double>("outOfTimeThresholdGain61mEE", 1000);
    desc.add<double>("amplitudeThresholdEB", 10);
    desc.add<double>("amplitudeThresholdEE", 10);
    desc.addUntracked<std::vector<uint32_t>>("kernelMinimizeThreads", {32, 1, 1});
    desc.add<bool>("shouldRunTimingComputation", true);
    confDesc.addWithDefaultLabel(desc);
  }
  
  EcalUncalibRecHitProducerPortable::EcalUncalibRecHitProducerPortable(const edm::ParameterSet& ps)
      : digisTokenEB_{consumes(ps.getParameter<edm::InputTag>("digisLabelEB"))},
        digisTokenEE_{consumes(ps.getParameter<edm::InputTag>("digisLabelEE"))},
        uncalibRecHitsTokenEB_{produces(ps.getParameter<std::string>("recHitsLabelEB"))},
        uncalibRecHitsTokenEE_{produces(ps.getParameter<std::string>("recHitsLabelEE"))},
        multifitConditionsToken_{esConsumes()} {
    std::pair<double, double> EBtimeFitLimits, EEtimeFitLimits;
    EBtimeFitLimits.first = ps.getParameter<double>("EBtimeFitLimits_Lower");
    EBtimeFitLimits.second = ps.getParameter<double>("EBtimeFitLimits_Upper");
    EEtimeFitLimits.first = ps.getParameter<double>("EEtimeFitLimits_Lower");
    EEtimeFitLimits.second = ps.getParameter<double>("EEtimeFitLimits_Upper");
  
    auto EBtimeConstantTerm = ps.getParameter<double>("EBtimeConstantTerm");
    auto EEtimeConstantTerm = ps.getParameter<double>("EEtimeConstantTerm");
    auto EBtimeNconst = ps.getParameter<double>("EBtimeNconst");
    auto EEtimeNconst = ps.getParameter<double>("EEtimeNconst");
  
    auto outOfTimeThreshG12pEB = ps.getParameter<double>("outOfTimeThresholdGain12pEB");
    auto outOfTimeThreshG12mEB = ps.getParameter<double>("outOfTimeThresholdGain12mEB");
    auto outOfTimeThreshG61pEB = ps.getParameter<double>("outOfTimeThresholdGain61pEB");
    auto outOfTimeThreshG61mEB = ps.getParameter<double>("outOfTimeThresholdGain61mEB");
    auto outOfTimeThreshG12pEE = ps.getParameter<double>("outOfTimeThresholdGain12pEE");
    auto outOfTimeThreshG12mEE = ps.getParameter<double>("outOfTimeThresholdGain12mEE");
    auto outOfTimeThreshG61pEE = ps.getParameter<double>("outOfTimeThresholdGain61pEE");
    auto outOfTimeThreshG61mEE = ps.getParameter<double>("outOfTimeThresholdGain61mEE");
    auto amplitudeThreshEB = ps.getParameter<double>("amplitudeThresholdEB");
    auto amplitudeThreshEE = ps.getParameter<double>("amplitudeThresholdEE");
  
    // switch to run timing computation kernels
    configParameters_.shouldRunTimingComputation = ps.getParameter<bool>("shouldRunTimingComputation");
  
    // minimize kernel launch conf
    auto threadsMinimize = ps.getUntrackedParameter<std::vector<uint32_t>>("kernelMinimizeThreads");
    configParameters_.kernelMinimizeThreads[0] = threadsMinimize[0];
    configParameters_.kernelMinimizeThreads[1] = threadsMinimize[1];
    configParameters_.kernelMinimizeThreads[2] = threadsMinimize[2];
  
    //
    // configuration and physics parameters: done once
    // assume there is a single device
    // use sync copying
    //
  
    // time fit parameters and limits
    configParameters_.timeFitLimitsFirstEB = EBtimeFitLimits.first;
    configParameters_.timeFitLimitsSecondEB = EBtimeFitLimits.second;
    configParameters_.timeFitLimitsFirstEE = EEtimeFitLimits.first;
    configParameters_.timeFitLimitsSecondEE = EEtimeFitLimits.second;
  
    // time constant terms
    configParameters_.timeConstantTermEB = EBtimeConstantTerm;
    configParameters_.timeConstantTermEE = EEtimeConstantTerm;
  
    // time N const
    configParameters_.timeNconstEB = EBtimeNconst;
    configParameters_.timeNconstEE = EEtimeNconst;
  
    // amplitude threshold for time flags
    configParameters_.amplitudeThreshEB = amplitudeThreshEB;
    configParameters_.amplitudeThreshEE = amplitudeThreshEE;
  
    // out of time thresholds gain-dependent
    configParameters_.outOfTimeThreshG12pEB = outOfTimeThreshG12pEB;
    configParameters_.outOfTimeThreshG12pEE = outOfTimeThreshG12pEE;
    configParameters_.outOfTimeThreshG61pEB = outOfTimeThreshG61pEB;
    configParameters_.outOfTimeThreshG61pEE = outOfTimeThreshG61pEE;
    configParameters_.outOfTimeThreshG12mEB = outOfTimeThreshG12mEB;
    configParameters_.outOfTimeThreshG12mEE = outOfTimeThreshG12mEE;
    configParameters_.outOfTimeThreshG61mEB = outOfTimeThreshG61mEB;
    configParameters_.outOfTimeThreshG61mEE = outOfTimeThreshG61mEE;
  }
  
  void EcalUncalibRecHitProducerPortable::produce(edm::StreamID sid, device::Event& event, device::EventSetup const& setup) const {
    //DurationMeasurer<std::chrono::milliseconds> timer{std::string{"produce duration"}};
  
    // get device collections from event
    auto const& ebDigisDev = event.get(digisTokenEB_);
    auto const& eeDigisDev = event.get(digisTokenEE_);

    // FIXME would be better to have the actual number of digis here
    const auto neb = ebDigisDev.const_view().metadata().size();
    const auto nee = eeDigisDev.const_view().metadata().size();
  
    // output device collections
    OutputProduct uncalibRecHitsDevEB{neb, event.queue()};
    OutputProduct uncalibRecHitsDevEE{nee, event.queue()};
    // reset the size scalar of the SoA
    // memset takes an alpaka view that is created from the scalar in a view to the portable device collection
    auto uncalibRecHitSizeViewEB = cms::alpakatools::make_device_view<uint32_t>(alpaka::getDev(event.queue()), uncalibRecHitsDevEB.view().size());
    auto uncalibRecHitSizeViewEE = cms::alpakatools::make_device_view<uint32_t>(alpaka::getDev(event.queue()), uncalibRecHitsDevEE.view().size());
    alpaka::memset(event.queue(), uncalibRecHitSizeViewEB, 0);
    alpaka::memset(event.queue(), uncalibRecHitSizeViewEE, 0);

    // stop here if there are no digis
    if (neb + nee > 0) {
      // conditions
      auto const& multifitConditionsDev = setup.getData(multifitConditionsToken_);
      //auto const& multifitParameters = setup.getData(multifitParametersToken_);
  
      //// assign ptrs/values: this is done not to change how things look downstream
      //configParameters_.amplitudeFitParametersEB = multifitParameters.amplitudeFitParametersEB.get();
      //configParameters_.amplitudeFitParametersEE = multifitParameters.amplitudeFitParametersEE.get();
      //configParameters_.timeFitParametersEB = multifitParameters.timeFitParametersEB.get();
      //configParameters_.timeFitParametersEE = multifitParameters.timeFitParametersEE.get();
      //configParameters_.timeFitParametersSizeEB = multifitParametersData.getValues()[2].get().size();
      //configParameters_.timeFitParametersSizeEE = multifitParametersData.getValues()[3].get().size();
  
      //
      // schedule algorithms
      //
      ecal::multifit::entryPoint(
          ebDigisDev, eeDigisDev, uncalibRecHitsDevEB, uncalibRecHitsDevEE,
          multifitConditionsDev,
          configParameters_, event.queue());
    }
  
    // put into the event
    event.emplace(uncalibRecHitsTokenEB_, std::move(uncalibRecHitsDevEB));
    event.emplace(uncalibRecHitsTokenEE_, std::move(uncalibRecHitsDevEE));
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(EcalUncalibRecHitProducerPortable);
