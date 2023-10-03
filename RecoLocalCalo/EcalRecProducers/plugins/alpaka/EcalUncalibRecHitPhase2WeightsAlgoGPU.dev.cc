#include "DataFormats/EcalDigi/interface/EcalDataFrame_Ph2.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

#include "EcalUncalibRecHitPhase2WeightsAlgoGPU.h"

#include "DataFormats/EcalDigi/interface/alpaka/EcalDigiPhase2DeviceCollection.h" 
#include "DataFormats/EcalRecHit/interface/alpaka/EcalUncalibratedRecHitDeviceCollection.h" 
#include <alpaka/alpaka.hpp> 
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"



namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace ecal {
    namespace weights {

      using namespace cms::alpakatools;

      class Phase2WeightsKernel{
        public:
          template <typename TAcc>
          ALPAKA_FN_ACC void operator()(TAcc const &acc, 
                                        double const *weightsdata,
                                        double const *timeWeightsdata,
                                        EcalDigiPhase2DeviceCollection::ConstView digisDev,
                                        EcalUncalibratedRecHitDeviceCollection::View recHitsDev
                                        ) const;
        };

      template <typename TAcc>
      ALPAKA_FN_ACC void Phase2WeightsKernel::operator()(TAcc const &acc, 
                                                        double const *weightsData,
                                                        double const *timeWeightsdata,                                                        
                                                        EcalDigiPhase2DeviceCollection::ConstView digisDev,
                                                        EcalUncalibratedRecHitDeviceCollection::View recHitsDev)  
                                                        const{
        
        constexpr int nsamples = EcalDataFrame_Ph2::MAXSAMPLES;                                    
        //unsigned int nchannels_per_block = alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u];    //unused currently
        auto const nchannels = digisDev.size();
        // one thread sets the output collection size scalar
        if (alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u] == 0) {
          recHitsDev.size() = digisDev.size();
        }

        auto* amplitude = recHitsDev.amplitude();                 // nchannels_per_block elements
        auto* jitter = recHitsDev.jitter();
        const auto* digis = &digisDev.data()->array;              // nchannels_per_block elements

        //unsigned int const threadx = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u];
        //unsigned int const blockx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];               //also unused currently

        const auto first = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u] +
                           alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u] *
                           alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u];

        const auto stride = alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u] * alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0u];
        for (auto tx = first; tx < nchannels; tx += stride) {
          bool g1=false;
          auto const did = DetId{digisDev.id()[tx]};
          amplitude[tx] = 0;
          jitter[tx] = 0;
          for (int sample = 0; sample < nsamples; ++sample) {
            const auto digi = digis[tx][sample];
            amplitude[tx] += ((static_cast<float>(ecalLiteDTU::adc(digi))) * 
                                 ecalPh2::gains[ecalLiteDTU::gainId(digi)] * *(weightsData + sample));
            jitter[tx] += ((static_cast<float>(ecalLiteDTU::adc(digi))) * 
                              ecalPh2::gains[ecalLiteDTU::gainId(digi)] * *(timeWeightsdata + sample));
            if (ecalLiteDTU::gainId(digi)== 1)
              g1=true;
            recHitsDev.outOfTimeAmplitudes()[tx].array[sample]= 0.;
          }
          recHitsDev.amplitudeError()[tx] = 1.0f;
          recHitsDev.id()[tx] = did.rawId();
          recHitsDev.flags()[tx] = 0;
          recHitsDev.pedestal()[tx] = 0.; 
          recHitsDev.jitterError()[tx] = 0.;
          recHitsDev.chi2()[tx] = 0.;
          recHitsDev.OOTchi2()[tx]= 0.;
          recHitsDev.aux()[tx] = 0;
          if (g1) {
            recHitsDev.flags()[tx] = 0x1 << EcalUncalibratedRecHit::kHasSwitchToGain1;
          }
        }  //if within nchannels
      }  //kernel}



      void phase2Weights(EcalDigiPhase2DeviceCollection const &digis,
                         EcalUncalibratedRecHitDeviceCollection &recHits,
                         cms::alpakatools::host_buffer<double[]> &weights,
                         cms::alpakatools::host_buffer<double[]> &timeWeights_,
                         Queue  &queue)
      {

        auto weights_d = make_device_buffer<double[]>(queue,ecalPh2::sampleSize);
        auto timeWeights_d = make_device_buffer<double[]>(queue,ecalPh2::sampleSize);
        alpaka::memcpy(queue, weights_d, weights);
        alpaka::memcpy(queue, timeWeights_d, timeWeights_);

        // use 64 items per group (this value is arbitrary, but it's a reasonable starting point)
        uint32_t items = 64;
        // use as many groups as needed to cover the whole problem
        uint32_t groups = divide_up_by(digis->metadata().size(), items);

        auto workDiv = make_workdiv<Acc1D>(groups, items);
        
        alpaka::exec<Acc1D>(queue, workDiv, Phase2WeightsKernel{}, weights_d.data(), timeWeights_d.data(), digis.const_view(),recHits.view()); 
      }

    }  // namespace weights
  }  // namespace ecal
}  //namespace ALPAKA_ACCELERATOR_NAMESPACE
