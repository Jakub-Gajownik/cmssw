#include "DataFormats/EcalDigi/interface/EcalDataFrame_Ph2.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

#include "EcalUncalibRecHitPhase2WeightsKernels.h"
#include "EcalUncalibRecHitPhase2WeightsAlgoGPU.h"


#include "DataFormats/EcalDigi/interface/alpaka/EcalDigiPhase2DeviceCollection.h" //
#include "DataFormats/EcalRecHit/interface/alpaka/EcalUncalibratedRecHitDeviceCollection.h" //
#include <alpaka/alpaka.hpp> //
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"//
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"//
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"//


/*
namespace ecal {
  namespace weights {

    void phase2Weights(ecal::DigisCollection<calo::common::DevStoragePolicy> const& digis,
                       EventOutputDataGPU& eventOutputGPU,
                       cms::cuda::device::unique_ptr<double[]>& weights_d,
                       cudaStream_t cudaStream) {
      unsigned int const totalChannels = digis.size;
      // 64 threads per block best occupancy from Nsight compute profiler
      unsigned int const threads_1d = 64;
      unsigned int const blocks_1d = (totalChannels + threads_1d - 1) / threads_1d;
      int shared_bytes = EcalDataFrame_Ph2::MAXSAMPLES * sizeof(double) +
                         threads_1d * (EcalDataFrame_Ph2::MAXSAMPLES * (sizeof(uint16_t)) + sizeof(float));
      Phase2WeightsKernel<<<blocks_1d, threads_1d, shared_bytes, cudaStream>>>(
          digis.data.get(),
          digis.ids.get(),
          eventOutputGPU.recHits.amplitude.get(),
          eventOutputGPU.recHits.amplitudeError.get(),
          eventOutputGPU.recHits.did.get(),
          totalChannels,
          weights_d.get(),
          eventOutputGPU.recHits.flags.get());
      cudaCheck(cudaGetLastError());
    }

  }  // namespace weights
}  // namespace ecal
*/

/* size of the digis, chceck alpaka for blocks threads, pass the collections to the kernel, 
comment off shared bytes for now*/

//porting
namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace ecal {
    namespace weights {

      using namespace cms::alpakatools;
      void phase2Weights(ecal::DigiPhase2DeviceCollection const& digis,
                         ecal::UncalibratedRecHitDeviceCollection & recHits,
                         cms::alpakatools::host_buffer<double[]> & weights_,
                         Queue  &queue)
      {

        auto weights_d = make_device_buffer<double[]>(queue,ecalPh2::sampleSize);  
        alpaka::memcpy(queue, weights_d, weights_);

        // use 64 items per group (this value is arbitrary, but it's a reasonable starting point)
        uint32_t items = 64;
        // use as many groups as needed to cover the whole problem
        uint32_t groups = divide_up_by(digis.const_view().size(), items);

        auto workDiv = make_workdiv<Acc1D>(groups, items); //what are groups and items here?
        
        alpaka::exec<Acc1D>(queue, workDiv, Phase2WeightsKernel{}, weights_d.data(), digis.const_view(),recHits.view()); 

        /*Phase2WeightsKernel<<<blocks_1d, threads_1d, shared_bytes, cudaStream>>>(
            digis.data.get(),
            digis.ids.get(),
            eventOutputGPU.recHits.amplitude.get(),
            eventOutputGPU.recHits.amplitudeError.get(),
            eventOutputGPU.recHits.did.get(),
            totalChannels,
            weights_d.get(),
            eventOutputGPU.recHits.flags.get());
        cudaCheck(cudaGetLastError());*/
      }

    }  // namespace weights
  }  // namespace ecal
}  //namespace ALPAKA_ACCELERATOR_NAMESPACE