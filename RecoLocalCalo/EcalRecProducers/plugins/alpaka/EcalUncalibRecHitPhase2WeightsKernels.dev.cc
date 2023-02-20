#include <cuda.h>

#include "FWCore/Utilities/interface/CMSUnrollLoop.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalDigi/interface/EcalLiteDTUSample.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"

#include "EcalUncalibRecHitPhase2WeightsKernels.h"

#include "DataFormats/EcalDigi/interface/alpaka/EcalDigiPhase2DeviceCollection.h" //
#include "DataFormats/EcalRecHit/interface/alpaka/EcalUncalibratedRecHitDeviceCollection.h" //

#include <alpaka/alpaka.hpp>
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"






/*namespace ecal {
  namespace weights {

    __global__ void Phase2WeightsKernel(uint16_t const* digis_in,
                                        uint32_t const* __restrict__ dids,
                                        ::ecal::reco::StorageScalarType* __restrict__ amplitude,
                                        ::ecal::reco::StorageScalarType* __restrict__ amplitudeError,
                                        uint32_t* __restrict__ dids_out,
                                        int const nchannels,
                                        double const* __restrict__ weights,
                                        uint32_t* __restrict__ flags) {
      constexpr int nsamples = EcalDataFrame_Ph2::MAXSAMPLES;
      unsigned int nchannels_per_block = blockDim.x;

      // copy data from global to shared memory
      extern __shared__ char shared_mem[];
      double* shr_weights = reinterpret_cast<double*>(shared_mem);                       // nsamples elements
      float* shr_amp = reinterpret_cast<float*>(shr_weights + nsamples);                 // nchannels_per_block elements
      uint16_t* shr_digis = reinterpret_cast<uint16_t*>(shr_amp + nchannels_per_block);  // nchannels_per_block elements
      for (int i = 0; i < nsamples; ++i)
        shr_weights[i] = weights[i];

      unsigned int const threadx = threadIdx.x;
      unsigned int const blockx = blockIdx.x;

      for (int sample = 0; sample < nsamples; ++sample) {
        int const idx = threadx * nsamples + sample;
        shr_digis[idx] = digis_in[blockx * nchannels_per_block * nsamples + idx];
      }
      shr_amp[threadx] = 0.;

      __syncthreads();

      const auto first = threadIdx.x + blockIdx.x * blockDim.x;
      const auto stride = blockDim.x * gridDim.x;
      for (auto tx = first; tx < nchannels; tx += stride) {
        auto const did = DetId{dids[tx]};
        CMS_UNROLL_LOOP
        for (int sample = 0; sample < nsamples; ++sample) {
          const unsigned int idx = threadx * nsamples + sample;
          const auto shr_digi = shr_digis[idx];
          shr_amp[threadx] += (static_cast<float>(ecalLiteDTU::adc(shr_digi)) *
                               ecalPh2::gains[ecalLiteDTU::gainId(shr_digi)] * shr_weights[sample]);
        }
        const unsigned int tdx = threadx * nsamples;
        amplitude[tx] = shr_amp[threadx];
        amplitudeError[tx] = 1.0f;
        dids_out[tx] = did.rawId();
        flags[tx] = 0;
        if (ecalLiteDTU::gainId(shr_digis[tdx + nsamples - 1])) {
          flags[tx] = 0x1 << EcalUncalibratedRecHit::kHasSwitchToGain1;
        }
      }  //if within nchannels
    }    //kernel
  }      //namespace weights
}  //namespace ecal
*/

//porting
namespace ALPAKA_ACCELERATOR_NAMESPACE{
  namespace ecal {
    namespace weights {
      
      using namespace cms::alpakatools;
      template <typename TAcc>
      ALPAKA_FN_ACC void Phase2WeightsKernel::operator()(TAcc const& acc, 
                                                        double const* weightsData,                            
                                                        DigiPhase2DeviceCollection::ConstView digisDev,
                                                        UncalibratedRecHitDeviceCollection::View recHitsDev)  
                                                        const{ 

        constexpr int nsamples = EcalDataFrame_Ph2::MAXSAMPLES;                                    
        unsigned int nchannels_per_block = alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u];
        auto const nchannels = digisDev.size();

        auto* amplitude = recHitsDev.amplitude();                 // nchannels_per_block elements
        const auto* digis = digisDev.data()->data();              // nchannels_per_block elements

        unsigned int const threadx = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u];
        unsigned int const blockx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];

        const auto first = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u] +
                           alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u] *
                           alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u];

        const auto stride = alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u] * alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0u];
        for (auto tx = first; tx < nchannels; tx += stride) {
          auto const did = DetId{digisDev.id()[tx]};
          amplitude[tx] = 0;                                                          
          CMS_UNROLL_LOOP                                                                            //possibly to drop/find equiv
          for (int sample = 0; sample < nsamples; ++sample) {
            const auto digi = digis[tx][sample];
            amplitude[tx] += (static_cast<float>(ecalLiteDTU::adc(digi)) *
                                 ecalPh2::gains[ecalLiteDTU::gainId(digi)] * *(weightsData + sample));
          }
          const unsigned int tdx = threadx * nsamples;
          recHitsDev.amplitudeError()[tx] = 1.0f;
          recHitsDev.id()[tx] = did.rawId();
          recHitsDev.flags()[tx] = 0;
          if (ecalLiteDTU::gainId(digis[tx])) {
            recHitsDev.flags()[tx] = 0x1 << EcalUncalibratedRecHit::kHasSwitchToGain1;
          }
        }  //if within nchannels
      }  //kernel
    }  //namespace weights
  }  //namespace ecal
}  //namespace ALPAKA_ACCELERATOR_NAMESPACE

