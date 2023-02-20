#ifndef RecoLocalCalo_EcalRecProducers_plugins_EcalUncalibRecHitPhase2WeightsKernels_h
#define RecoLocalCalo_EcalRecProducers_plugins_EcalUncalibRecHitPhase2WeightsKernels_h

#include "DeclsForKernelsPhase2.h"

#include <alpaka/alpaka.hpp> //
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"//
#include "DataFormats/EcalDigi/interface/alpaka/EcalDigiPhase2DeviceCollection.h" //
#include "DataFormats/EcalRecHit/interface/alpaka/EcalUncalibratedRecHitDeviceCollection.h" // should the collections be also included here?



namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace ecal {
    namespace weights {

      class Phase2WeightsKernel{
        public:
          template <typename TAcc>
          ALPAKA_FN_ACC void operator()(TAcc const& acc, 
                                        double const* weightsdata,
                                        DigiPhase2DeviceCollection::ConstView digisDev,
                                        UncalibratedRecHitDeviceCollection::View recHitsDev
                                        ) const;
        };

      /*__global__ void Phase2WeightsKernel(uint16_t const* digis_in_eb,
                                          uint32_t const* dids_eb,
                                          ::ecal::reco::StorageScalarType* amplitudeEB,
                                          ::ecal::reco::StorageScalarType* amplitudeErrorEB,
                                          uint32_t* dids_outEB,
                                          int const nchannels,
                                          double const* weights_d,
                                          uint32_t* flagsEB);*/
    }  //namespace weights
  }  //namespace ecal
}  //namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif


/*questions 
shared memory ditch or keep
accessing the collection elements
passing multiple elements to a kernel(comas brackets or what)*/