#alpaka with no switch producer
import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Accelerators_cff import *
from HeterogeneousCore.AlpakaCore.ProcessAcceleratorAlpaka_cfi import ProcessAcceleratorAlpaka

from RecoLocalCalo.EcalRecProducers.ecalPhase2DigiToPortableProducer_cfi import ecalPhase2DigiToPortableProducer as _ecalPhase2DigiToPortableProducer
ecalPhase2DigiToPortableProducer = _ecalPhase2DigiToPortableProducer.clone()

# portable weights
from RecoLocalCalo.EcalRecProducers.ecalUncalibRecHitPhase2Portable_cfi import ecalUncalibRecHitPhase2Portable as _ecalUncalibRecHitPhase2Portable
ecalUncalibRecHitPhase2Portable = _ecalUncalibRecHitPhase2Portable.clone(
  digisLabelEB = ('ecalPhase2DigiToPortableProducer', 'ebDigis')
)

from RecoLocalCalo.EcalRecProducers.ecalUncalibRecHitConvertPortable2CPUFormat_cfi import ecalUncalibRecHitConvertPortable2CPUFormat as _ecalUncalibRecHitConvertPortable2CPUFormat
ecalUncalibRecHitPhase2 = _ecalUncalibRecHitConvertPortable2CPUFormat.clone(
    isPhase2 = True,
    recHitsLabelPortableEB = ('ecalUncalibRecHitPhase2Portable', 'EcalUncalibRecHitsEB'),
    recHitsLabelPortableEE = None,  # remove unneeded Phase1 parameters
    recHitsLabelCPUEE = None
)


ecalUncalibRecHitPhase2Task = cms.Task(
  # convert phase2 digis to Portable Collection
  ecalPhase2DigiToPortableProducer, 
  # ECAL weights running on Portable
  ecalUncalibRecHitPhase2Portable,
  # ECAL multifit running on CPU, or convert the uncalibrated rechits from Portable Collection to legacy format
  ecalUncalibRecHitPhase2
)

#previous CUDA switch producer
'''import FWCore.ParameterSet.Config as cms
from HeterogeneousCore.CUDACore.SwitchProducerCUDA import SwitchProducerCUDA
from Configuration.ProcessModifiers.gpu_cff import gpu


from RecoLocalCalo.EcalRecProducers.ecalUncalibRecHitPhase2_cfi import ecalUncalibRecHitPhase2 as _ecalUncalibRecHitPhase2
ecalUncalibRecHitPhase2 = SwitchProducerCUDA(
          cpu = _ecalUncalibRecHitPhase2.clone()
          )

# cpu weights
ecalUncalibRecHitPhase2Task = cms.Task(ecalUncalibRecHitPhase2)

# conditions used on gpu


from Configuration.StandardSequences.Accelerators_cff import *
from HeterogeneousCore.AlpakaCore.ProcessAcceleratorAlpaka_cfi import ProcessAcceleratorAlpaka

from RecoLocalCalo.EcalRecProducers.ecalPhase2DigiToGPUProducer_cfi import ecalPhase2DigiToGPUProducer as _ecalPhase2DigiToGPUProducer
ecalPhase2DigiToGPUProducer = _ecalPhase2DigiToGPUProducer.clone()

# gpu weights
from RecoLocalCalo.EcalRecProducers.ecalUncalibRecHitPhase2GPU_cfi import ecalUncalibRecHitPhase2GPU as _ecalUncalibRecHitPhase2GPU
ecalUncalibRecHitPhase2GPU = _ecalUncalibRecHitPhase2GPU.clone(
          digisLabelEB = ('ecalPhase2DigiToGPUProducer', 'ebDigis')
          )

from RecoLocalCalo.EcalRecProducers.ecalUncalibRecHitConvertPortable2CPUFormat_cfi import ecalUncalibRecHitConvertPortable2CPUFormat as _ecalUncalibRecHitConvertPortable2CPUFormat
gpu.toModify(ecalUncalibRecHitPhase2,
            cuda = _ecalUncalibRecHitConvertPortable2CPUFormat.clone(
                        isPhase2 = True,
                        recHitsLabelGPUEB = ('ecalUncalibRecHitPhase2GPU', 'EcalUncalibRecHitsEB'),
                        recHitsLabelGPUEE = None,  # remove unneeded Phase1 parameters
                        recHitsLabelCPUEE = None
                )
            )

gpu.toReplaceWith(ecalUncalibRecHitPhase2Task, cms.Task(
# convert phase2 digis to GPU SoA
ecalPhase2DigiToGPUProducer, 
# ECAL weights running on GPU
ecalUncalibRecHitPhase2GPU,
# ECAL multifit running on CPU, or convert the uncalibrated rechits from SoA to legacy format
ecalUncalibRecHitPhase2,
))
'''
