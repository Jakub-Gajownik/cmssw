import FWCore.ParameterSet.Config as cms
from HeterogeneousCore.CUDACore.SwitchProducerCUDA import SwitchProducerCUDA
from Configuration.ProcessModifiers.gpu_cff import gpu

# ECAL multifit running on CPU
from RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi import ecalMultiFitUncalibRecHit as _ecalMultiFitUncalibRecHit
ecalMultiFitUncalibRecHit = SwitchProducerCUDA(
  cpu = _ecalMultiFitUncalibRecHit.clone()
)

ecalMultiFitUncalibRecHitTask = cms.Task(
  # ECAL multifit running on CPU
  ecalMultiFitUncalibRecHit
)

from Configuration.StandardSequences.Accelerators_cff import *
from HeterogeneousCore.AlpakaCore.ProcessAcceleratorAlpaka_cfi import ProcessAcceleratorAlpaka

# ECAL conditions used by the multifit running on GPU
from RecoLocalCalo.EcalRecProducers.ecalMultifitConditionsPortableESProducer_cfi import ecalMultifitConditionsPortableESProducer
from RecoLocalCalo.EcalRecProducers.ecalMultifitParametersGPUESProducer_cfi import ecalMultifitParametersGPUESProducer

# ECAL multifit running on GPU
from RecoLocalCalo.EcalRecProducers.ecalUncalibRecHitProducerPortable_cfi import ecalUncalibRecHitProducerPortable as _ecalUncalibRecHitProducerPortable
ecalMultiFitUncalibRecHitPortable = _ecalUncalibRecHitProducerPortable.clone(
  digisLabelEB = cms.InputTag('ecalDigisPortable', 'ebDigis'),
  digisLabelEE = cms.InputTag('ecalDigisPortable', 'eeDigis'),
)

# convert the uncalibrated rechits from SoA to legacy format
from RecoLocalCalo.EcalRecProducers.ecalUncalibRecHitConvertPortable2CPUFormat_cfi import ecalUncalibRecHitConvertPortable2CPUFormat as _ecalUncalibRecHitConvertPortable2CPUFormat
gpu.toModify(ecalMultiFitUncalibRecHit,
  cuda = _ecalUncalibRecHitConvertPortable2CPUFormat.clone(
    recHitsLabelGPUEB = cms.InputTag('ecalMultiFitUncalibRecHitPortable', 'EcalUncalibRecHitsEB'),
    recHitsLabelGPUEE = cms.InputTag('ecalMultiFitUncalibRecHitPortable', 'EcalUncalibRecHitsEE'),
  )
)

gpu.toReplaceWith(ecalMultiFitUncalibRecHitTask, cms.Task(
  # ECAL conditions used by the multifit running on GPU
  ecalMultifitConditionsPortableESProducer,
  ecalMultifitParametersGPUESProducer,
  # ECAL multifit running on device
  ecalMultiFitUncalibRecHitPortable,
  # ECAL multifit running on CPU, or convert the uncalibrated rechits from SoA to legacy format
  ecalMultiFitUncalibRecHit,
))
