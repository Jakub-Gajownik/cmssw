import FWCore.ParameterSet.Config as cms

from DQM.EcalMonitorTasks.EcalMonitorTask_cfi import *

# Customization to run the CPU vs GPU comparison task if the job runs on a GPU enabled machine
from Configuration.ProcessModifiers.gpuValidationEcal_cff import gpuValidationEcal
from DQM.EcalMonitorTasks.ecalGpuTask_cfi import ecalGpuTask

gpuValidationEcal.toModify(ecalGpuTask.params, runGpuTask = True)
gpuValidationEcal.toModify(ecalMonitorTask.workers, func = lambda workers: workers.append("GpuTask"))
gpuValidationEcal.toModify(ecalMonitorTask, workerParameters = dict(GpuTask = ecalGpuTask))

# customisation for Phase 2
from Configuration.Eras.Modifier_phase2_ecal_devel_cff import phase2_ecal_devel
phase2_ecal_devel.toModify(ecalGpuTask.params, enableDigi = False)
phase2_ecal_devel.toModify(ecalGpuTask.params, enableRecHit = False)
phase2_ecal_devel.toModify(ecalMonitorTask.collectionTags, EBCpuUncalibRecHit = "ecalUncalibRecHitPhase2@cpu:EcalUncalibRecHitsEB")
phase2_ecal_devel.toModify(ecalMonitorTask.collectionTags, EBGpuUncalibRecHit = "ecalUncalibRecHitPhase2@cuda:EcalUncalibRecHitsEB")
