# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

process = cms.Process("HESourcing")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

import sys
runNumber = sys.argv[2]
runType = int(sys.argv[3])
irunNumber = int(runNumber)
emap = "emap_v1.txt"
if irunNumber >= 287582:
	emap = "emap_v2.txt"
print emap

process.source = cms.Source("HcalTBSource",
    fileNames = cms.untracked.vstring('file:./Data/USC_'+runNumber+'.root')
)

process.hcalAnalyzer = cms.EDAnalyzer('HESourcing',
        OutFileName = cms.untracked.string('N_'+runNumber+'.root'),
	RunType = cms.int32(runType),
	histoFED =  cms.int32(65),
	driverFED = cms.int32(12)
)

process.load('Configuration.Geometry.GeometryIdeal_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.autoCond import autoCond

from CondCore.CondDB.CondDB_cfi import *

process.GlobalTag.globaltag = autoCond['startup'] 

process.es_ascii = cms.ESSource('HcalTextCalibrations',
        input = cms.VPSet(
               cms.PSet(
                object = cms.string('ElectronicsMap'),
                file = cms.FileInPath('HCALCommissioning2017/HESourcing/test/'+emap)
               )
        )
)

process.es_prefer = cms.ESPrefer('HcalTextCalibrations', 'es_ascii')
process.p = cms.Path(process.hcalAnalyzer)

