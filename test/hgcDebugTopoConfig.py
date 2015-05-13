import FWCore.ParameterSet.Config as cms

process = cms.Process("HGCSimHitsAnalysis")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')    
process.load('FWCore.MessageService.MessageLogger_cfi')
#v5 geometry
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("RecoParticleFlow.PFClusterProducer.particleFlowRecHitHGCEE_cfi")



process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False)
                                        ) 

import os,sys
from UserCode.HGCanalysis.storeTools_cff import fillFromStore

fileNames = open("sample/0PU/HGG_SLHC25.txt","r")

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(fileNames),
                            skipEvents=cms.untracked.uint32(0))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

#load the analyzer
process.TFileService = cms.Service("TFileService", fileName = cms.string('debugTopo_neighbours.root'))


process.hgg = cms.EDAnalyzer("debugTopo",
                             #geometrySource   = cms.untracked.vstring('HGCalEESensitive','HGCalHESiliconSensitive',  'HGCalHEScintillatorSensitive')
                             debug = cms.uint32(0),
                             geometrySource = cms.untracked.vstring('HGCalEESensitive','HGCalHESiliconSensitive', 'HGCalHEScintillatorSensitive')
                          )


#run it
process.p = cms.Path(#process.analysis
    process.hgg

    )

