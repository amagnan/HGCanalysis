import FWCore.ParameterSet.Config as cms

process = cms.Process("HGCSimHitsAnalysis")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')    
process.load('FWCore.MessageService.MessageLogger_cfi')
#v5 geometry
#process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
#v4 geometry
process.load('Configuration.Geometry.GeometryExtended2023HGCalV4MuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalV4Muon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False)
                                        #SkipEvent = cms.untracked.vstring('ProductNotFound')
                                        ) 

# configure from command line
# cmsRun test/runHGCHitsAnalyzer_cfg.py tag
# where tag can be any sub-directory under /store/cmst3/group/hgcal/CMSSW
#           or any upgrade relval sample (may need tweaking for new releases...)
ffile=0
step=-1
preFix='Single13_CMSSW_6_2_0_SLHC18'
doFullAnalysis=True
import os,sys
#if(len(sys.argv)<3):
#    print '\ncmsRun runHGCHitsAnalyzer_cfg.py doFullAnalysis tag first_file step\n'
#    print '\ttag - process tag'
#    print '\tfirst_file - first file to process'
#    print '\tstep - number of files to process\n'
#    sys.exit()

#preFix=sys.argv[2]
#if(len(sys.argv)>3):
#    if(sys.argv[3].isdigit()) : ffile=int(sys.argv[3])
#if(len(sys.argv)>4):
#    if(sys.argv[4].isdigit()) : step=int(sys.argv[4])
#print '[runHGCHitsAnalyzer] processing %d files of %s, starting from %d'%(step,preFix,ffile)

#configure the source (list all files in directory within range [ffile,ffile+step[
from UserCode.HGCanalysis.storeTools_cff import fillFromStore
#process.source = cms.Source("PoolSource",                            
#                            fileNames=cms.untracked.vstring()
#                            )
#if preFix.find('/store')>=0 :
#    process.source.fileNames=fillFromStore(preFix,ffile,step)
#else :
#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("file:/afs/cern.ch/user/l/lcorpe/work/public/HGCAL/SingleElectronPt35_PU0_RECO_1.root",
#"file:/afs/cern.ch/user/l/lcorpe/work/public/HGCAL/SingleElectronPt35_PU0_RECO_2.root",
#"file:/afs/cern.ch/user/l/lcorpe/work/public/HGCAL/SingleElectronPt35_PU0_RECO_3.root",
#"file:/afs/cern.ch/user/l/lcorpe/work/public/HGCAL/SingleElectronPt35_PU0_RECO_4.root",
#"file:/afs/cern.ch/user/l/lcorpe/work/public/HGCAL/SingleElectronPt35_PU0_RECO_5.root",
#"file:/afs/cern.ch/user/l/lcorpe/work/public/HGCAL/SingleElectronPt35_PU0_RECO_6.root",
#"file:/afs/cern.ch/user/l/lcorpe/work/public/HGCAL/SingleElectronPt35_PU0_RECO_7.root",
#"file:/afs/cern.ch/user/l/lcorpe/work/public/HGCAL/SingleElectronPt35_PU0_RECO_8.root",
#"file:/afs/cern.ch/user/l/lcorpe/work/public/HGCAL/SingleElectronPt35_PU0_RECO_9.root"))

#fileNames = open("LCFilenames.txt","r")
fileNames = open("newRecoFiles3.txt","r")

process.source = cms.Source("PoolSource",
                            #fileNames=cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/relval/CMSSW_6_2_0_SLHC22/RelValH130GGgluonfusion_14TeV/GEN-SIM-RECO/PH2_1K_FB_V6_UPGHGCalV5-v1/00000/1CC2630B-6A8F-E411-95D3-0025905A48BA.root"),
                            #fileNames=cms.untracked.vstring(fileNames),
                            fileNames=cms.untracked.vstring("file:../../../sample/0PU/Hgg0PU-2kEvents_0_KEEP.root"),
                            #fileNames=cms.untracked.vstring("file:/afs/cern.ch/user/l/lcorpe/work/private/HGCALreco3/CMSSW_6_2_0_SLHC22/src/Hgg0PU-1kEvents_1.root"),
                            skipEvents=cms.untracked.uint32(0))

#process.source.fileNames=fillFromStore('/store/cmst3/group/hgcal/CMSSW/%s'%preFix,ffile,step)
#process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#load the analyzer
import getpass
whoami=getpass.getuser()
outputTag=preFix.replace('/','_')
#process.TFileService = cms.Service("TFileService", fileName = cms.string('/tmp/%s/%s_Hits_%d.root'%(whoami,outputTag,ffile)))
process.TFileService = cms.Service("TFileService", fileName = cms.string('BasicHggPatch2.root'))
process.load('UserCode.HGCanalysis.hgcHitsAnalyzer_cfi')

process.hgg = cms.EDAnalyzer("BasicHggAnalyser",
                          #geometrySource   = cms.untracked.vstring('HGCalEESensitive','HGCalHESiliconSensitive',  'HGCalHEScintillatorSensitive')
												endcapRecHitCollection = cms.untracked.InputTag("HGCalRecHit:HGCEERecHits"),
												endcapSuperClusterCollection = cms.untracked.InputTag("particleFlowSuperClusterHGCEE"),
												endcapClusterCollection = cms.untracked.InputTag("particleFlowClusterHGCEE"),
												genParticlesTag =  cms.untracked.InputTag("genParticles"),
                          )


#run it
process.p = cms.Path(#process.analysis
										 process.hgg
)

