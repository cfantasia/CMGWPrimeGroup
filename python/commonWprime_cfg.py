import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
from UserCode.CMGWPrimeGroup.selectors_cff import *
from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import heepBarrelCuts, heepEndcapCuts

process = cms.Process("WPrimeAnalysis")
# get JSON file correctly parsed
JSONfile = 'UserCode/CMGWPrimeGroup/JSON/Cert_160404-167784_7TeV_PromptReco_Collisions11_JSON.txt'
myList = LumiList.LumiList (filename = JSONfile).getCMSSWString().split(',')

process.inputs = cms.PSet (
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    )
process.inputs.lumisToProcess.extend(myList)


process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring()
)

process.MessageLogger = cms.Service("MessageLogger")

process.WprimeAnalyzer = cms.PSet(
    ## common input for wrapped analyzers
    fileNames   = cms.vstring(),  ## keep empty!

    sample_cross_sections = cms.string("samples_cross_sections.txt"),
    outputFile  = cms.string('Wprime_analysis.root'),
    logFile = cms.string("Wprime_event_counts.txt"),

    maxEvents   = cms.int32(-1),                     
    reportAfter = cms.uint32(5000),                  

    useJSON = cms.bool(True),
    countGenEvts = cms.bool(False),
    eventCounters = cms.vstring(),

    doRecoilCorrectionForW = cms.bool(False),
    useAdjustedMET = cms.bool(False),

    ## enable analysis in individual channels
    runMuMETAnalysis = cms.bool(False),
    runElMETAnalysis = cms.bool(False),
    runWZAnalysis    = cms.bool(False),
    runHadVZAnalysis = cms.bool(False),
    runTBAnalysis    = cms.bool(False),
    runWgammaAnalysis = cms.bool(False),

    ## input specific for this analyzer
    muonReconstructor = cms.int32(3), ## see TeVMuon_tracking.h

    muons = cms.string('userPatMuons'),
    electrons = cms.string('userPatElectrons'),
    met   = cms.string('patMETsPF'),
    particleFlow = cms.InputTag('selectedPatPFParticles'),
    genParticles = cms.InputTag('prunedGenParticles'),
    hltEventTag = cms.string('patTriggerEvent'),
    pileupTag  = cms.string('addPileupInfo'),

    inputs = process.inputs,

    #PileUp Inputs
    MCPUDistFile = cms.string('UserCode/CMGWPrimeGroup/root_macros/MCPUDist.root'),
    MCPUDistHist = cms.string('pileup'),
    DataPUDistFile = cms.string('UserCode/CMGWPrimeGroup/root_macros/DataPUDist.root'),
    DataPUDistHist = cms.string('pileup'),
    
    )

