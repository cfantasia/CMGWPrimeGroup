import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
from UserCode.CMGWPrimeGroup.selectors_cff import *
from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import heepBarrelCuts, heepEndcapCuts
from SHarper.HEEPAnalyzer.WP80SelectionCuts_cfi import wp80BarrelCuts, wp80EndcapCuts

process = cms.Process("WPrimeAnalysis")
# get JSON file correctly parsed
goldenJSONfile = 'UserCode/CMGWPrimeGroup/JSON/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt'
MuonPhysJSONfile = 'UserCode/CMGWPrimeGroup/JSON/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt'
goldenJSONList = LumiList.LumiList (filename = goldenJSONfile).getCMSSWString().split(',')
MuonPhysJSONList = LumiList.LumiList (filename = MuonPhysJSONfile).getCMSSWString().split(',')

process.inputs = cms.PSet (
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    )
process.inputs.lumisToProcess.extend(goldenJSONList)


process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring()
)

process.MessageLogger = cms.Service("MessageLogger")

process.WprimeAnalyzer = cms.PSet(
    ## common input for wrapped analyzers
    fileNames   = cms.vstring(),  ## keep empty!

    sample_cross_sections = cms.string("samples_cross_sections.txt"),
    outputFile  = cms.string("Wprime_analysis.root"),
    logFile = cms.string("Wprime_event_counts.txt"),
    candEvtFile = cms.string("Wprime_candEvts.txt"),

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
    runHadVWAnalysis = cms.bool(False),
    runTBAnalysis    = cms.bool(False),
    runWgammaAnalysis = cms.bool(False),

    ## input specific for this analyzer
    muonReconstructor = cms.int32(3), ## see TeVMuon_tracking.h

    muons = cms.InputTag('userPatMuons'),
    electrons = cms.InputTag('userPatElectrons'),
    jets = cms.InputTag('selectedPatJets'),
    met   = cms.InputTag('patMETsPF'),
    particleFlow = cms.InputTag('selectedPatPFParticles'),
    genParticles = cms.InputTag('prunedGenParticles'),
    hltEventTag = cms.InputTag('patTriggerEvent'),
    pileupTag  = cms.InputTag('addPileupInfo'),
    vertexTag  = cms.InputTag('offlinePrimaryVertices'),

    electronSelectors = electronSelectors,
    muonSelectors = muonSelectors,
    jetSelectors = jetSelectors,
    photonSelectors = photonSelectors,

    Cuts = cms.vstring(),
    debugme = cms.bool(False),
    preselect = cms.bool(False),
    triggersToUse = cms.vstring(),
    
    inputs = process.inputs,

    #PileUp Inputs
    MCPUDistFile = cms.string('UserCode/CMGWPrimeGroup/root_macros/MCSummer11S4PUDist.root'),
    MCPUDistHist = cms.string('probdistFlat10'),
    DataPUDistFile = cms.string('UserCode/CMGWPrimeGroup/root_macros/Cert_160404-178078_7TeV_PromptReco_Collisons11_JSON.pileupTruth_v2.root'),
    DataPUDistHist = cms.string('pileup'),

    vEventsToDebug = cms.VEventID(#'1:1:1', '2:2:22','3:33:3'
    ),
    
    
    )

