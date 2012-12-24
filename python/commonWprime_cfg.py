import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
from UserCode.CMGWPrimeGroup.selectors_cff import *
from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import heepBarrelCuts, heepEndcapCuts
from SHarper.HEEPAnalyzer.WP80SelectionCuts_cfi import wp80BarrelCuts, wp80EndcapCuts

process = cms.Process("WPrimeAnalysis")
# get JSON file correctly parsed 
goldenJSONfile = 'UserCode/CMGWPrimeGroup/JSON/json_190456_207898_analysis.txt'
dcsJSONfile = 'UserCode/CMGWPrimeGroup/JSON/json_190456-194912.txt'
MuonPhysJSONfile = 'UserCode/CMGWPrimeGroup/JSON/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt'
goldenJSONList = LumiList.LumiList (filename = goldenJSONfile).getCMSSWString().split(',')
dcsJSONList = LumiList.LumiList (filename = dcsJSONfile).getCMSSWString().split(',')
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
    reportAfter = cms.int32(5000),                  

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
    runTTbarAnalysis = cms.bool(False),

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
    debug = cms.bool(False),
    preselect = cms.bool(False),
    triggersToUse = cms.vstring(),
    
    inputs = process.inputs,

    #PileUp Inputs
    MCPUDistFile = cms.string('UserCode/CMGWPrimeGroup/pileup/MCPUDist.root'),
    MCPUDistHist = cms.string('Summer12Dist'),
#    DataPUDistFile = cms.string('UserCode/CMGWPrimeGroup/pileup/MCPUDist.root'),
#    DataPUDistHist = cms.string('Summer12Dist'),
    DataPUDistFile = cms.string('UserCode/CMGWPrimeGroup/pileup/DataPileupHistogram.root'),
    DataPUDistHist = cms.string('pileup'),
    puScale = cms.double(1.0),

    vEventsToDebug = cms.VEventID(#'1:1:1', '2:2:22','3:33:3'
    ),
    
    
    )

