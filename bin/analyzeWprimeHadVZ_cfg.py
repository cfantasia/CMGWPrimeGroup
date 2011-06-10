import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("WPrimeAnalysis")
# get JSON file correctly parced
#JSONfile = 'UserCode/CMGWPrimeGroup/JSON/Cert_160404-165542_7TeV_PromptReco_Collisions11_JSON.txt'
#JSONfile = 'UserCode/CMGWPrimeGroup/JSON/json_160404-166011_DCSonly.txt'

#myList = LumiList.LumiList (filename = JSONfile).getCMSSWString().split(',')

#process.inputs = cms.PSet (
#    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
#    )
#process.inputs.lumisToProcess.extend(myList)

process.inputs = cms.PSet()

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
#    'file:patTuple.root'
  )
)

process.MessageLogger = cms.Service("MessageLogger")

from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import heepBarrelCuts, heepEndcapCuts


process.WprimeAnalyzer = cms.PSet(
    ## common input for wrapped analyzers
    fileNames   = cms.vstring(),  ## keep empty!
    # fileNames   = cms.vstring('file:patTuple.root'),  ## mandatory
    outputFile  = cms.string('test.root'),## mandatory
    maxEvents   = cms.int32(-1),                      ## optional
    reportAfter = cms.uint32(1),                     ## optional
    logFile     = cms.string("test.log"),
    candEvtFile = cms.string("testCandEvtFile.txt"),
    debugme     = cms.bool(False),
    preselect   = cms.bool(False),
    doRecoilCorrectionForW = cms.bool(False),
    sample_cross_sections = cms.string("samples_cross_sections_HadVZ.txt"),
    ## enable analysis in individual channels
    runMuMETAnalysis = cms.bool(False),
    runElMETAnalysis = cms.bool(False),
    runWZAnalysis    = cms.bool(False),
    runHadVZAnalysis = cms.bool(True),
    runTBAnalysis    = cms.bool(False),
    runWgammaAnalysis = cms.bool(False),
    ## input specific for this analyzer
    muons = cms.string('selectedPatMuons'),
    muonAlgo = cms.uint32(0),
    jets = cms.string('selectedPatJets'),
    genParticles = cms.InputTag('prunedGenParticles'),
    hltEventTag = cms.string('patTriggerEvent'),
    pileupTag  = cms.string('addPileupInfo'),
    triggersToUse = cms.vstring(),
    inputs = process.inputs,
    #
    MCPUDistFile = cms.string('UserCode/CMGWPrimeGroup/root_macros/MCPUDist.root'),
    MCPUDistHist = cms.string('pileup'),
    DataPUDistFile = cms.string('UserCode/CMGWPrimeGroup/root_macros/DataPUDist.root'),
    DataPUDistHist = cms.string('pileup'),
    #
    maxNumZs = cms.uint32(1),
    minNLeptons =cms.uint32(2),
    minNumJets = cms.uint32(0), # Larger than
    maxNumJets = cms.uint32(9999), # Smaller than 
    minLeadPt = cms.double(20.0),
    maxAngleBetweenJets = cms.double(9999.9),
    #
    minZpt = cms.double(0.0), # All units in GeV
    minZmass = cms.double(80.0),
    maxZmass = cms.double(100.0),
    minHadVPt = cms.double(0.0),
    minHadVmass = cms.double(0.0),
    maxHadVmass = cms.double(9999.9),
    # +++++++++++++++++++Muon General Cuts
    maxMuonEta = cms.double(2.5),
    minMuonLoosePt = cms.double(10.),
    minMuonTightPt = cms.double(20.),
    #VBTF Recommended Cuts
    maxMuonDxy = cms.double(0.2),
    maxMuonNormChi2 = cms.double(10.),
    minMuonNPixHit = cms.int32(0),
    minMuonNTrkHit = cms.int32(10),
    minMuonStations = cms.int32(1),
    minMuonHitsUsed = cms.int32(0),
    # +++++++++++++++++++Jet General Cuts
    minJetPt = cms.double(30.0),
    maxJetEta = cms.double(2.4),
    maxJetNHF = cms.double(0.99),
    maxJetNEF = cms.double(0.99),
    minJetnumConst = cms.uint32(1),
    minJetCHF = cms.double(0.0),
    minJetcMult = cms.uint32(0),
    maxJetCEF = cms.double(0.99)
) 
