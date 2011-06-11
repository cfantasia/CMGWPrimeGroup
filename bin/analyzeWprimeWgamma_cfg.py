import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("WPrimeAnalysis")
# get JSON file correctly parced
JSONfile = 'UserCode/CMGWPrimeGroup/JSON/Cert_160404-163757_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt'
#JSONfile = 'UserCode/CMGWPrimeGroup/JSON/json_160404-163869_DCSonly.txt'
myList = LumiList.LumiList (filename = JSONfile).getCMSSWString().split(',')

process.inputs = cms.PSet (
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    )
process.inputs.lumisToProcess.extend(myList)


process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    #'rfio:/castor/cern.ch/user/c/cleonido/wprime/V210/Data_Run2011A_7.426ipb.root ',
  )
)

process.MessageLogger = cms.Service("MessageLogger")

process.WprimeAnalyzer = cms.PSet(
    ## common input for wrapped analyzers
    fileNames   = cms.vstring(),  ## keep empty!
    ##fileNames   = cms.vstring('rfio:/castor/cern.ch/user/c/cleonido/wprime/V210/Data_Run2011A_7.426ipb.root'),  ## mandatory
    outputFile  = cms.string('Wgamma_analysis.root'),## mandatory
    maxEvents   = cms.int32(-1),                      ## optional
    reportAfter = cms.uint32(15000),                     ## optional
    doRecoilCorrectionForW = cms.bool(False),
    sample_cross_sections = cms.string("samples_cross_sections_MuMET.txt"),
    logFile = cms.string("Wprime_event_counts.txt"),
    ## enable analysis in individual channels
    runMuMETAnalysis = cms.bool(False),
    runElMETAnalysis = cms.bool(False),
    runWZAnalysis    = cms.bool(False),
    runHadronicVZAnalysis = cms.bool(False),
    runTBAnalysis    = cms.bool(False),
    runWgammaAnalysis = cms.bool(True),
    ## input specific for this analyzer
    muons = cms.InputTag('selectedPatMuons'),
    mets   = cms.InputTag('patMETsPF'),
    particleFlow = cms.InputTag('selectedPatPFParticles'),
    photons = cms.InputTag('selectedPatPhotons'),
    genParticles = cms.InputTag('prunedGenParticles'),
    #
    muonReconstructor = cms.int32(3), ## see TeVMuon_tracking.h
    # do not consider muons below this pt-threshold
    muonPtThreshold   = cms.double(10), ## in GeV
    oneMuPtTrackCut   = cms.double(25), ## in GeV
    chi2Cut           = cms.double(10),
    muonEtaCut        = cms.double(2.1),
    combRelCut        = cms.double(0.15),
    highestPtMuonOnly = cms.bool(True),
    highestPtPhotonOnly = cms.bool(False),
    dumpHighPtMuons   = cms.bool(True),
    dumpHighPtMuonThreshold = cms.double(200),
    dumpHighPtPhotons = cms.bool(True),
    dumpHighPtPhotonThreshold = cms.double(100),
    BarrelJurrasicECALIsoConst = cms.double(4.2),
    BarrelJurrasicECALIsoSlope = cms.double(0.006),
    BarrelTowerHCALIsoConst = cms.double(2.2),
    BarrelTowerHCALIsoSlope = cms.double(0.0025),
    BarrelMaxHadronicOverEm = cms.double(0.05),
    BarrelHollowConeTrkIsoConst = cms.double(3.5),
    BarrelHollowConeTrkIsoSlope = cms.double(0.001),
    BarrelMaxSigmaIetaIeta = cms.double(99999.9),
    EndcapJurrasicECALIsoConst = cms.double(4.2),
    EndcapJurrasicECALIsoSlope = cms.double(0.006),
    EndcapTowerHCALIsoConst = cms.double(2.2),
    EndcapTowerHCALIsoSlope = cms.double(0.0025),
    EndcapMaxHadronicOverEm = cms.double(0.05),
    EndcapHollowConeTrkIsoConst = cms.double(3.5),
    EndcapHollowConeTrkIsoSlope = cms.double(0.001),
    EndcapMaxSigmaIetaIeta = cms.double(99999.9),
    ApplyTrackVeto = cms.bool(True),
    minPt = cms.double(10),
    maxEta = cms.double(5),
    inputs = process.inputs,

    MCPUDistFile = cms.string('UserCode/CMGWPrimeGroup/root_macros/MCPUDist.root'),
    MCPUDistHist = cms.string('pileup'),
    DataPUDistFile = cms.string('UserCode/CMGWPrimeGroup/root_macros/DataPUDist.root'),
    DataPUDistHist = cms.string('pileup'),
)
