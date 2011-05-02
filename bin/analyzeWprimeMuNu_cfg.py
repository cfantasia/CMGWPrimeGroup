import FWCore.ParameterSet.Config as cms

process = cms.Process("WPrimeAnalysis")

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
##    'file:patTuple.root'
  )
)

process.MessageLogger = cms.Service("MessageLogger")

process.WprimeAnalyzer = cms.PSet(
    ## common input for wrapped analyzers
 fileNames   = cms.vstring(),  ## keep empty!
   # fileNames   = cms.vstring('file:patTuple.root'),  ## mandatory
    outputFile  = cms.string('Wprime_analysis.root'),## mandatory
    maxEvents   = cms.int32(-1),                      ## optional
    reportAfter = cms.uint32(5000),                     ## optional
    doRecoilCorrectionForW = cms.bool(False),
    sample_cross_sections = cms.string("samples_cross_sections_MuMET.txt"),
    ## enable analysis in individual channels
    runMuMETAnalysis = cms.bool(True),
    runElMETAnalysis = cms.bool(False),
    runWZAnalysis    = cms.bool(False),
    runTBAnalysis    = cms.bool(False),
    runWgammaAnalysis = cms.bool(False),
    ## input specific for this analyzer
    muons = cms.InputTag('selectedPatMuons'),
    met   = cms.InputTag('patMETsPF'),
    particleFlow = cms.InputTag('selectedPatPFParticles'),
    genParticles = cms.InputTag('prunedGenParticles'),
    #
    muonReconstructor = cms.int32(3), ## see mumet_histo_constants.h
    muonPtThreshold   = cms.double(10), ## in GeV
    oneMuPtTrackCut   = cms.double(25), ## in GeV
    chi2Cut           = cms.double(10),
    muonEtaCut        = cms.double(2.1),
    combRelCut        = cms.double(0.15),
    highestPtMuonOnly = cms.bool(False),
    dumpHighPtMuons   = cms.bool(True),
    dumpHighPtMuonThreshold = cms.double(200)
)

