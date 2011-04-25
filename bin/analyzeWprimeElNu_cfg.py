import FWCore.ParameterSet.Config as cms

process = cms.Process("WPrimeAnalysis")

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
##    'file:patTuple.root'
  )
)

process.MessageLogger = cms.Service("MessageLogger")

from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import heepBarrelCuts, heepEndcapCuts


process.WprimeAnalyzer = cms.PSet(
    ## common input for wrapped analyzers
 fileNames   = cms.vstring(),  ## keep empty!
   # fileNames   = cms.vstring('file:patTuple.root'),  ## mandatory
    outputFile  = cms.string('Wprime_analysis.root'),## mandatory
    maxEvents   = cms.int32(-1),                      ## optional
    reportAfter = cms.uint32(1000),                     ## optional
    doRecoilCorrectionForW = cms.bool(True),
    sample_cross_sections = cms.string("samples_cross_sections_ElMET.txt"),
    ## enable analysis in individual channels
    runMuMETAnalysis = cms.bool(False),
    runElMETAnalysis = cms.bool(True),
    runWZAnalysis    = cms.bool(False),
    runTBAnalysis    = cms.bool(False),
    runWgammaAnalysis = cms.bool(False),
    ## input specific for this analyzer
    electrons = cms.InputTag('selectedPatElectrons'),
    met   = cms.InputTag('patMETsPF'),
    genParticles = cms.InputTag('prunedGenParticles'),
    #
    oneEleEtCut   = cms.double(25), ## in GeV
    highestEtElectronOnly = cms.bool(False),
    dumpHighEtElectrons   = cms.bool(False),
    dumpHighEtElectronThreshold = cms.double(200),
    barrelCuts = heepBarrelCuts,
    endcapCuts = heepEndcapCuts
)

