import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("WPrimeAnalysis")
# get JSON file correctly parced
JSONfile = 'UserCode/CMGWPrimeGroup/JSON/Cert_160404-165542_7TeV_PromptReco_Collisions11_JSON.txt'
#JSONfile = 'UserCode/CMGWPrimeGroup/JSON/json_160404-166011_DCSonly.txt'

myList = LumiList.LumiList (filename = JSONfile).getCMSSWString().split(',')

process.inputs = cms.PSet (
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    )
process.inputs.lumisToProcess.extend(myList)

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
    reportAfter = cms.uint32(15000),                     ## optional
    doRecoilCorrectionForW = cms.bool(False),
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
    dumpHighEtElectrons   = cms.bool(True),
    dumpHighEtElectronThreshold = cms.double(200),
    barrelCuts = heepBarrelCuts,
    endcapCuts = heepEndcapCuts,
    inputs = process.inputs
 )

