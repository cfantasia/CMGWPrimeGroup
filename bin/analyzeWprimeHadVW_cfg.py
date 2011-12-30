from UserCode.CMGWPrimeGroup.commonWprime_cfg import *

process.WprimeAnalyzer.outputFile  = "HadVW.root"
process.WprimeAnalyzer.logFile     = "HadVW.log"
process.WprimeAnalyzer.candEvtFile = "HadVW.evt"
process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_HadVZ.txt"
#process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_HadVW.txt"

process.WprimeAnalyzer.maxEvents   = -1
process.WprimeAnalyzer.reportAfter = -10
process.WprimeAnalyzer.useJSON = False
process.WprimeAnalyzer.debug = False
process.WprimeAnalyzer.preselect = False

## enable analysis in individual channels
process.WprimeAnalyzer.runHadVWAnalysis = True
process.WprimeAnalyzer.triggersToUse = cms.vstring()

## input specific for this analyzer
process.WprimeAnalyzer.useAdjustedMET = False
process.WprimeAnalyzer.muonReconstructor = 3

process.WprimeAnalyzer.muons = 'selectedPatMuons'
process.WprimeAnalyzer.electrons = 'selectedPatElectrons'
process.WprimeAnalyzer.jets = 'selectedPatJetsAK5PF'

process.WprimeAnalyzer.minNLeptons =cms.untracked.uint32(1)
process.WprimeAnalyzer.minNJets = cms.untracked.uint32(1) 
#process.WprimeAnalyzer.maxAngleBetweenJets = cms.double(9999.9)
#
process.WprimeAnalyzer.minMET = cms.untracked.double(30.0) # All units in GeV
process.WprimeAnalyzer.minWtransMass = cms.untracked.double(30.0) # All units in GeV
process.WprimeAnalyzer.minWpt = cms.untracked.double(50.0) # All units in GeV

process.WprimeAnalyzer.minVmass = cms.untracked.double(70.0)
process.WprimeAnalyzer.maxVmass = cms.untracked.double(110.0)
process.WprimeAnalyzer.minVpt = cms.untracked.double(50.0)

process.WprimeAnalyzer.Cuts = cms.vstring(
    "NoCuts", 
    "MinNLeptons",
    "MinNJets",

    #    "HLT", 

    "ValidV", 
    "ValidW", 
    "ValidVWCand", 

    "MET",
    "WTransMass",
    "VMass",
    "Wpt",
    "Vpt",
    
    "AllCuts")

process.WprimeAnalyzer.LooseElectronType = cms.untracked.string("WZLoose")
process.WprimeAnalyzer.LooseMuonType = cms.untracked.string("HadVZLoose")
process.WprimeAnalyzer.LooseJetType = cms.untracked.string("Base")

