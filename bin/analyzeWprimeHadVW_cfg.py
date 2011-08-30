from UserCode.CMGWPrimeGroup.commonWprime_cfg import *

process.WprimeAnalyzer.outputFile  = 'HadVW_analysis.root'
process.WprimeAnalyzer.logFile     = "HadVW_event_counts.txt"
process.WprimeAnalyzer.candEvtFile = cms.string("HadVW_CandEvts.txt")
process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_HadVW.txt"
process.WprimeAnalyzer.maxEvents   = 100
process.WprimeAnalyzer.reportAfter = 1000
process.WprimeAnalyzer.useJSON = False
process.WprimeAnalyzer.debugme = cms.bool(True)
process.WprimeAnalyzer.preselect = cms.bool(False)

## enable analysis in individual channels
process.WprimeAnalyzer.runHadVWAnalysis = True
process.WprimeAnalyzer.triggersToUse = cms.vstring()

## input specific for this analyzer
process.WprimeAnalyzer.useAdjustedMET = False
process.WprimeAnalyzer.muonReconstructor = 7

process.WprimeAnalyzer.muons = 'selectedPatMuons'
process.WprimeAnalyzer.electrons = 'selectedPatElectrons'
process.WprimeAnalyzer.jets = 'selectedPatJets'

process.WprimeAnalyzer.minNLeptons =cms.untracked.uint32(1)
process.WprimeAnalyzer.minNJets = cms.untracked.uint32(1) 
process.WprimeAnalyzer.maxAngleBetweenJets = cms.double(9999.9)
#
process.WprimeAnalyzer.minHt = cms.double(0.0) # All units in GeV

process.WprimeAnalyzer.minMET = cms.untracked.double(30.0) # All units in GeV
process.WprimeAnalyzer.minWtransMass = cms.untracked.double(0.0) # All units in GeV
process.WprimeAnalyzer.minWpt = cms.untracked.double(0.0) # All units in GeV

process.WprimeAnalyzer.minVmass = cms.untracked.double(70.0)
process.WprimeAnalyzer.maxVmass = cms.untracked.double(110.0)
process.WprimeAnalyzer.minVpt = cms.untracked.double(0.0)

process.WprimeAnalyzer.Cuts = cms.vstring(
    "NoCuts", 
    "MinNLeptons",
    "MinNJets",

    #    "HLT", 

    "ValidV", 
    "ValidW", 
    "MET",

    
    "AllCuts")

process.WprimeAnalyzer.LooseElectronType = cms.string("WZLoose")
process.WprimeAnalyzer.TightElectronType = cms.string("WZTight")
process.WprimeAnalyzer.LooseMuonType = cms.string("WZLoose")
process.WprimeAnalyzer.TightMuonType = cms.string("WZTight")
process.WprimeAnalyzer.LooseJetType = cms.string("Base")

process.WprimeAnalyzer.effectiveElecArea = cms.vdouble(0.0997,0.1123)#Not using Recommended PI*0.3*0.3
process.WprimeAnalyzer.effectiveMuonArea = cms.vdouble(0.1057,0.0769)
