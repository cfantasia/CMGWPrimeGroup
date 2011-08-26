from UserCode.CMGWPrimeGroup.commonWprime_cfg import *

process.WprimeAnalyzer.outputFile  = cms.string('HadVW_analysis.root')## mandatory
process.WprimeAnalyzer.logFile     = cms.string("HadVW_event_counts.txt")
process.WprimeAnalyzer.candEvtFile = cms.string("HadVW_CandEvts.txt")
process.WprimeAnalyzer.sample_cross_sections = cms.string("samples_cross_sections_HadVW.txt")
process.WprimeAnalyzer.maxEvents   = cms.int32(100)                      ## optional
process.WprimeAnalyzer.reportAfter = cms.uint32(1000)                     ## optional
process.WprimeAnalyzer.useJSON = cms.bool(False)
process.WprimeAnalyzer.debugme = cms.bool(True)
process.WprimeAnalyzer.preselect = cms.bool(False)

## enable analysis in individual channels
process.WprimeAnalyzer.runHadVWAnalysis = cms.bool(True)
process.WprimeAnalyzer.triggersToUse = cms.vstring()

## input specific for this analyzer
process.WprimeAnalyzer.useAdjustedMET = cms.bool(False)
process.WprimeAnalyzer.muonReconstructor = cms.uint32(7)

process.WprimeAnalyzer.muons = 'userPatMuons'
process.WprimeAnalyzer.electrons = 'userPatElectrons'

process.WprimeAnalyzer.minNLeptons =cms.uint32(1)
process.WprimeAnalyzer.maxNLeptons =cms.uint32(9999999)
process.WprimeAnalyzer.minNJets = cms.uint32(1) 
process.WprimeAnalyzer.maxAngleBetweenJets = cms.double(9999.9)
#
process.WprimeAnalyzer.minHt = cms.double(0.0) # All units in GeV

process.WprimeAnalyzer.minMET = cms.double(30.0) # All units in GeV
process.WprimeAnalyzer.minWtransMass = cms.double(0.0) # All units in GeV
process.WprimeAnalyzer.minWpt = cms.double(0.0) # All units in GeV

process.WprimeAnalyzer.minVmass = cms.double(70.0)
process.WprimeAnalyzer.maxVmass = cms.double(110.0)
process.WprimeAnalyzer.minVpt = cms.double(0.0)

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
