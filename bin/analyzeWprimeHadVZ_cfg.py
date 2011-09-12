from UserCode.CMGWPrimeGroup.commonWprime_cfg import *

process.WprimeAnalyzer.outputFile  = "HadVZAnalyzer.root"
process.WprimeAnalyzer.logFile     = "HadVZAnalyzer.log"
process.WprimeAnalyzer.candEvtFile = "HadVZAnalyzer.lst"
process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_HadVZ.txt"

process.WprimeAnalyzer.maxEvents   = -1
process.WprimeAnalyzer.reportAfter = 10000
process.WprimeAnalyzer.useJSON = False
process.WprimeAnalyzer.debugme = False
process.WprimeAnalyzer.preselect = False

## enable analysis in individual channels
process.WprimeAnalyzer.runHadVZAnalysis = True
process.WprimeAnalyzer.triggersToUse = ''

## input specific for this analyzer
process.WprimeAnalyzer.electrons = 'selectedPatElectrons'
process.WprimeAnalyzer.muons = 'selectedPatMuons'
process.WprimeAnalyzer.muonReconstructor = 3
process.WprimeAnalyzer.jets = 'selectedPatJets'
#
process.WprimeAnalyzer.minNLeptons =cms.untracked.uint32(2)
process.WprimeAnalyzer.minNJets = cms.untracked.uint32(1) # Larger than
process.WprimeAnalyzer.maxNJets = cms.uint32(99999) # Smaller than
process.WprimeAnalyzer.maxAngleBetweenJets = cms.double(9999.9)
#
process.WprimeAnalyzer.minZpt = cms.untracked.double(50.0) # All units in GeV
process.WprimeAnalyzer.minZmass = cms.untracked.double(70.0)
process.WprimeAnalyzer.maxZmass = cms.untracked.double(110.0)
process.WprimeAnalyzer.minVPt = cms.untracked.double(50.0)
process.WprimeAnalyzer.minVmass = cms.untracked.double(60.0)
process.WprimeAnalyzer.maxVmass = cms.untracked.double(110.0)

process.WprimeAnalyzer.Cuts = cms.vstring(
        "NoCuts",

        "MinNLeptons",
        "MinNJets",

        "ValidZ",
        "ZMass",
        "Zpt",

        "ValidV",
        "VMass",
        "VPt",

        "ValidVZ",
              
        #    "HLT",
        "AllCuts")

process.WprimeAnalyzer.LooseElectronType = cms.untracked.string("WZLoose")
process.WprimeAnalyzer.TightElectronType = cms.untracked.string("WZTight")

process.WprimeAnalyzer.LooseMuonType = cms.untracked.string("HadVZLoose")
process.WprimeAnalyzer.TightMuonType = cms.untracked.string("HadVZTight")

process.WprimeAnalyzer.LooseJetType = cms.untracked.string("Base")
