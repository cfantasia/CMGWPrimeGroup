from UserCode.CMGWPrimeGroup.commonWprime_cfg import *

process.WprimeAnalyzer.outputFile  = "TB.root"
process.WprimeAnalyzer.logFile     = "TB.log"
process.WprimeAnalyzer.candEvtFile = "TB.lst"
process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_HadVZ.txt"
#process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_TB.txt"

process.WprimeAnalyzer.maxEvents   = 1000
process.WprimeAnalyzer.reportAfter = 10000
process.WprimeAnalyzer.useJSON = False
process.WprimeAnalyzer.debug = True
process.WprimeAnalyzer.preselect = False

## enable analysis in individual channels
process.WprimeAnalyzer.runTBAnalysis = True
process.WprimeAnalyzer.triggersToUse = cms.vstring()

## input specific for this analyzer
process.WprimeAnalyzer.electrons = 'selectedPatElectrons'
process.WprimeAnalyzer.muons = 'selectedPatMuons'
process.WprimeAnalyzer.muonReconstructor = 3
process.WprimeAnalyzer.jets = 'selectedPatJetsAK5PF'
#
process.WprimeAnalyzer.minNLeptons =cms.untracked.uint32(1)
process.WprimeAnalyzer.minNJets = cms.untracked.uint32(1)
process.WprimeAnalyzer.minNBJets = cms.untracked.uint32(1)
process.WprimeAnalyzer.minMET = cms.untracked.double(30)

process.WprimeAnalyzer.minWtransMass = cms.untracked.double(30)
process.WprimeAnalyzer.minWpt = cms.untracked.double(30)

process.WprimeAnalyzer.minTpt = cms.untracked.double(50)
process.WprimeAnalyzer.minTMass = cms.untracked.double(150)
process.WprimeAnalyzer.maxTMass = cms.untracked.double(200)

process.WprimeAnalyzer.minBDisc = cms.untracked.double(1.7)##???
process.WprimeAnalyzer.BDisc = cms.untracked.string("trackCountingHighEffBJetTags")
process.WprimeAnalyzer.maxBMass = cms.untracked.double(10)
process.WprimeAnalyzer.minBpt = cms.untracked.double(50)

#

process.WprimeAnalyzer.Cuts = cms.vstring(
        "NoCuts",

        "MinNLeptons",
        "MinNJets",
        "MinNBJets",
        "MET",

#        "ValidB",
        "ValidW",

        "ValidT",
        "TMass",
        "TPt",

        "ValidB",
        "BPt2",

        "ValidTB",
              
        #    "HLT",
        "AllCuts")

process.WprimeAnalyzer.LooseElectronType = cms.untracked.string("WZLoose")
process.WprimeAnalyzer.TightElectronType = cms.untracked.string("WZTight")

process.WprimeAnalyzer.LooseMuonType = cms.untracked.string("HadVZLoose")
process.WprimeAnalyzer.TightMuonType = cms.untracked.string("HadVZTight")

process.WprimeAnalyzer.LooseJetType = cms.untracked.string("Base")
