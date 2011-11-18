from UserCode.CMGWPrimeGroup.commonWprime_cfg import *

process.WprimeAnalyzer.outputFile  = "TB.root"
process.WprimeAnalyzer.logFile     = "TB.log"
process.WprimeAnalyzer.candEvtFile = "TB.lst"
process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_TB.txt"

process.WprimeAnalyzer.maxEvents   = -1
process.WprimeAnalyzer.reportAfter = 10000
process.WprimeAnalyzer.useJSON = False
process.WprimeAnalyzer.debug = False
process.WprimeAnalyzer.preselect = False

## enable analysis in individual channels
process.WprimeAnalyzer.runTBAnalysis = True
process.WprimeAnalyzer.triggersToUse = cms.vstring()

## input specific for this analyzer
process.WprimeAnalyzer.electrons = 'selectedPatElectrons'
process.WprimeAnalyzer.muons = 'selectedPatMuons'
process.WprimeAnalyzer.muonReconstructor = 3
process.WprimeAnalyzer.jets = 'selectedPatJets'
#
process.WprimeAnalyzer.minNLeptons =cms.untracked.uint32(0)
process.WprimeAnalyzer.minNJets = cms.untracked.uint32(1) # Larger than
#

process.WprimeAnalyzer.Cuts = cms.vstring(
        "NoCuts",

        "MinNLeptons",
        "MinNJets",

        "ValidB",
        "ValidW",

        "ValidT",
#        "TMass",
#        "TPt",

        "ValidB",

        "ValidTB",
              
        #    "HLT",
        "AllCuts")

process.WprimeAnalyzer.LooseElectronType = cms.untracked.string("WZLoose")
process.WprimeAnalyzer.TightElectronType = cms.untracked.string("WZTight")

process.WprimeAnalyzer.LooseMuonType = cms.untracked.string("HadVZLoose")
process.WprimeAnalyzer.TightMuonType = cms.untracked.string("HadVZTight")

process.WprimeAnalyzer.LooseJetType = cms.untracked.string("Base")
