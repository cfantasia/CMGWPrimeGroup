from UserCode.CMGWPrimeGroup.commonWprime_cfg import *

process.WprimeAnalyzer.outputFile  = "HadVZAnalyzer.root"
process.WprimeAnalyzer.logFile     = "HadVZAnalyzer.log"
process.WprimeAnalyzer.candEvtFile = "HadVZAnalyzer.lst"
process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_HadVZ.txt"

process.WprimeAnalyzer.maxEvents   = -1
process.WprimeAnalyzer.reportAfter = 100000
process.WprimeAnalyzer.useJSON = False
process.WprimeAnalyzer.debug = False
process.WprimeAnalyzer.preselect = False

## enable analysis in individual channels
process.WprimeAnalyzer.runHadVZAnalysis = True
process.WprimeAnalyzer.triggersToUse = cms.vstring(
    'HLT_Mu11',          # from 2010 data
    'HLT_Mu15_v*',
    'HLT_Mu17_v*',
    'HLT_Mu20_v*',
    'HLT_Mu24_v*',
    'HLT_Mu30_v*',
    'HLT_IsoMu15_v*',
    'HLT_IsoMu17_v*',
    'HLT_IsoMu24_v*',
    'HLT_IsoMu30_v*',
    'HLT_IsoMu40_v*',
    'HLT_DoubleMu5_v*',
    'HLT_DoubleMu7_v*',
    'HLT_TripleMu5_v*',
    'HLT_Mu13_Mu8_v*', #1e33 unprescaled
    'HLT_Mu17_Mu8_v*' #3e33 unprescaled
    'HLT_Ele17_SW_L1R',  # from 2010 data
    'HLT_Ele17_SW_Isol_L1R_v*',
    'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*',
    'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*',
    'HLT_Ele42_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*' #2e33(?) unprescaled
    'HLT_DoubleEle17_SW_L1R_v*',
    'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*',
    'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',
    )

## input specific for this analyzer
process.WprimeAnalyzer.electrons = 'selectedPatElectrons'
process.WprimeAnalyzer.muons = 'selectedPatMuons'
process.WprimeAnalyzer.muonReconstructor = 3
process.WprimeAnalyzer.jets = 'selectedPatJetsAK7PF'
#
process.WprimeAnalyzer.minNLeptons =cms.untracked.uint32(2)
process.WprimeAnalyzer.minNJets = cms.untracked.uint32(1) # Larger than
process.WprimeAnalyzer.maxNJets = cms.uint32(3) # Smaller than
process.WprimeAnalyzer.maxAngleBetweenJets = cms.double(9999.9)
#
process.WprimeAnalyzer.minZpt = cms.untracked.double(150.0) # All units in GeV
process.WprimeAnalyzer.minZmass = cms.untracked.double(70.0)
process.WprimeAnalyzer.maxZmass = cms.untracked.double(110.0)
process.WprimeAnalyzer.minVpt = cms.untracked.double(290.0)
process.WprimeAnalyzer.minVmass = cms.untracked.double(70.0)
process.WprimeAnalyzer.maxVmass = cms.untracked.double(120.0)
#

process.WprimeAnalyzer.Cuts = cms.vstring(
        "NoCuts",

        "MinNLeptons",
        "MinNJets",

        "ValidZ",
        "ZMass",

        "ValidV",
        "VMass",

        "ValidVZ",
              
        "HLT",

        "Zpt",
        "VPt",
        
        "AllCuts")

process.WprimeAnalyzer.LooseElectronType = cms.untracked.string("WZLoose")
process.WprimeAnalyzer.TightElectronType = cms.untracked.string("WZTight")

process.WprimeAnalyzer.LooseMuonType = cms.untracked.string("HadVZLoose")
process.WprimeAnalyzer.TightMuonType = cms.untracked.string("HadVZTight")

process.WprimeAnalyzer.LooseJetType = cms.untracked.string("Base")
