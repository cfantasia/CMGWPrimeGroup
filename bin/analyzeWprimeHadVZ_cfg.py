from UserCode.CMGWPrimeGroup.commonWprime_cfg import *

process.WprimeAnalyzer.outputFile  = "HadVZAnalyzer_Fall11.root"
process.WprimeAnalyzer.logFile     = "HadVZAnalyzer_Fall11.log"
process.WprimeAnalyzer.candEvtFile = "HadVZAnalyzer_Fall11.lst"
process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_HadVZ.txt"

process.WprimeAnalyzer.maxEvents   = -1
process.WprimeAnalyzer.reportAfter = -10
process.WprimeAnalyzer.useJSON = False
process.WprimeAnalyzer.debug = False
process.WprimeAnalyzer.preselect = False

## Add proper Fall11 PU distribution
process.WprimeAnalyzer.MCPUDistHist = 'Fall11Dist'

## enable analysis in individual channels
process.WprimeAnalyzer.runHadVZAnalysis = True
process.WprimeAnalyzer.triggersToUse = cms.vstring(
    #single mu
    'HLT_Mu15_v*',
    'HLT_Mu17_v*',
    'HLT_Mu20_v*',
    'HLT_Mu24_v*',
    'HLT_Mu30_v*',
    'HLT_Mu40_v*',
    #double mu
    'HLT_DoubleMu7_v*',#MC?
    'HLT_Mu13_Mu8_v*', #1e33 unprescaled
    'HLT_Mu17_Mu8_v*', #3e33 unprescaled
    'HLT_Mu17_TkMu8_v*', #3e33 unprescaled

    #Double E (1 for mc, 1 for data)
    'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*',
    'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*',#MC
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
#Nominal V Mass
process.WprimeAnalyzer.minVmass = cms.untracked.double(70.0)
process.WprimeAnalyzer.maxVmass = cms.untracked.double(120.0)
#Sideband V Mass
#process.WprimeAnalyzer.minVmass = cms.untracked.double(30.0)
#process.WprimeAnalyzer.maxVmass = cms.untracked.double(70.0)
#

process.WprimeAnalyzer.Cuts = cms.vstring(
        "NoCuts",

        "MinNLeptons",
        "MinNJets",

        "ValidZ",
        "ZMass",
        "Zpt",
        
        "ValidV",
        "VMass",

        "ValidVZ",
              
        "HLT",


        "VPt",
        
        "AllCuts")

process.WprimeAnalyzer.LooseElectronType = cms.untracked.string("WZLoose")
process.WprimeAnalyzer.TightElectronType = cms.untracked.string("WZTight")

process.WprimeAnalyzer.LooseMuonType = cms.untracked.string("HadVZLoose")
process.WprimeAnalyzer.TightMuonType = cms.untracked.string("HadVZTight")

process.WprimeAnalyzer.LooseJetType = cms.untracked.string("Base")
