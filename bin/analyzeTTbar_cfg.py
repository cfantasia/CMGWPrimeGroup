from UserCode.CMGWPrimeGroup.commonWprime_cfg import *

process.WprimeAnalyzer.outputFile  = "TTbarAnalyzer.root"
process.WprimeAnalyzer.logFile     = "TTbarAnalyzer.log"
process.WprimeAnalyzer.candEvtFile = "TTbarAnalyzer.lst"
process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_TTbar.txt"

process.WprimeAnalyzer.maxEvents   = -1
process.WprimeAnalyzer.reportAfter = -10
process.WprimeAnalyzer.useJSON = False
process.WprimeAnalyzer.debug = False
process.WprimeAnalyzer.preselect = False

## Add proper Fall11 PU distribution
process.WprimeAnalyzer.MCPUDistHist = 'Fall11Dist'

## To run PU systematics jobs
##process.WprimeAnalyzer.puScale = 1.08

## enable analysis in individual channels
process.WprimeAnalyzer.runTTbarAnalysis = True
process.WprimeAnalyzer.triggersToUse = cms.vstring(
    #single mu
    'HLT_Mu15_v*',
    'HLT_Mu17_v*',
    'HLT_Mu20_v*',
    'HLT_Mu24_v*',
    'HLT_Mu30_v*',
    'HLT_Mu40_v*',
    'HLT_Mu40_eta2p1_v*',
    #double mu
#    'HLT_DoubleMu7_v*',#MC?
#    'HLT_Mu13_Mu8_v*', #1e33 unprescaled
#    'HLT_Mu17_Mu8_v*', #3e33 unprescaled
#    'HLT_Mu17_TkMu8_v*', #3e33 unprescaled

    #SingleE
    'HLT_Ele17_SW_L1R',  # from 2010 data
    'HLT_Ele17_SW_Isol_L1R_v*',
    'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*',
    'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*',
    'HLT_Ele42_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*',
    'HLT_Ele52_CaloIdVT_TrkIdT_v*',
    'HLT_Ele65_CaloIdVT_TrkIdT_v*',
    'HLT_Ele80_CaloIdVT_TrkIdT_v*',
    


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

process.WprimeAnalyzer.minVpt = cms.untracked.double(250.0)
#Nominal V Mass
process.WprimeAnalyzer.minVmass = cms.untracked.double(65.0)
process.WprimeAnalyzer.maxVmass = cms.untracked.double(120.0)
#Sideband V Mass
#process.WprimeAnalyzer.minVmass = cms.untracked.double(50.0)
#process.WprimeAnalyzer.maxVmass = cms.untracked.double(65.0)

process.WprimeAnalyzer.Cuts = cms.vstring(
        "NoCuts",

        "MinNLeptons",
        "MinNJets",
        "MinNBJets",
        "MinMET",

        "ValidW",
        "ValidT",
        "TMass",
        
        "ValidB2",
        "BPt2",

        "ValidT2",
        "w2pt",
        "T2Mass",
        "HadWMass",
        "w2pt_gt_200",
        "w2pt_gr_250",
        "AllCuts")

process.WprimeAnalyzer.LooseElectronType = cms.untracked.string("WZLoose")
process.WprimeAnalyzer.TightElectronType = cms.untracked.string("WZTight")

process.WprimeAnalyzer.LooseMuonType = cms.untracked.string("HadVZLoose")
process.WprimeAnalyzer.TightMuonType = cms.untracked.string("HadVZTight")

process.WprimeAnalyzer.LooseJetType = cms.untracked.string("Base")
