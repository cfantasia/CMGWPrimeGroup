from UserCode.CMGWPrimeGroup.commonWprime_cfg import *

#process.WprimeAnalyzer.maxEvents   = cms.int32(500)
#process.WprimeAnalyzer.debugme = cms.bool(True)

process.WprimeAnalyzer.debugme = cms.bool(False)
process.WprimeAnalyzer.reportAfter = cms.uint32(25000)
process.WprimeAnalyzer.runWZAnalysis    = cms.bool(True)
process.WprimeAnalyzer.sample_cross_sections = cms.string("samples_cross_sections_WZ.txt")

process.WprimeAnalyzer.doRecoilCorrectionForW = cms.bool(False)
process.WprimeAnalyzer.useAdjustedMET = cms.bool(False)
process.WprimeAnalyzer.muonReconstructor = cms.uint32(7)

process.WprimeAnalyzer.useJSON = cms.bool(False)
process.WprimeAnalyzer.countGenEvts = cms.bool(False)
process.WprimeAnalyzer.eventCounters = cms.vstring(
    'nEventsTotal',
    'nEventsHLT',
    'nEventsFiltered',
    'nEventsPat')


## input specific for this analyzer
process.WprimeAnalyzer.muons = 'userPatMuons'
process.WprimeAnalyzer.electrons = 'userPatElectrons'
process.WprimeAnalyzer.met   = 'patMETsPF'
process.WprimeAnalyzer.particleFlow = 'selectedPatPFParticles'
process.WprimeAnalyzer.genParticles = 'prunedGenParticles'
process.WprimeAnalyzer.hltEventTag = 'patTriggerEvent'

process.WprimeAnalyzer.preselect = cms.bool(False)

process.WprimeAnalyzer.minDeltaR = cms.double(0.)
process.WprimeAnalyzer.maxZMassDiff = cms.double(999999.)

process.WprimeAnalyzer.effectiveElecArea = cms.vdouble(0.0997,0.1123)#Not using Recommended PI*0.3*0.3
process.WprimeAnalyzer.effectiveMuonArea = cms.vdouble(0.1057,0.0769)

###Analysis Cuts
EWKWZCuts = cms.vstring(
    "NoCuts", 
    "MinNLeptons",
    
    "ValidZ", 
#    "ZMass", 
#    "ZLepTrigMatch",
#    "ZLepPt",
    "HLT", 
    "NumZs", 
    
    "ValidW",    
    
    "MET",

    "AllCuts")
WprimeWZCuts = cms.vstring(
    "NoCuts", 
    "MinNLeptons",
    
    "ValidZ", 
#    "ZMass", 
#    "ZLepTrigMatch",
#    "ZLepPt",
    "HLT", 
    "NumZs", 
    
    "ValidW", 
    
    "MET",

    "ValidWZCand",   
    "Ht", 
#    "Zpt", 
#    "Wpt",

    "AllCuts")
WZFakeElecCuts = cms.vstring(
    "NoCuts", 
    "HLT",
    "FakeEvt",
    
    "ValidW", 
    "WFlavorMuon",
    
    "WTransMass",
    "MET",
    
    "FakeLepProbeLoose",
    "FakeLepProbeTight",
    )
WZFakeMuonCuts = cms.vstring(
    "NoCuts", 
    "HLT",
    "FakeEvt",
    
    "ValidW", 
    "WFlavorElec",
    
    "WTransMass",
    "MET",
    
    "FakeLepProbeLoose",
    "FakeLepProbeTight",
    )
WZEffCuts = cms.vstring(
    "NoCuts", 
    "MinNLeptons",
    "ValidZ", 
    "AllCuts")

process.WprimeAnalyzer.LooseElectronType = cms.string("WZLoose")
process.WprimeAnalyzer.TightElectronType = cms.string("WZTight")
process.WprimeAnalyzer.LooseMuonType = cms.string("WZLoose")
process.WprimeAnalyzer.TightMuonType = cms.string("WZTight")
process.WprimeAnalyzer.LooseJetType = cms.string("Base")

####Triggers
DoubleTriggers = cms.vstring(
    'HLT_DoubleMu5_v*', #For MC
    'HLT_DoubleMu7_v*',
    'HLT_Mu13_Mu8_v*', #1e33 unprescaled

    'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*',
    )
SingleElecTriggers = cms.vstring('HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*',
                                 'HLT_Ele17_CaloIdL_CaloIsoVL_v*','HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*',
                                 'HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v*',
                                 'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*',
                                 'HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL_v*','HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*',
                                 'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*',
                                 'HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v*','HLT_Ele45_CaloIdVT_TrkIdT_v*',
                                 'HLT_Ele17_SW_L1R_v*','HLT_Ele17_SW_Isol_L1R_v*',
                                 'HLT_Ele17_SW_TighterEleIdIsol_L1R_v*','HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v*',
                                 'HLT_Ele22_SW_L1R_v*','HLT_Ele22_SW_TighterCaloIdIsol_L1R_v*',
                                 'HLT_Ele22_SW_TighterEleId_L1R_v*','HLT_Ele32_SW_TighterEleId_L1R_v*'
                                 )

SingleMuonTriggers = cms.vstring('HLT_Mu15_v*','HLT_Mu17_v*','HLT_Mu19_v*','HLT_Mu20_v*','HLT_Mu21_v*','HLT_Mu24_v*','HLT_Mu25_v*','HLT_Mu30_v*')
process.WprimeAnalyzer.triggersToUse = DoubleTriggers

####################

# +++++++++++++++++++General Cut values
process.WprimeAnalyzer.maxNumZs = cms.uint32(1)
process.WprimeAnalyzer.minNLeptons = cms.untracked.uint32(3)
process.WprimeAnalyzer.minLeadPt = cms.double(35.)
process.WprimeAnalyzer.minMET = cms.untracked.double(30.)

# +++++++++++++++++++W Cuts
process.WprimeAnalyzer.minWtransMass = cms.untracked.double(0.)#Cory: Removed cut
process.WprimeAnalyzer.minWlepPt = cms.double(20.)

process.WprimeAnalyzer.cutWenuWPRelIsoMask = cms.int32(2)#Cory: Iso only
process.WprimeAnalyzer.cutElecWPTightType = cms.string("simpleEleId80relIso")

# +++++++++++++++++++Z Cuts
process.WprimeAnalyzer.minZeePt1 =  cms.double(20.)
process.WprimeAnalyzer.minZeePt2 =  cms.double(10.)
process.WprimeAnalyzer.minZmmPt1 =  cms.double(15.)
process.WprimeAnalyzer.minZmmPt2 =  cms.double(15.)

process.WprimeAnalyzer.minZmass =  cms.untracked.double(60.)
process.WprimeAnalyzer.maxZmass =  cms.untracked.double(120.)

# +++++++++++++++++++Analysis Cuts
process.WprimeAnalyzer.minHt = cms.double(250.)
process.WprimeAnalyzer.minZpt =  cms.untracked.double(0.)#90 for TC300) 110 for W'400
process.WprimeAnalyzer.minWpt = cms.untracked.double(0.)#90 for TC300) 110 for W'400
