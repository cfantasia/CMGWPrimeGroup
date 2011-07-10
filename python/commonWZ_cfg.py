from UserCode.CMGWPrimeGroup.commonWprime_cfg import *

process.WprimeAnalyzer.debugme = cms.bool(False)
process.WprimeAnalyzer.runWZAnalysis    = cms.bool(True)
process.WprimeAnalyzer.sample_cross_sections = cms.string("samples_cross_sections_WZ.txt")

process.WprimeAnalyzer.doRecoilCorrectionForW = cms.bool(False)
process.WprimeAnalyzer.useAdjustedMET = cms.bool(False)
process.WprimeAnalyzer.muonAlgo = cms.uint32(0)

process.WprimeAnalyzer.useJSON = cms.bool(False)
process.WprimeAnalyzer.countGenEvts = cms.bool(True)
process.WprimeAnalyzer.eventCounters = cms.vstring(
    'nEventsTotal',
    'nEventsHLT',
    'nEventsFiltered',
    'nEventsPat')


## input specific for this analyzer
process.WprimeAnalyzer.muons = cms.string('userPatMuons')
process.WprimeAnalyzer.electrons = cms.string('userPatElectrons')
process.WprimeAnalyzer.met   = cms.string('patMETsPF')
process.WprimeAnalyzer.particleFlow = cms.InputTag('selectedPatPFParticles')
process.WprimeAnalyzer.genParticles = cms.InputTag('prunedGenParticles')
process.WprimeAnalyzer.hltEventTag = cms.string('patTriggerEvent')

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
    "ZMass", 
    "NumZs", 
    "ZLepPt",
    "ZLepTrigMatch",
    "HLT", 
    
    "ValidW",    
    "EvtSetup",
    
    "MET",

    "AllCuts")
WprimeWZCuts = cms.vstring(
    "NoCuts", 
    "MinNLeptons",
    
    "ValidZ", 
    "ZMass", 
    "NumZs", 
    "ZLepPt",
    "ZLepTrigMatch",
    "HLT", 
    
    "ValidW", 
    
    "MET",

    "ValidWZCand",
    "EvtSetup",
    
    "Ht", 
    "Zpt", 
    "Wpt",

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
    ),
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
    ),
WZEffCuts = cms.vstring(
    "NoCuts", 
    "MinNLeptons",

    "ValidZ", 
    "ZLepPt",
    
    "AllCuts"),    

process.WprimeAnalyzer.LooseElectronType = cms.string("WZLoose")
process.WprimeAnalyzer.TightElectronType = cms.string("WZTight")
process.WprimeAnalyzer.LooseMuonType = cms.string("WZLoose")
process.WprimeAnalyzer.TightMuonType = cms.string("WZTight")

####Triggers
DoubleTriggers = cms.vstring(
    'HLT_DoubleMu5_v*', #For MC
    'HLT_DoubleMu7_v*',
    'HLT_Mu13_Mu8_v*', #1e33 unprescaled

    'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*',
    )

SingleElecTriggers = cms.vstring('HLT_Ele10_SW_L1R_v2','HLT_Ele12_SW_TighterEleId_L1R_v2',
                                 'HLT_Ele17_SW_L1R_v2','HLT_Ele17_SW_Isol_L1R_v2',
                                 'HLT_Ele17_SW_TighterEleIdIsol_L1R_v3','HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v2',
                                 'HLT_Ele22_SW_L1R_v2','HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2',
                                 'HLT_Ele22_SW_TighterEleId_L1R_v3','HLT_Ele32_SW_TighterEleId_L1R_v2')

SingleMuonTriggers = cms.vstring('HLT_Mu3_v2','HLT_Mu5','HLT_Mu7','HLT_Mu9','HLT_Mu11',
                                 'HLT_Mu13_v1','HLT_Mu15_v1','HLT_Mu17_v1','HLT_Mu19_v1',
                                 'HLT_Mu21_v1','HLT_Mu25_v1')

process.WprimeAnalyzer.triggersToUse = DoubleTriggers

####################

# +++++++++++++++++++General Cut values
process.WprimeAnalyzer.maxNumZs = cms.uint32(1)
process.WprimeAnalyzer.minNLeptons = cms.uint32(3)
process.WprimeAnalyzer.maxNLeptons = cms.uint32(3)
process.WprimeAnalyzer.minLeadPt = cms.double(35.)
process.WprimeAnalyzer.minMET = cms.double(30.)

# +++++++++++++++++++Ht Cuts
process.WprimeAnalyzer.minHt = cms.double(190.)#150 for TC300) 190 for W'400

# +++++++++++++++++++W Cuts
process.WprimeAnalyzer.minWtransMass = cms.double(0.)#Cory: Removed cut
process.WprimeAnalyzer.minWpt = cms.double(110.)#90 for TC300) 110 for W'400

process.WprimeAnalyzer.minWlepPt = cms.double(20.)

process.WprimeAnalyzer.cutWenuWPRelIsoMask = cms.int32(2)#Cory: Iso only
process.WprimeAnalyzer.cutElecWPTightType = cms.string("simpleEleId80relIso")

# +++++++++++++++++++Z Cuts
process.WprimeAnalyzer.minZeePt1 =  cms.double(20.)
process.WprimeAnalyzer.minZeePt2 =  cms.double(10.)
process.WprimeAnalyzer.minZmmPt1 =  cms.double(15.)
process.WprimeAnalyzer.minZmmPt2 =  cms.double(15.)
process.WprimeAnalyzer.minZpt =  cms.double(110.)#90 for TC300) 110 for W'400
process.WprimeAnalyzer.minZmass =  cms.double(60.)
process.WprimeAnalyzer.maxZmass =  cms.double(120.)
