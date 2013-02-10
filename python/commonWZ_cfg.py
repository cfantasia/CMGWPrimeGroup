from UserCode.CMGWPrimeGroup.commonWprime_cfg import *

#process.WprimeAnalyzer.maxEvents   = 100
#process.WprimeAnalyzer.debug = True

process.WprimeAnalyzer.reportAfter = -5
process.WprimeAnalyzer.runWZAnalysis    = True
process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_WZ.txt"

process.WprimeAnalyzer.doRecoilCorrectionForW = False
process.WprimeAnalyzer.useAdjustedMET = False
process.WprimeAnalyzer.muonReconstructor = 7 #kPat
#process.WprimeAnalyzer.muonReconstructor = 3 #kCocktail

process.WprimeAnalyzer.useJSON = True
process.WprimeAnalyzer.countGenEvts = False
process.WprimeAnalyzer.eventCounters = cms.vstring(
    'nEventsTotal',
    'nEventsHLT',
    'nEventsFiltered',
    'nEventsPat')
#process.WprimeAnalyzer.puScale = 1.08

## input specific for this analyzer
process.WprimeAnalyzer.muons = 'selectedPatMuons'
process.WprimeAnalyzer.electrons = 'selectedPatElectrons'
process.WprimeAnalyzer.met   = "pfType1CorrectedMet"
process.WprimeAnalyzer.particleFlow = 'selectedPatPFParticles'
process.WprimeAnalyzer.genParticles = 'prunedGenParticles'
process.WprimeAnalyzer.hltEventTag = 'patTriggerEvent'
process.WprimeAnalyzer.rhoFastJet = cms.InputTag('kt6PFJets:rho')

process.WprimeAnalyzer.preselect = True

process.WprimeAnalyzer.minDeltaR = cms.double(0.)
process.WprimeAnalyzer.maxZMassDiff = cms.double(999999.)

process.WprimeAnalyzer.effectiveElecArea = cms.vdouble(0.0997,0.1123)#Not using Recommended PI*0.3*0.3
process.WprimeAnalyzer.effectiveMuonArea = cms.vdouble(0.1057,0.0769)

###Analysis Cuts
EWKWZCuts = cms.vstring(
    "NoCuts", 
    "MinNLeptons",
    
    "ValidZ", 
    "HLT", 
    "NumZs", 
    
    "ValidW",    
    
    "MET",
    )
WprimeWZCuts = cms.vstring(
    "NoCuts", 
    "MinNLeptons",
    
    "ValidZ", 
    "HLT", 
    "NumZs", 
    
    "ValidW", 
    
    "MET",

#    "ValidWZCand",   
    "Lt", 
    )
WZFakeElecCuts = cms.vstring(
    "NoCuts", 
    "MaxNVLLeptons",
    "MinNLeptons",
    "MinNTightLeptons",
    "FakeEvt",
    "HLT",
    
    "ValidW", 
    "WFlavorMuon",
    
    "WTransMass",
    "MaxWTransMass",
#    "MaxNJets",
    "MET",
    "MaxDeltaWTransMass",
    
    "FakeLepProbe",
    )
WZFakeMuonCuts = cms.vstring(
    "NoCuts", 
    "MaxNVLLeptons",
    "MinNLeptons",
    "MinNTightLeptons",
    "FakeEvt",
    "HLT",
    
    "ValidW", 
    "WFlavorElec",
    
    "WTransMass",
    "MaxWTransMass",
#    "MaxNJets",
    "MET",
    "MaxDeltaWTransMass",
    
    "FakeLepProbe",
    )
WZEffCuts = cms.vstring(
    "NoCuts", 
    "MinNLeptons",
    "ValidZ", 
    )

process.WprimeAnalyzer.LooseZElectronType = cms.untracked.string("EWKWZLoose")
process.WprimeAnalyzer.TightZElectronType = cms.untracked.string("EWKWZTight")
process.WprimeAnalyzer.LooseWElectronType = cms.untracked.string("EWKWZLoose")
process.WprimeAnalyzer.TightWElectronType = cms.untracked.string("EWKWZTight")
process.WprimeAnalyzer.LooseZMuonType = cms.untracked.string("EWKWZLoose")
process.WprimeAnalyzer.TightZMuonType = cms.untracked.string("EWKWZTight")
process.WprimeAnalyzer.LooseWMuonType = cms.untracked.string("EWKWZLoose")
process.WprimeAnalyzer.TightWMuonType = cms.untracked.string("EWKWZTight")
process.WprimeAnalyzer.LooseJetType = cms.untracked.string("Pat")

####Triggers
DoubleTriggers = cms.vstring(
    'HLT_DoubleMu7_v*',
    'HLT_Mu13_Mu8_v*', #1e33 unprescaled
    'HLT_Mu17_Mu8_v*', #3e33 unprescaled
    'HLT_Mu17_TkMu8_v*', #2011B
    'HLT_Mu22_TkMu8_v*', #2012 single l1 seeded
    'HLT_Mu22_TkMu22_v*', #2012

    'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*', #2011A
    'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*',#MC
    'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*'#Data
    'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*',
    'HLT_DoubleEle33_CaloIdT_v*',
    
    )
SingleElecTriggers = cms.vstring( 'HLT_Ele17_SW_L1R',  # from 2010 data
                                  'HLT_Ele17_SW_Isol_L1R_v*',
                                  'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*',
                                  'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*',
                                  'HLT_Ele42_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*', #1e33 unprescaled (removed after end of June TS)
                                  'HLT_Ele52_CaloIdVT_TrkIdT_v*',                  #1e33 unprescaled
                                  'HLT_Ele65_CaloIdVT_TrkIdT_v*',                  #3e33 unprescaled
                                  'HLT_Ele80_CaloIdVT_TrkIdT_v*',
                                  'HLT_Ele27_WP80_v*',
                                  'HLT_Ele27_WP80_PFMET_MT50_v*',
                                  'HLT_Ele80_CaloIdVT_TrkIdT_v*',
                                  )

SingleMuonTriggers = cms.vstring(    'HLT_Mu15_v*',
                                     'HLT_Mu17_v*',
                                     'HLT_Mu19_v*',
                                     'HLT_Mu20_v*',
                                     'HLT_Mu21_v*',
                                     'HLT_Mu24_v*',
                                     'HLT_Mu25_v*',
                                     'HLT_Mu30_v*',
                                     'HLT_Mu40_v*',         # 2.5e33 unprescaled
                                     'HLT_Mu40_eta2p1_v*',
                                     'HLT_IsoMu15_v*',
                                     'HLT_IsoMu17_v*',
                                     'HLT_IsoMu24_v*',
                                     'HLT_IsoMu30_v*',
                                     'HLT_IsoMu40_v*',
                                     'HLT_IsoMu24_eta2p1_v*',
                                     )
process.WprimeAnalyzer.triggersToUse = DoubleTriggers

####################
#process.WprimeAnalyzer.removeTauEvents = cms.untracked.bool(True)
process.WprimeAnalyzer.adjustMETPhi = cms.untracked.bool(True)
process.WprimeAnalyzer.elScaleFactor = cms.double(0.)
process.WprimeAnalyzer.muScaleFactor = cms.double(0.)

# +++++++++++++++++++General Cut values
process.WprimeAnalyzer.maxNumZs = cms.uint32(1)
process.WprimeAnalyzer.minNLeptons = cms.untracked.uint32(3)
process.WprimeAnalyzer.minLeadPt = cms.double(35.)
process.WprimeAnalyzer.minMET = cms.untracked.double(30.)

# +++++++++++++++++++W Cuts
process.WprimeAnalyzer.minWleptPt = cms.untracked.double(0.)
process.WprimeAnalyzer.minWtransMass = cms.untracked.double(0.)#Cory: Removed cut
process.WprimeAnalyzer.minWlepPt = cms.double(20.)

process.WprimeAnalyzer.cutWenuWPRelIsoMask = cms.int32(2)#Cory: Iso only
process.WprimeAnalyzer.cutElecWPTightType = cms.string("simpleEleId80relIso")

# +++++++++++++++++++Z Cuts
process.WprimeAnalyzer.minZeePt1 =  cms.double(35.)
process.WprimeAnalyzer.minZeePt2 =  cms.double(35.)
process.WprimeAnalyzer.minZmmPt1 =  cms.double(25.)
process.WprimeAnalyzer.minZmmPt2 =  cms.double(10.)

process.WprimeAnalyzer.minZmass =  cms.untracked.double( 71.188)
process.WprimeAnalyzer.maxZmass =  cms.untracked.double(111.188)

