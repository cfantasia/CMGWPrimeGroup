import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("WPrimeAnalysis")
# get JSON file correctly parced
JSONfile = 'UserCode/CMGWPrimeGroup/JSON/json_160404-166011_DCSonly.txt'
myList = LumiList.LumiList (filename = JSONfile).getCMSSWString().split(',')

process.inputs = cms.PSet (
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    )
process.inputs.lumisToProcess.extend(myList)


process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
##    'file:patTuple.root'
  )
)

process.MessageLogger = cms.Service("MessageLogger")

process.WprimeAnalyzer = cms.PSet(
    ## common input for wrapped analyzers
    fileNames   = cms.vstring(),  ## keep empty!
   # fileNames   = cms.vstring('file:patTuple.root'),  ## mandatory
    outputFile  = cms.string('WZFakeRateMuon_analysis.root'),## mandatory
    maxEvents   = cms.int32(-1),                      ## optional
    reportAfter = cms.uint32(5000),                     ## optional
    useJSON = cms.bool(False),
    doRecoilCorrectionForW = cms.bool(False),
    sample_cross_sections = cms.string("samples_cross_sections_WZDilepton.txt"),
    debugme = cms.bool(False),
    preselect = cms.bool(False),
    logFile = cms.string("WZFakeRateMuon_event_counts.txt"),
    candEvtFile = cms.string("WZFakeRateMuon_CandEvts.txt"),
    ## enable analysis in individual channels
    runMuMETAnalysis = cms.bool(False),
    runElMETAnalysis = cms.bool(False),
    runWZAnalysis    = cms.bool(True),
    runHadVZAnalysis = cms.bool(False),
    runTBAnalysis    = cms.bool(False),
    runWgammaAnalysis = cms.bool(False),
    ## input specific for this analyzer
    muons = cms.string('userPatMuons'),
    electrons = cms.string('userPatElectrons'),
    met   = cms.string('patMETsPF'),
    particleFlow = cms.InputTag('selectedPatPFParticles'),
    genParticles = cms.InputTag('prunedGenParticles'),
    hltEventTag = cms.string('patTriggerEvent'),
    pileupTag  = cms.string('addPileupInfo'),
    inputs = process.inputs,

    MCPUDistFile = cms.string('UserCode/CMGWPrimeGroup/root_macros/MCPUDist.root'),
    MCPUDistHist = cms.string('pileup'),
    DataPUDistFile = cms.string('UserCode/CMGWPrimeGroup/root_macros/DataPUDist.root'),
    DataPUDistHist = cms.string('pileup'),

    useAdjustedMET = cms.bool(False),
    muonAlgo = cms.uint32(0),
    minDeltaR = cms.double(-999.),
    effectiveElecArea = cms.vdouble(0.0997,0.1123),#Not using Recommended PI*0.3*0.3
    effectiveMuonArea = cms.vdouble(0.1057,0.0769),
    triggersToUse = cms.vstring(#"HLT_Mu9",
                                #"HLT_Mu11",
                                #"HLT_Mu15_v*",
                                #'HLT_Mu17_v*',
                                #"HLT_Mu20_v*",
                                #"HLT_Mu24_v*",
                                #"HLT_Mu30_v*",

                                "HLT_Ele15_LW_L1R",
                                "HLT_Ele15_SW_L1R",
                                "HLT_Ele15_SW_CaloEleId_L1R",
                                "HLT_Ele17_SW_CaloEleId_L1R",
                                "HLT_Ele17_SW_TightEleId_L1R",
                                "HLT_Ele17_SW_TighterEleIdIsol_L1R_v*",
                                "HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*",
                                "HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*",
                                "HLT_Ele45_CaloIdVT_TrkIdT_v*",
                                ),


    Cuts = cms.vstring("NoCuts", 
                       "HLT",
                       
                       "ValidW", 
#                       "WFlavorMuon",
                       "WFlavorElec",

                       "FakeEvt",
                       "FakeLepTag",
                       "FakeLepProbeLoose",
                       "FakeLepProbeTight",
                       "EvtSetup",

                       "AllCuts"),
    LooseElecCuts = cms.vstring("ElecEta",
                                "ElecTightEt",
                                "ElecTightNMiss",
                                "ElecTightDistDCot",     
                                "ElecTightSigmaNN",
                                "ElecTightDeltaPhi",
                                "ElecTightDeltaEta",
                                ),
    TightElecCuts = cms.vstring("ElecEta",
                                "ElecTightEt",
                                "ElecTightNMiss",
                                "ElecTightDistDCot",     
                                "ElecTightSigmaNN",
                                "ElecTightDeltaPhi",
                                "ElecTightDeltaEta",
                                "ElecTightCombRelIso",
                                ),
    LooseMuonCuts = cms.vstring("MuonEta",
                                "MuonTightPt",
                                "MuonGlobal",
                                "MuonDxy",
                                "MuonNpxl",
                                "MuonNtrk",
                                "MuonNormChi2",
                                "MuonHitsUsed",
                                "MuonStations",
                                ),
    TightMuonCuts = cms.vstring("MuonEta",
                                "MuonTightPt",
                                "MuonGlobal",
                                "MuonDxy",
                                "MuonNpxl",
                                "MuonNtrk",
                                "MuonNormChi2",
                                "MuonHitsUsed",
                                "MuonStations",
                                "MuonIso",
                                ),

####################

    # +++++++++++++++++++General Cut values
    maxNumZs = cms.uint32(1),
    minNLeptons = cms.uint32(2),
    maxNLeptons = cms.uint32(3),
    minLeadPt = cms.double(35.),
    minMET = cms.double(30.),
    
    # +++++++++++++++++++Ht Cuts
    minHt = cms.double(190.),#150 for TC300), 190 for W'400
    
    # +++++++++++++++++++W Cuts
    minWtransMass = cms.double(0.),#Cory: Removed cut
    minWpt = cms.double(110.),#90 for TC300), 110 for W'400
    
    minWlepPt = cms.double(20.),

    cutWenuWPRelIsoMask = cms.int32(2),#Cory: Iso only
    cutElecWPTightType = cms.string("simpleEleId80relIso"),

    maxElecTightNMissingHits = cms.double(0.),
    minElecTightDist = cms.double(0.02),
    minElecTightDeltaCotTheta = cms.double(0.02),
    maxElecTightSigmaIetaIeta = cms.vdouble(0.01,0.031),
    maxElecTightDeltaPhi  = cms.vdouble(0.027,0.021),
    maxElecTightDeltaEta  = cms.vdouble(0.005,0.006),
    maxElecTightHOverE    = cms.vdouble(0.,0.),#Not used in 2011
#    maxElecTightCombRelIso = cms.vdouble(0.040,0.033),#2011 Rec
    maxElecTightCombRelIso = cms.vdouble(0.070,0.06), #2010 Rec

    maxMuonTightCombRelIso = cms.double(0.1),
      
    # +++++++++++++++++++Z Cuts
    minZeePt1 =  cms.double(20.),
    minZeePt2 =  cms.double(10.),
    minZmmPt1 =  cms.double(15.),
    minZmmPt2 =  cms.double(15.),
    minZpt =  cms.double(110.),#90 for TC300), 110 for W'400
    minZmass =  cms.double(60.),
    maxZmass =  cms.double(120.),
    
    # +++++++++++++++++++Electron General Cuts
    #VBTF Recommended Cuts
    minElecLooseEt = cms.double(10.),
    minElecTightEt = cms.double(20.),
    cutElecWPLooseMask = cms.int32(5),#Cory: No Iso
    cutElecWPLooseType = cms.string("simpleEleId95relIso"),
    
    maxElecNMissingHits = cms.double(0.),
    minElecDist = cms.double(0.),
    minElecDeltaCotTheta = cms.double(0.),
    maxElecSigmaIetaIeta = cms.vdouble(0.012,0.031),
    maxElecDeltaPhi  = cms.vdouble(0.800,0.7),
    maxElecDeltaEta  = cms.vdouble(0.007,0.011),
    maxElecHOverE    = cms.vdouble(0.,0.),#Not used in 2011
    maxElecCombRelIso = cms.vdouble(0.150,0.100),
    
    # +++++++++++++++++++Muon General Cuts
    maxMuonEta = cms.double(2.5),
    minMuonLoosePt = cms.double(15.),
    minMuonTightPt = cms.double(20.),
    #VBTF Recommended Cuts
    maxMuonDxy = cms.double(0.2),
    maxMuonNormChi2 = cms.double(10.),
    minMuonNPixHit = cms.int32(0),
    minMuonNTrkHit = cms.int32(10),
    minMuonStations = cms.int32(1),
    minMuonHitsUsed = cms.int32(0),
    maxMuonLooseCombRelIso = cms.double(0.15),

    )

