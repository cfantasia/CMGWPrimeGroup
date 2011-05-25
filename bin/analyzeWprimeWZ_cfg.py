import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("WPrimeAnalysis")
# get JSON file correctly parced
#JSONfile = 'UserCode/CMGWPrimeGroup/JSON/Cert_160404-163757_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt'
#JSONfile = 'UserCode/CMGWPrimeGroup/JSON/json_160404-163869_DCSonly.txt'
#myList = LumiList.LumiList (filename = JSONfile).getCMSSWString().split(',')

#process.inputs = cms.PSet (
#    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
#    )
#process.inputs.lumisToProcess.extend(myList)


process = cms.Process("WPrimeAnalysis")

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
    outputFile  = cms.string('Wprime_analysis.root'),## mandatory
    maxEvents   = cms.int32(-1),                      ## optional
    reportAfter = cms.uint32(1000),                     ## optional
    doRecoilCorrectionForW = cms.bool(False),
    sample_cross_sections = cms.string("samples_cross_sections_WZ.txt"),
    debugme = cms.bool(False),
    preselect = cms.bool(True),
    LogFile = cms.string("Wprime_event_counts.txt"),
    CandEvtFile = cms.string("Wprime_CandEvts.txt"),
    ## enable analysis in individual channels
    runMuMETAnalysis = cms.bool(False),
    runElMETAnalysis = cms.bool(False),
    runWZAnalysis    = cms.bool(True),
    runTBAnalysis    = cms.bool(False),
    runWgammaAnalysis = cms.bool(False),
    ## input specific for this analyzer
    muons = cms.string('userPatMuons'),
    electrons = cms.string('userPatElectrons'),
    met   = cms.string('patMETsPF'),
    particleFlow = cms.string('selectedPatPFParticles'),
    genParticles = cms.string('prunedGenParticles'),
    hltEventTag = cms.string('patTriggerEvent'),
    pileupTag  = cms.string('addPileupInfo'),
    inputs = process.inputs,

    muonAlgo = cms.int32(0),
    #

    triggersToUse = cms.vstring("HLT_Mu9",
                                "HLT_Mu11",
                                "HLT_Mu15_v?",
                                "HLT_Mu20_v?",
                                "HLT_Mu24_v?",
                                "HLT_Mu30_v?",
                                
                                "HLT_Ele15_LW_L1R",
                                "HLT_Ele15_SW_L1R",
                                "HLT_Ele15_SW_CaloEleId_L1R",
                                "HLT_Ele17_SW_CaloEleId_L1R",
                                "HLT_Ele17_SW_TightEleId_L1R",
                                "HLT_Ele17_SW_TighterEleIdIsol_L1R_v?",
                                "HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v?",
                                "HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v?",
                                "HLT_Ele45_CaloIdVT_TrkIdT_v?",

                                ),


    Cuts = cms.vstring("NoCuts", 
                       #"HLT", 
                       "ValidWandZ", 
                       "ValidWZCand",
                       "NumZs", 
                       "ZMass", 
                       #"WTransMass", 
                       "MET",
                       "Ht", 
                       "Zpt", 
                       "Wpt",
                       "AllCuts"),
    LooseElecCuts = cms.vstring("ElecLooseEt", 
                                "ElecLooseID"),
    TightElecCuts = cms.vstring("ElecLoose",
                                "ElecTightEt",
                                "ElecIso"),
    LooseMuonCuts = cms.vstring("MuonLoosePt", 
                                "MuonEta",
                                "MuonGlobal",
                                "MuonDxy",
                                "MuonNpxl",
                                "MuonNtrk",
                                "MuonNormChi2",
                                "MuonHitsUsed",
                                "MuonStations",
                                ),
    TightMuonCuts = cms.vstring("MuonLoose",
                                "MuonTightPt",
                                "MuonIso"),

####################

# +++++++++++++++++++General Cut values
  maxNumZs = cms.int32(2),
  minNumLeptons = cms.int32(3),
  minMET = cms.double(30.),

# +++++++++++++++++++Ht Cuts
  minHt = cms.double(190.),#150 for TC300), 190 for W'400

# +++++++++++++++++++W Cuts
  minWtransMass = cms.double(0.),#Cory: Removed cut
  minWpt = cms.double(110.),#90 for TC300), 110 for W'400

  maxWmunuCombRelIso = cms.double(0.15),

  cutWenuWPRelIsoMask = cms.int32(2),#Cory: Iso only

  maxWenuTrkRelIso   = cms.vdouble(0.30,0.20),
  maxWenuECalRelIso  = cms.vdouble(0.20,0.15),
  maxWenuHCalRelIso  = cms.vdouble(0.15,0.12),

# +++++++++++++++++++Z Cuts
  minZpt =  cms.double(110.),#90 for TC300), 110 for W'400
  minZmass =  cms.double(80.),
  maxZmass =  cms.double(100.),

# +++++++++++++++++++Electron General Cuts
#VBTF Recommended Cuts
  minElecLooseEt = cms.double(10.),
  minElecTightEt = cms.double(20.),
  cutElecWPLooseMask = cms.int32(5),#Cory: No Iso

  maxElecSigmaiEtaiEta = cms.vdouble(0.01,0.03),
  maxElecDeltaPhiIn  = cms.vdouble(0.08,0.7),
  maxElecDeltaEtaIn  = cms.vdouble(0.007,0.01),
  maxElecHOverE      = cms.vdouble(0.15,0.07),

# +++++++++++++++++++Muon General Cuts
  maxMuonEta = cms.double(2.5),
  minMuonLoosePt = cms.double(10.),
  minMuonTightPt = cms.double(20.),
#VBTF Recommended Cuts
  maxMuonDxy = cms.double(0.2),
  maxMuonNormChi2 = cms.double(10.),
  minMuonNPixHit = cms.int32(0),
  minMuonNTrkHit = cms.int32(10),
  minMuonStations = cms.int32(1),
  minMuonHitsUsed = cms.int32(0),
)

