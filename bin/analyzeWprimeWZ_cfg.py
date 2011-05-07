import FWCore.ParameterSet.Config as cms

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
    ## enable analysis in individual channels
    runMuMETAnalysis = cms.bool(False),
    runElMETAnalysis = cms.bool(False),
    runWZAnalysis    = cms.bool(True),
    runTBAnalysis    = cms.bool(False),
    runWgammaAnalysis = cms.bool(False),
    ## input specific for this analyzer
    muons = cms.InputTag('userPatMuons'),
    electrons = cms.InputTag('userPatElectrons'),
    met   = cms.InputTag('patMETsPF'),
    particleFlow = cms.InputTag('selectedPatPFParticles'),
    genParticles = cms.InputTag('prunedGenParticles'),
    #
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

