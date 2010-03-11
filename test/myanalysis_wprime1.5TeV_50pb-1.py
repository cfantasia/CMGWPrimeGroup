import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo2")
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring("file:/home/cleonido/wprime/wprime_1.5/2212_RECO_50PB-1/WPrimeReDigi2212M1500_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_50PBMU_V1_111.root",
                                      "file:/home/cleonido/wprime/wprime_1.5/2212_RECO_50PB-1/WPrimeReDigi2212M1500_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_50PBMU_V1_107.root",
                                      "file:/home/cleonido/wprime/wprime_1.5/2212_RECO_50PB-1/WPrimeReDigi2212M1500_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_50PBMU_V1_91.root",
                                      "file:/home/cleonido/wprime/wprime_1.5/2212_RECO_50PB-1/WPrimeReDigi2212M1500_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_50PBMU_V1_93.root",
                                      "file:/home/cleonido/wprime/wprime_1.5/2212_RECO_50PB-1/WPrimeReDigi2212M1500_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_50PBMU_V1_17.root",
                                      "file:/home/cleonido/wprime/wprime_1.5/2212_RECO_50PB-1/WPrimeReDigi2212M1500_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_50PBMU_V1_27.root"
                                      )
#    fileNames = cms.untracked.vstring(
#    "rfio:/castor/cern.ch/user/g/goys/Wprime/2_2_12/Wnew_1/skim/outfilter_1_1.root",
#    "rfio:/castor/cern.ch/user/g/goys/Wprime/2_2_12/Wnew_1/skim/outfilter_1_2.root",
#    )
)

# Number of events to process
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
#    input = cms.untracked.int32(100)    
)

# configuration for Wprime muon reconstruction
process.StdMu = cms.EDFilter("Wprime_muonreco",

    # input tags defined here
    MuonTag = cms.InputTag("muons"),
    MetTag = cms.InputTag("met"),
    JetTag = cms.InputTag("iterativeCone5CaloJets"),
    TkIsoMapTag = cms.InputTag("muIsoDepositTk"),
    EcalIsoMapTag = cms.InputTag("muIsoDepositCalByAssociatorTowers","ecal"),
    HcalIsoMapTag = cms.InputTag("muIsoDepositCalByAssociatorTowers","hcal"),

    HLTriggerResults = cms.InputTag("TriggerResults","","HLT"),
     #	InputTag HLTriggerResults = "TriggerResults::HLT2"

    # muon-detector eta acceptance
    Detmu_acceptance = cms.double(2.4),


    # sample description
    description = cms.string('Wprime 1.5 TeV'),
    #     # of produced events (before filtering)
    Nprod_evt = cms.int32(6000),

    # "golden" single-muon trigger name
    SingleMuHLT_20x = cms.string('HLT1MuonNonIso9'),
    SingleMuL1 = cms.string('HLT_L1Mu'),
    SingleMuHLT_21x = cms.string('HLT_Mu9'),

    # generic trigger name containing single muons
    AnyMuHLT_20x = cms.string('HLT1Muon'),
    AnyMuHLT_21x = cms.string('Mu'),

)

process.TFileService = cms.Service("TFileService",
    #       fileName = cms.string('TeVMuon_startupV4.root')
     fileName = cms.string('/home/cleonido/wprime/v66_50pb-1_hack/wprime_1.5TeV.root')
#  fileName = cms.string('/afs/cern.ch/user/c/cleonido/scratch0/w_50pb-1_1.root')  
)

process.p = cms.Path(process.StdMu)
