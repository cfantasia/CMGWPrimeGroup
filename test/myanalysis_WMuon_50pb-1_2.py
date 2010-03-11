import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo2")
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
#    fileNames = cms.untracked.vstring('file:/home/cleonido/wprime/wprime_1.0/218_RECO_IDEAL/0055C643-B0B1-DD11-AEAB-001E0B47E400.root')
    fileNames = cms.untracked.vstring(
    "rfio:/castor/cern.ch/user/g/goys/Wprime/2_2_12/Wnew_1/skim/outfilter_2_1.root",
    "rfio:/castor/cern.ch/user/g/goys/Wprime/2_2_12/Wnew_1/skim/outfilter_2_2.root",
    "rfio:/castor/cern.ch/user/g/goys/Wprime/2_2_12/Wnew_1/skim/outfilter_2_3.root",
    "rfio:/castor/cern.ch/user/g/goys/Wprime/2_2_12/Wnew_1/skim/outfilter_2_4.root",
    "rfio:/castor/cern.ch/user/g/goys/Wprime/2_2_12/Wnew_1/skim/outfilter_2_5.root",
    "rfio:/castor/cern.ch/user/g/goys/Wprime/2_2_12/Wnew_1/skim/outfilter_2_6.root",
    "rfio:/castor/cern.ch/user/g/goys/Wprime/2_2_12/Wnew_1/skim/outfilter_2_7.root"
    )
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
    description = cms.string('W'),
    #     # of produced events (before filtering)
    Nprod_evt = cms.int32(2280000),

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
    # fileName = cms.string('wprime_1TeV.root')
 fileName = cms.string('/afs/cern.ch/user/c/cleonido/scratch0/Wmunu_2212_50pb-1_2.root')  
)

process.p = cms.Path(process.StdMu)

