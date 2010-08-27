import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo2")
process.source = cms.Source("PoolSource",
    # replace input file with one you want to use
 fileNames = cms.untracked.vstring('file:/home/cleonido/wprime/Summer09MC/W/F866EB68-E39C-DE11-8F01-00145EDD7879.root')
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
    pfMetTag = cms.InputTag("pfMet"),
    JetTag = cms.InputTag("iterativeCone5CaloJets"),
    TkIsoMapTag = cms.InputTag("muIsoDepositTk"),
    EcalIsoMapTag = cms.InputTag("muIsoDepositCalByAssociatorTowers","ecal"),
    HcalIsoMapTag = cms.InputTag("muIsoDepositCalByAssociatorTowers","hcal"),

    HLTriggerResults = cms.InputTag("TriggerResults","","HLT"),
     #	InputTag HLTriggerResults = "TriggerResults::HLT2"

    # muon-detector eta acceptance
    Detmu_acceptance = cms.double(2.4),

)

process.TFileService = cms.Service("TFileService",
    #       fileName = cms.string('TeVMuon_startupV4.root')
   fileName = cms.string('Wmu_31x_IDEAL_7TeV_1.root')
#   fileName = cms.string('wprime_1TeV.root')
)

process.p = cms.Path(process.StdMu)

