import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo2")

#The global tag is needed if the option extractL1Prescales is set to True
#Just switch to False if this is not needed or if there is some problem
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR10_H_V9::All'
process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.limit = 10
process.MessageLogger.cerr.FwkReport.reportEvery = 100


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
#    fileNames = cms.untracked.vstring('file:/home/cleonido/wprime/wprime_1.0/218_RECO_IDEAL/0055C643-B0B1-DD11-AEAB-001E0B47E400.root')
    fileNames = cms.untracked.vstring(
    "rfio:/castor/cern.ch/user/g/goys/Wprime/2_2_12/Wnew_1/skim/outfilter_1_1.root",
    "rfio:/castor/cern.ch/user/g/goys/Wprime/2_2_12/Wnew_1/skim/outfilter_1_2.root",
    "rfio:/castor/cern.ch/user/g/goys/Wprime/2_2_12/Wnew_1/skim/outfilter_1_3.root",
    "rfio:/castor/cern.ch/user/g/goys/Wprime/2_2_12/Wnew_1/skim/outfilter_1_4.root",
    "rfio:/castor/cern.ch/user/g/goys/Wprime/2_2_12/Wnew_1/skim/outfilter_1_5.root",
    "rfio:/castor/cern.ch/user/g/goys/Wprime/2_2_12/Wnew_1/skim/outfilter_1_6.root",
    "rfio:/castor/cern.ch/user/g/goys/Wprime/2_2_12/Wnew_1/skim/outfilter_1_7.root",
    "rfio:/castor/cern.ch/user/g/goys/Wprime/2_2_12/Wnew_1/skim/outfilter_1_8.root",
    "rfio:/castor/cern.ch/user/g/goys/Wprime/2_2_12/Wnew_1/skim/outfilter_1_9.root"
    )
)

# Number of events to process
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
#    input = cms.untracked.int32(100)    
)

# configuration for Wprime muon reconstruction
process.primaryVertexFilter = cms.EDFilter("VertexSelector",
                                           src = cms.InputTag("offlinePrimaryVertices"),
                                           cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"), # tracksSize() > 3 for the older cut
                                           filter = cms.bool(True),
                                           # otherwise it won't filter the events, just produce an empty vertex collection.
                                           )

process.StdMu = cms.EDAnalyzer("Wprime_muonreco",

    # input tags defined here
    pvTag = cms.InputTag("offlinePrimaryVertices"),
    pvBSTag = cms.InputTag("offlinePrimaryVerticesWithBS"),
    MuonTag = cms.InputTag("muons"),
    tevMuonLabel = cms.string("tevMuons"),
    pfMetTag = cms.InputTag("pfMet"),
    pfJetTag = cms.InputTag("ak5PFJets"),
    caloJetTag = cms.InputTag("ak5CaloJets"),
    TkIsoMapTag = cms.InputTag("muIsoDepositTk"),
    EcalIsoMapTag = cms.InputTag("muIsoDepositCalByAssociatorTowers","ecal"),
    HcalIsoMapTag = cms.InputTag("muIsoDepositCalByAssociatorTowers","hcal"),

    HLTriggerResults = cms.InputTag("TriggerResults","","HLT"),
     #	InputTag HLTriggerResults = "TriggerResults::HLT2"
    #this needs a global tag above.
    #By default, the HLT prescales are extracted
    extractL1Prescales = cms.bool(True),
    triggerConditions = cms.vstring('HLT_L1MuOpen',
                                    'HLT_L1MuOpen_v*',
                                    'HLT_L1Mu',
                                    'HLT_L1Mu_v*',
                                    'HLT_L2Mu5',
                                    'HLT_L2Mu5_v*',
                                    'HLT_L2Mu9',
                                    'HLT_L2Mu9_v*',
                                    'HLT_L2Mu11',
                                    'HLT_L2Mu11_v*',
                                    'HLT_L2Mu15',
                                    'HLT_L2Mu15_v*',
                                    'HLT_L2Mu25',
                                    'HLT_L2Mu25_v*',
                                    'HLT_Mu3',
                                    'HLT_Mu3_v*',
                                    'HLT_Mu5',
                                    'HLT_Mu5_v*',
                                    'HLT_Mu7',
                                    'HLT_Mu7_v*',
                                    'HLT_Mu9',
                                    'HLT_Mu9_v*',
                                    'HLT_Mu11',
                                    'HLT_Mu11_v*',
                                    'HLT_Mu13',
                                    'HLT_Mu13_v*',
                                    'HLT_Mu15',
                                    'HLT_Mu15_v*',
                                    'HLT_Mu17',
                                    'HLT_Mu17_v*',
                                    'HLT_Mu19',
                                    'HLT_Mu19_v*',
                                    'HLT_Mu21',
                                    'HLT_Mu21_v*',
                                    'HLT_IsoMu9',
                                    'HLT_IsoMu9_v*',
                                    'HLT_IsoMu11',
                                    'HLT_IsoMu11_v*',
                                    'HLT_IsoMu13',
                                    'HLT_IsoMu13_v*',
                                    'HLT_IsoMu15',
                                    'HLT_IsoMu15_v*',
                                    'HLT_IsoMu17',
                                    'HLT_IsoMu17_v*',
                                    'HLT_IsoMu19',
                                    'HLT_IsoMu19_v*',
                                    'HLT_IsoMu21',
                                    'HLT_IsoMu21_v*',
                                    ),   

    # muon-detector eta acceptance
    Detmu_acceptance = cms.double(2.4),


    # sample description
    description = cms.string('W'),

)

process.TFileService = cms.Service("TFileService",
    #       fileName = cms.string('TeVMuon_startupV4.root')
    # fileName = cms.string('wprime_1TeV.root')
 fileName = cms.string('/afs/cern.ch/user/c/cleonido/scratch0/Wmunu_2212_50pb-1_1.root')  
)

process.p = cms.Path(process.StdMu)

