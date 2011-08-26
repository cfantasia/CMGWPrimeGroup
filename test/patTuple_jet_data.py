from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *

process.load("UserCode.CMGWPrimeGroup.patTuple_jet_cfg")
from UserCode.CMGWPrimeGroup.patTuple_jet_cfg import jetExtra_config
from UserCode.CMGWPrimeGroup.patTuple_mumet_cfg import *

## remove MC matching from the default sequence when running on data
removeMCMatching(process, ['All'])

# 2nd argument: message-logger frequency
# 3rd argument: # of events to process
mumet_config(process, 1000, -1)
jetExtra_config(process, 1000, -1)

# keep all events with 2 global muons above 10 GeV, |eta| < 2.4
process.selectedPatMuons.cut = "pt > 10. & abs(eta) < 2.4 & isGlobalMuon"

process.lowPtMuonFilter = cms.EDFilter("CandViewCountFilter",
                                src = cms.InputTag("selectedPatMuons"),
                                minNumber = cms.uint32(2)
                                )

# keep all events with jet-pt above 30 GeV, |eta| < 2.4
process.selectedPatJets.cut = "pt > 30. & abs(eta) < 2.4"

process.lowPtJetFilter = cms.EDFilter("CandViewCountFilter",
                                src = cms.InputTag("selectedPatJets"),
                                minNumber = cms.uint32(1)
                                )


from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFJetsPFlow = kt4PFJets.clone(
        rParam = cms.double(0.6),
            src = cms.InputTag('pfNoElectron'),
            doAreaFastjet = cms.bool(True),
            doRhoFastjet = cms.bool(True)
            )
process.patJetCorrFactors.rho = cms.InputTag("kt6PFJetsPFlow", "rho")



getattr(process,"PF2PATmod").replace(
        getattr(process,"pfNoElectron"),
            getattr(process,"pfNoElectron")*process.kt6PFJetsPFlow )


## let it run
process.p = cms.Path(
#    process.muonMatch +
    process.patMuons *
    process.selectedPatMuons *
    process.lowPtMuonFilter +
    process.goodOfflinePrimaryVertices*
    getattr(process,"PF2PATmod")* 
    #process.PF2PATmod *
    (process.patJetCorrFactors +
     process.patJets +
     process.selectedPatJets +
     process.lowPtJetFilter
     )
)

process.source.fileNames = [
    '/store/data/Run2011A/SingleMu/AOD/PromptReco-v6/000/173/212/105751F0-8AC6-E011-AC0F-0030487C8CB6.root'
#    'file:/data/fladias/425_data_test.root' 
#    'file:/afs/cern.ch/user/t/tomei/tmp/V260/CMSSW_4_1_4/src/PYTHIA6_EXOTICA_RSGravZZ_kMpl005_M1000_7TeV_mumujj_cff_py_GEN_FASTSIM_HLT.root',
    ] 

process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )
