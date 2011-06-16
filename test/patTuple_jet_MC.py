from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *

process.load("UserCode.CMGWPrimeGroup.patTuple_jet_cfg")
from UserCode.CMGWPrimeGroup.patTuple_jet_cfg import jetExtra_config
from UserCode.CMGWPrimeGroup.patTuple_mumet_cfg import *
from UserCode.CMGWPrimeGroup.patTuple_mc_cfg import *

## remove MC matching from the default sequence when running on data
#removeMCMatching(process, ['All'])

# 2nd argument: message-logger frequency
# 3rd argument: # of events to process
mumet_config(process, 1, 5000)
jetExtra_config(process, 1, 5000)
mc_config(process, cms)

# keep all events with jet-pt above 30 GeV, |eta| < 2.4
process.selectedPatJets.cut = "pt > 30. & abs(eta) < 2.4"
process.lowPtJetFilter = cms.EDFilter("CandViewCountFilter",
                                src = cms.InputTag("selectedPatJets"),
                                minNumber = cms.uint32(1)
                                )

## let it run
process.p = cms.Path(
    process.muonMatch + 
    process.patMuons + 
    process.selectedPatMuons +
    process.PF2PATmod *
    (process.patJetCorrFactors + 
     process.patJets +
     process.selectedPatJets +
     process.lowPtJetFilter)
)

process.source.fileNames = [
    'file:/afs/cern.ch/user/t/tomei/public/PYTHIA6_EXOTICA_RSGravZZ_kMpl005_M1000_7TeV_mumujj_cff_py_GEN_FASTSIM_HLT.root',
    ] 

process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )
