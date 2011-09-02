from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *

from UserCode.CMGWPrimeGroup.patTuple_mumet_cfg import *

## remove MC matching from the default sequence when running on data
removeMCMatching(process, ['All'])

# 2nd argument: message-logger frequency
# 3rd argument: # of events to process
mumet_config(process, 1000, 10)

# keep all events with muon-pt above 25 GeV
process.selectedPatMuons.cut = "pt > 25. & abs(eta) < 2.5"

process.lowPtMuonFilter = cms.EDFilter("CandViewCountFilter",
                                src = cms.InputTag("selectedPatMuons"),
                                minNumber = cms.uint32(1)
                                )

## let it run
process.p = cms.Path(
    process.patDefaultSequence *
    process.lowPtMuonFilter
)

#                                         ##
process.source.fileNames = [          ##
#    '/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/161/312/F8AEC745-DF57-E011-8D23-001D09F290BF.root'
    'file:/tmp/cleonido/PromptRecoV6/7AE2B739-30C7-E011-9C8A-003048F117B4.root'
#    '/store/relval/CMSSW_3_8_6/RelValTTbar/GEN-SIM-RECO/START38_V13-v1/0068/98EA8C65-25E8-DF11-B0A0-0018F3D095F8.root'
    ]                          
#                                         ##
process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )

process.out.fileName = '/tmp/cleonido/patTuple_MuMET.root'

