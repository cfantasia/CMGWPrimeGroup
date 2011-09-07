from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *

from UserCode.CMGWPrimeGroup.patTuple_elmet_cfg import *

## remove MC matching from the default sequence when running on data
removeMCMatching(process, ['All'])

# 2nd argument: message-logger frequency
# 3rd argument: # of events to process
elmet_config(process, 1000, -1)

# keep all events with electron-pt above 100 GeV
process.selectedPatElectrons.cut = "pt > 100. & abs(eta) < 2.5"

process.lowPtElectronFilter = cms.EDFilter("CandViewCountFilter",
                                src = cms.InputTag("selectedPatElectrons"),
                                minNumber = cms.uint32(1)
                                )

## let it run
process.p = cms.Path(
    process.patDefaultSequence *
    process.lowPtElectronFilter
)

#                                         ##
process.source.fileNames = [          ##
#    '/store/data/Run2011A/SingleElectron/AOD/PromptReco-v1/000/161/312/90646AF9-F957-E011-B0DB-003048F118C4.root'
    'file:/tmp/cleonido/PromptRecoV6/688EA5AA-50C7-E011-831C-003048F1182E.root'
    ] 
#                                         ##

process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )

#process.out.fileName = '/tmp/cleonido/patTuple_ElMET.root'
