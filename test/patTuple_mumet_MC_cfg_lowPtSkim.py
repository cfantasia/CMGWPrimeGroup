from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *

from UserCode.CMGWPrimeGroup.patTuple_mumet_cfg import *
from UserCode.CMGWPrimeGroup.patTuple_mc_cfg import *

# 2nd argument: message-logger frequency
# 3rd argument: # of events to process
mumet_config(process, 1000, -1)
mc_config(process, cms)

# keep all events with muon-pt above 25 GeV, but less than 100 GeV
process.lowPtMuons = process.selectedPatMuons.clone()
process.lowPtMuons.cut = "pt > 25. & abs(eta) < 2.4"

process.lowPtMuonFilter = cms.EDFilter("CandViewCountFilter",
                                src = cms.InputTag("lowPtMuons"),
                                minNumber = cms.uint32(1)
                                )

## let it run
process.p = cms.Path(
    process.patDefaultSequence *
    process.lowPtMuons *
    process.highPtMuons *
    process.prunedGenParticles *
    process.lowPtMuonFilter *
    ~process.highPtMuonFilter
)


process.source.fileNames = [          ##
    'rfio:/castor/cern.ch/user/c/cleonido/414_Wleptonic/WPlusToMuNu_CT10_TuneZ2_7TeV_powheg-pythia-FE7FA83F-CA5B-E011-BC09-90E6BA0D09B6.root'
    ]                          
#                                         ##
process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )

#process.out.fileName = '/tmp/cleonido/patTuple.root'

