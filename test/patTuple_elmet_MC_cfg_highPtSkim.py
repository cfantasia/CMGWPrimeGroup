from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *

from UserCode.CMGWPrimeGroup.patTuple_elmet_cfg import *
from UserCode.CMGWPrimeGroup.patTuple_mc_cfg import *

# 2nd argument: message-logger frequency
# 3rd argument: # of events to process
elmet_config(process, 1000, -1)
mc_config(process, cms)

# keep all events with electron-pt above 100 GeV
process.selectedPatElectrons.cut = "pt > 100. & abs(eta) < 2.5"

process.highPtElectronFilter = cms.EDFilter("CandViewCountFilter",
                                src = cms.InputTag("selectedPatElectrons"),
                                minNumber = cms.uint32(1)
                                )

## let it run
process.p = cms.Path(
    process.patDefaultSequence *
    process.prunedGenParticles *
    process.highPtElectronFilter
)


process.source.fileNames = [          ##
    'rfio:/castor/cern.ch/user/c/cleonido/414_Wleptonic/WPlusToENu_CT10_TuneZ2_7TeV_powheg-pythia-4816882C-2D56-E011-B247-0017A4771030.root'
    ] 
#                                         ##
process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )

#process.out.fileName = '/tmp/cleonido/patTuple.root'
