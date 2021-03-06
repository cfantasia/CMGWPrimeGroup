from UserCode.CMGWPrimeGroup.patTuple_common_cfg import *
from PhysicsTools.PatAlgos.tools.jetTools import *

def met_config(process) :
    ## select events with pfMET
    process.pfMETSelector = cms.EDFilter(
        "CandViewSelector",
        src = cms.InputTag("pfMet"),
        cut = cms.string( "pt > 10." ),
        )

    process.pfMETCounter = cms.EDFilter(
        "CandViewCountFilter",
        src = cms.InputTag("pfMETSelector"),
        minNumber = cms.uint32(1),
        )

    process.pfMETSeq = cms.Sequence(process.pfMETSelector*
                                    process.pfMETCounter
                                    )

    process.out.outputCommands.append('keep *_patMETsPF_*_*')
    process.out.outputCommands.append('keep *_pfMet_*_*')
    process.out.outputCommands.append('keep *_pfType1CorrectedMet_*_*')
    process.out.outputCommands.append('keep *_pfType1p2CorrectedMet_*_*')
    

 
