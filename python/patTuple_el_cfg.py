from UserCode.CMGWPrimeGroup.patTuple_common_cfg import *

def el_config(process) :
    # keep all events with electron-pt above 25 GeV
    process.lowPtElectronFilter = cms.EDFilter("CandViewCountFilter",
                                           src = cms.InputTag("selectedPatElectrons"),
                                           minNumber = cms.uint32(1),
                                           cut = cms.string( "pt > 10. & abs(eta) < 2.5" )
                                           )
    
    process.selectedPatPFParticles.cut = "abs(pdgId())==11"
    process.out.outputCommands.append('keep *_selectedPatElectrons_*_*')
    
