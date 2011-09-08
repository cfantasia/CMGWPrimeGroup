from UserCode.CMGWPrimeGroup.patTuple_common_cfg import *

def el_config(process) :
    # event content to include all electrons within |eta|<2.5 with pt>20
    process.selectedPatElectrons.cut = "pt > 20. & abs(eta) < 2.5"
    
    process.selectedPatPFParticles.cut = "abs(pdgId())==11"
    process.out.outputCommands.append('keep *_selectedPatElectrons_*_*')
    
    # define filter with electron-pt above 100 GeV
    # NB: it is selectedPatElectrons that are saved in the event!
    # NB: highPtElectrons is only used for filtering when invoked
    process.highPtElectrons = process.selectedPatElectrons.clone()
    process.highPtElectrons.cut = "pt > 100. & abs(eta) < 2.5"
    
    process.highPtMuonFilter = cms.EDFilter("CandViewCountFilter",
                                            src = cms.InputTag("highPtElectrons"),
                                            minNumber = cms.uint32(1)
                                            )
