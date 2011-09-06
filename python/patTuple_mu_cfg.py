from UserCode.CMGWPrimeGroup.patTuple_common_cfg import *

def mu_config(process) :
    process.lowPtMuonFilter = cms.EDFilter("CandViewCountFilter",
                                           src = cms.InputTag("selectedPatMuons"),
                                           minNumber = cms.uint32(1),
                                           cut = cms.string( "pt > 10. & abs(eta) < 2.5" )
                                           )
    
    process.selectedPatPFParticles.cut = "abs(pdgId())==13"
    process.out.outputCommands.append('keep *_selectedPatMuons_*_*')
    
