from UserCode.CMGWPrimeGroup.patTuple_common_cfg import *

def mu_config(process) :

    process.patMuons.embedTrack = True
    process.patMuons.embedPFCandidate = True

    # event content to include all muons within |eta|<2.4 with pt>20
    process.selectedPatMuons.cut = "pt > 20. & abs(eta) < 2.4"

    process.selectedPatPFParticles.cut = "abs(pdgId())==13"
    process.out.outputCommands.append('keep *_selectedPatMuons_*_*')
    
    # define filter with muon-pt above 100 GeV
    # NB: it is selectedPatMuons that are saved in the event!
    # NB: highPtMuons is only used for filtering when invoked
    process.highPtMuons = process.selectedPatMuons.clone()
    process.highPtMuons.cut = "pt > 100. & abs(eta) < 2.4"
    # comment (CL): we should probably add a isGlobal criterion here
    # in a high-pt skim, about 80% do not have a valid global muon track...
    
    process.highPtMuonFilter = cms.EDFilter("CandViewCountFilter",
                                            src = cms.InputTag("highPtMuons"),
                                            minNumber = cms.uint32(1)
                                            )
