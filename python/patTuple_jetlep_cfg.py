from UserCode.CMGWPrimeGroup.patTuple_common_cfg import *
from UserCode.CMGWPrimeGroup.patTuple_el_cfg import *
from UserCode.CMGWPrimeGroup.patTuple_mu_cfg import *
from UserCode.CMGWPrimeGroup.patTuple_jet_cfg import *

def jetlep_config(process, reportEveryNum=100, maxEvents=-1) :
    process.load("UserCode.CMGWPrimeGroup.patTuple_jet_cfg")
    common_config(process, reportEveryNum, maxEvents)
    jet_config(process)
    el_config(process)
    mu_config(process)

    # keep all events with 2 leptons above 10 GeV
    process.countPatLeptons.electronSource = "selectedPatElectrons"
    process.countPatLeptons.muonSource     = "selectedPatMuons"
    process.countPatLeptons.minNumber = 2
    
    process.selectedPatElectrons.cut = "pt > 10. & abs(eta) < 2.5"
    process.selectedPatMuons.cut = "pt > 10. & abs(eta) < 2.4 & isGlobalMuon"
    
    # keep all events with jet-pt above 30 GeV, |eta| < 2.4
    process.selectedPatJets.cut = "pt > 30. & abs(eta) < 2.4"
    process.countPatJets.minNumber = 1
    process.countPatJets.src = "selectedPatJets"
    
    ## let it run
    process.p = cms.Path(
        process.patMuons *
        process.selectedPatMuons *
        process.patElectrons *
        process.selectedPatElectrons *
        process.countPatLeptons +
        process.goodOfflinePrimaryVertices*
        getattr(process,"PF2PATmod")* 
        #process.PF2PATmod *
        (process.patJetCorrFactors +
         process.patJets +
         process.selectedPatJets +
         process.countPatJets
         )*
        process.patTrigger
        )
    
    process.out.SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
        )

    from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
    process.kt6PFJetsPFlow = kt4PFJets.clone(
        rParam = cms.double(0.6),
        src = cms.InputTag('pfNoElectron'),
        doAreaFastjet = cms.bool(True),
        doRhoFastjet = cms.bool(True)
        )
    process.patJetCorrFactors.rho = cms.InputTag("kt6PFJetsPFlow", "rho")
    
    
    
    getattr(process,"PF2PATmod").replace(
        getattr(process,"pfNoElectron"),
        getattr(process,"pfNoElectron")*process.kt6PFJetsPFlow )
    
