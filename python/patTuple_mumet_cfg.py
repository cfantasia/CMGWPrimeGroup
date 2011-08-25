from UserCode.CMGWPrimeGroup.patTuple_common_cfg import *

def mumet_config(process, reportEveryNum=100, maxEvents=-1) :
    common_config(process, reportEveryNum, maxEvents)
    process.selectedPatPFParticles.cut = "abs(pdgId())==13"

    # RECO
    process.out.outputCommands.append('keep *_selectedPatMuons_*_*')
    process.out.outputCommands.append('keep *_patMETs*_*_*')
