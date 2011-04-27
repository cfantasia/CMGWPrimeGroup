from PhysicsTools.PatAlgos.tools.pfTools import addPFCandidates

from UserCode.CMGWPrimeGroup.patTuple_common_cfg import *

def mumet_config(process, reportEveryNum=100, maxEvents=-1) :
    common_config(process, reportEveryNum, maxEvents)
    addPFCandidates(process, 'particleFlow')
    # this is needed so we can correct the pfMET by adjusting the muon-pt
    # when switching to one of the dedicated TeV muon reconstructors
    process.selectedPatPFParticles.cut = "particleId == 3"

    # RECO
    process.out.outputCommands.append('keep *_selectedPatMuons_*_*')
    process.out.outputCommands.append('keep *_patMETs*_*_*')
    process.out.outputCommands.append('keep *_selectedPatPFParticles*_*_*')
