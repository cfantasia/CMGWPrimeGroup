from UserCode.CMGWPrimeGroup.patTuple_common_cfg import *

def elmet_config(process, reportEveryNum=100, maxEvents=-1) :
    common_config(process, reportEveryNum, maxEvents)
    # RECO
    process.out.outputCommands.append('keep *_selectedPatElectrons_*_*')
    process.out.outputCommands.append('keep *_patMETs*_*_*')