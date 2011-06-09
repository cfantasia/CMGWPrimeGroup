from UserCode.CMGWPrimeGroup.patTuple_common_cfg import *

def jet_config(process, reportEveryNum=100, maxEvents=-1) :
    common_config(process, reportEveryNum, maxEvents)
    # RECO
    process.out.outputCommands.append('keep patJets_*_*_*')
    

