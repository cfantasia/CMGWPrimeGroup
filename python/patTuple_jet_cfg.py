from UserCode.CMGWPrimeGroup.patTuple_common_cfg import *

def jet_config(process, reportEveryNum=100, maxEvents=-1) :
    common_config(process, reportEveryNum, maxEvents)
    # RECO
    process.out.outputCommands.append('keep *_selectedPatJets_*_*')
    process.out.outputCommands.append('keep *_patJets_*_*')
    
def jetExtra_config(process, reportEveryNum=100, maxEvents=-1) :
    # RECO
    process.out.outputCommands.append('keep *_selectedPatJets_*_*')
    process.out.outputCommands.append('keep *_patJets_*_*')
    
