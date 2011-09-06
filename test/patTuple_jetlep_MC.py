from UserCode.CMGWPrimeGroup.patTuple_jetlep_cfg import *
from UserCode.CMGWPrimeGroup.patTuple_mc_cfg import *

# 2nd argument: message-logger frequency
# 3rd argument: # of events to process
jetlep_config(process, 1000, -1)
mc_config(process, cms)

process.p.replace(process.patMuons, process.muonMatch+process.patMuons)
process.p.replace(process.patElectrons, process.electronMatch+process.patElectrons)

process.source.fileNames = [
    'file:/data/fladias/RSGravitonToZZToMuMuJJ_kMpl005_M_1000_pythia6_cff_py_GEN_FASTSIM_HLT.root'
#    'file:/afs/cern.ch/user/t/tomei/public/PYTHIA6_EXOTICA_RSGravZZ_kMpl005_M1000_7TeV_mumujj_cff_py_GEN_FASTSIM_HLT.root',
    ] 
