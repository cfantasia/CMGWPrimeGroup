from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *

from UserCode.CMGWPrimeGroup.patTuple_jetlep_cfg import *
from UserCode.CMGWPrimeGroup.patTuple_mc_cfg import *

# 2nd argument: message-logger frequency
# 3rd argument: # of events to process
jetlep_config(process, 1000, -1)
mc_config(process, cms)

process.p.replace(process.patMuons, process.muonMatch+process.patMuons+process.prunedGenParticles)
process.p.replace(process.patElectrons, process.electronMatch+process.patElectrons)

process.p.replace(process.pfNoPileUpIso, process.pfNoPileUpIso*process.newAK7PF)

process.patJetsAK7PF.jetSource = cms.InputTag("newAK7PF")


process.source.fileNames = [
    '/store/mc/Fall11/WprimeToWZTo2Q2L_M-1000_7TeV-pythia6-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/B8733CA2-DA28-E111-B278-00215E21D702.root'
#    'file:/home/fladias/WPrime/RSGravitonToZZToEEJJ_kMpl005_M_750_pythia6_cff_py_GEN_FASTSIM_HLT.root'
#    'file:/data/fladias/RSGravitonToZZToMuMuJJ_kMpl005_M_1000_pythia6_cff_py_GEN_FASTSIM_HLT.root'
#    'file:/afs/cern.ch/user/t/tomei/public/PYTHIA6_EXOTICA_RSGravZZ_kMpl005_M1000_7TeV_mumujj_cff_py_GEN_FASTSIM_HLT.root',
    ] 
