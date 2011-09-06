from UserCode.CMGWPrimeGroup.patTuple_jetlep_cfg import *

## remove MC matching from the default sequence when running on data
removeMCMatching(process, ['All'])

# 2nd argument: message-logger frequency
# 3rd argument: # of events to process
jetlep_config(process, 1000, -1)

print process.out.outputCommands

process.source.fileNames = [
    '/store/data/Run2011A/DoubleMu/AOD/PromptReco-v6/000/173/659/E2C7E831-B8CD-E011-9CA2-BCAEC5329725.root'
#    'file:/data/fladias/425_data_test.root' 
#    'file:/afs/cern.ch/user/t/tomei/tmp/V260/CMSSW_4_1_4/src/PYTHIA6_EXOTICA_RSGravZZ_kMpl005_M1000_7TeV_mumujj_cff_py_GEN_FASTSIM_HLT.root',
    ] 

