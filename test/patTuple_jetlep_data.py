from UserCode.CMGWPrimeGroup.patTuple_jetlep_cfg import *
from UserCode.CMGWPrimeGroup.patTuple_trigger_cfg import *

## remove MC matching from the default sequence when running on data
removeMCMatching(process, ['All'])

# 2nd argument: message-logger frequency
# 3rd argument: # of events to process
jetlep_config(process, 1000, 100)


#addHLTFilter(process, 'HLT', "singlemu")
#addHLTFilter(process, 'HLT', "singleelectron")
addHLTFilter(process, 'HLT', "doublemu")
#addHLTFilter(process, 'HLT', "doubleelectron")

process.patJetCorrFactorsAK7PF.levels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual')
process.patJetCorrFactorsAK5PF.levels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual')


process.p.replace(process.patTrigger, process.patTrigger+process.hltFilter)

process.source.fileNames = [
    '/store/data/Run2011A/DoubleMu/AOD/PromptReco-v6/000/173/659/E2C7E831-B8CD-E011-9CA2-BCAEC5329725.root'
#    '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v6/000/173/659/4645CC17-B8CD-E011-8BA0-001D09F2906A.root'
#    'file:/data/fladias/425_data_test.root' 
#    'file:/afs/cern.ch/user/t/tomei/tmp/V260/CMSSW_4_1_4/src/PYTHIA6_EXOTICA_RSGravZZ_kMpl005_M1000_7TeV_mumujj_cff_py_GEN_FASTSIM_HLT.root',
    ] 

