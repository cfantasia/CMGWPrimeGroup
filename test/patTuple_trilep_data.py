from UserCode.CMGWPrimeGroup.patTuple_trilep_cfg import *
from UserCode.CMGWPrimeGroup.patTuple_trigger_cfg import *


# 2nd argument: message-logger frequency
# 3rd argument: # of events to process
trilep_config(process, 100, 100, True)


from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.register ('globalTag',
                  "DONOTEXIST::All",
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.string,          # string, int, or float
                  "GlobalTag")
options.parseArguments()
process.GlobalTag.globaltag = options.globalTag

## remove MC matching from the default sequence when running on data
#removeMCMatching(process, ['All'])

addHLTFilter(process, 'HLT', "off")#For initial testing
#addHLTFilter(process, 'HLT', "singlemu")
#addHLTFilter(process, 'HLT', "singleelectron")
#addHLTFilter(process, 'HLT', "doublemu")
#addHLTFilter(process, 'HLT', "doubleelectron")

process.p.replace(process.patTrigger, process.patTrigger+process.hltFilter)

process.source.fileNames = [
    '/store/data/Run2012C/DoubleMu/AOD/PromptReco-v2/000/200/160/8457C3D4-14DF-E111-A67A-00237DDBE0E2.root'
    #'/store/data/Run2012C/DoubleElectron/AOD/PromptReco-v2/000/200/160/B4E531D5-14DF-E111-B75C-002481E0D73C.root'
    ] 

#process.out.outputCommands.append('keep *_*_*_*')

#print process.p

print process.GlobalTag.globaltag
