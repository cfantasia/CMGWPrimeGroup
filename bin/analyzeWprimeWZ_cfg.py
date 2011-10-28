from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "WprimeWZ.root"
process.WprimeAnalyzer.logFile     = "WprimeWZ.dat"
process.WprimeAnalyzer.candEvtFile = "WprimeWZ.evt"

process.WprimeAnalyzer.Cuts = WprimeWZCuts

process.WprimeAnalyzer.minMET = 50.

# +++++++++++++++++++Analysis Cuts
process.WprimeAnalyzer.minHt = cms.double(300.) #W'600 Cut
process.WprimeAnalyzer.minZpt =  cms.untracked.double(0.)
process.WprimeAnalyzer.minWpt = cms.untracked.double(0.)
