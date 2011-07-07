from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = cms.string('Wprime_analysis.root')
process.WprimeAnalyzer.logFile = cms.string("Wprime_event_counts.txt")
process.WprimeAnalyzer.candEvtFile = cms.string("Wprime_CandEvts.txt")

process.WprimeAnalyzer.Cuts = WprimeWZCuts
