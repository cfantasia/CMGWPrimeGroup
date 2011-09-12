from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "Wprime_analysis.root"
process.WprimeAnalyzer.logFile     = "Wprime_event_counts.txt"
process.WprimeAnalyzer.candEvtFile = "Wprime_CandEvts.txt"

process.WprimeAnalyzer.Cuts = WprimeWZCuts

#process.WprimeAnalyzer.maxEvents   = 100
