from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "WprimeWZ_analysis.root"
process.WprimeAnalyzer.logFile     = "WprimeWZ_event_counts.txt"
process.WprimeAnalyzer.candEvtFile = "WprimeWZ_CandEvts.txt"

process.WprimeAnalyzer.Cuts = WprimeWZCuts

process.WprimeAnalyzer.minHt = 300. #W'600 Cut

#process.WprimeAnalyzer.maxEvents   = 100
