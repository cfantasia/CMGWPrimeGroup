from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = cms.string('EWKWZ_analysis.root')
process.WprimeAnalyzer.logFile = cms.string("EWKWZ_event_counts.txt")
process.WprimeAnalyzer.candEvtFile = cms.string("EWKWZ_CandEvts.txt")

process.WprimeAnalyzer.Cuts          = EWKWZCuts
