from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "EWKWZ.root"
process.WprimeAnalyzer.logFile     = "EWKWZ.dat"
process.WprimeAnalyzer.candEvtFile = "EWKWZ.evt"

process.WprimeAnalyzer.Cuts          = EWKWZCuts
