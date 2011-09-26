from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "WprimeWZ.root"
process.WprimeAnalyzer.logFile     = "WprimeWZ.dat"
process.WprimeAnalyzer.candEvtFile = "WprimeWZ.evt"

process.WprimeAnalyzer.Cuts = WprimeWZCuts

process.WprimeAnalyzer.minHt = 300. #W'600 Cut
