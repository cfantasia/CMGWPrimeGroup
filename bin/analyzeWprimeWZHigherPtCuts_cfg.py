import sys
sys.path.append("UserCode/CMGWPrimeGroup/bin")
from analyzeWprimeWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "WprimeWZHigherPtCuts.root"
process.WprimeAnalyzer.logFile     = "WprimeWZHigherPtCuts.dat"
process.WprimeAnalyzer.candEvtFile = "WprimeWZHigherPtCuts.evt"

process.WprimeAnalyzer.minZmmPt1 =  35.
process.WprimeAnalyzer.minZmmPt2 =  35.
