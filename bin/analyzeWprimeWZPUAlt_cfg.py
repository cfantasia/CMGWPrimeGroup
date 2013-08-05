import sys
sys.path.append("UserCode/CMGWPrimeGroup/bin")
from analyzeWprimeWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "WprimeWZPUAlt.root"
process.WprimeAnalyzer.logFile     = "WprimeWZPUAlt.dat"
process.WprimeAnalyzer.candEvtFile = "WprimeWZPUAlt.evt"

process.WprimeAnalyzer.DataPUDistFile = 'UserCode/CMGWPrimeGroup/pileup/DataPileupHistogram-69p3mb.root'
