import sys
sys.path.append("UserCode/CMGWPrimeGroup/bin")
from analyzeWprimeWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "WprimeWZPUSys.root"
process.WprimeAnalyzer.logFile     = "WprimeWZPUSys.dat"
process.WprimeAnalyzer.candEvtFile = "WprimeWZPUSys.evt"

process.WprimeAnalyzer.DataPUDistFile = 'UserCode/CMGWPrimeGroup/pileup/DataPileupHistogram-73p5mb.root'
