import sys
sys.path.append("UserCode/CMGWPrimeGroup/bin")
from analyzeWprimeWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "WprimeWZTrilep.root"
process.WprimeAnalyzer.logFile     = "WprimeWZTrilep.dat"
process.WprimeAnalyzer.candEvtFile = "WprimeWZTrilep.evt"

process.WprimeAnalyzer.minDeltaR = cms.double(0.)

