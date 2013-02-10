import sys
sys.path.append("UserCode/CMGWPrimeGroup/bin")
from analyzeWprimeWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "WprimeWZMuPtRes.root"
process.WprimeAnalyzer.logFile     = "WprimeWZMuPtRes.dat"
process.WprimeAnalyzer.candEvtFile = "WprimeWZMuPtRes.evt"

# +++++++++++++++++++Options
process.WprimeAnalyzer.muScaleFactor = 1.006
