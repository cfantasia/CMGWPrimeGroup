import sys
sys.path.append("UserCode/CMGWPrimeGroup/bin")
from analyzeWprimeWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "WprimeWZMuPtScale.root"
process.WprimeAnalyzer.logFile     = "WprimeWZMuPtScale.dat"
process.WprimeAnalyzer.candEvtFile = "WprimeWZMuPtScale.evt"

# +++++++++++++++++++Options
process.WprimeAnalyzer.muScaleFactor = 1.002
