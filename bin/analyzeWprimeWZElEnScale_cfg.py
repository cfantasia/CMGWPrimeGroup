import sys
sys.path.append("UserCode/CMGWPrimeGroup/bin")
from analyzeWprimeWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "WprimeWZElEnScale.root"
process.WprimeAnalyzer.logFile     = "WprimeWZElEnScale.dat"
process.WprimeAnalyzer.candEvtFile = "WprimeWZElEnScale.evt"

# +++++++++++++++++++Options
process.WprimeAnalyzer.elScaleFactor = 1.02
