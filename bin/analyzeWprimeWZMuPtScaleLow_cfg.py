import sys
sys.path.append("UserCode/CMGWPrimeGroup/bin")
from analyzeWprimeWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "WprimeWZMuPtScaleLow.root"
process.WprimeAnalyzer.logFile     = "WprimeWZMuPtScaleLow.dat"
process.WprimeAnalyzer.candEvtFile = "WprimeWZMuPtScaleLow.evt"

# +++++++++++++++++++Options
process.WprimeAnalyzer.muScaleFactor = 1.002
