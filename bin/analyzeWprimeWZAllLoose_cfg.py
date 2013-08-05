import sys
sys.path.append("UserCode/CMGWPrimeGroup/bin")
from analyzeWprimeWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "WprimeWZAllLoose.root"
process.WprimeAnalyzer.logFile     = "WprimeWZAllLoose.dat"
process.WprimeAnalyzer.candEvtFile = "WprimeWZAllLoose.evt"

process.WprimeAnalyzer.LooseWElectronType = cms.untracked.string("CiC2012LoosePt20")
process.WprimeAnalyzer.TightWElectronType = cms.untracked.string("CiC2012LoosePt20")

