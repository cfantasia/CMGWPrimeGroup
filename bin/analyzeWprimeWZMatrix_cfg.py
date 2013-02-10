import sys
sys.path.append("UserCode/CMGWPrimeGroup/bin")
from analyzeWprimeWZ_cfg import *

process.WprimeAnalyzer.outputFile  = 'WprimeWZMatrix.root'
process.WprimeAnalyzer.logFile     = "WprimeWZMatrix.dat"
process.WprimeAnalyzer.candEvtFile = "WprimeWZMatrix.evt"

process.WprimeAnalyzer.LooseElectronType = "WZLoose"
process.WprimeAnalyzer.TightElectronType = "WZRelaxed"
process.WprimeAnalyzer.LooseMuonType = "WZLoose"
process.WprimeAnalyzer.TightMuonType = "WZRelaxed"

# +++++++++++++++++++Options
process.WprimeAnalyzer.minNLeptons = cms.untracked.uint32(2)
