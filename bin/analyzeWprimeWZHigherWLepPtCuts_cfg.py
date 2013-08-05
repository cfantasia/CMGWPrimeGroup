import sys
sys.path.append("UserCode/CMGWPrimeGroup/bin")
from analyzeWprimeWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "WprimeWZHigherWLepPtCuts.root"
process.WprimeAnalyzer.logFile     = "WprimeWZHigherWLepPtCuts.dat"
process.WprimeAnalyzer.candEvtFile = "WprimeWZHigherWLepPtCuts.evt"

process.WprimeAnalyzer.minZmmPt1 =  35.
process.WprimeAnalyzer.minZmmPt2 =  35.

process.WprimeAnalyzer.electronSelectors.CiC2012MediumPt20.joint.minPt = 35.
process.WprimeAnalyzer.muonSelectors.PFIsoHighPtTightPt20.minPt = 25.
