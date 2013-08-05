import sys
sys.path.append("UserCode/CMGWPrimeGroup/bin")
from analyzeWprimeWZ_cfg import *

process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_WZDilepton.txt"

process.WprimeAnalyzer.outputFile  = 'WprimeWZMatrix.root'
process.WprimeAnalyzer.logFile     = "WprimeWZMatrix.dat"
process.WprimeAnalyzer.candEvtFile = "WprimeWZMatrix.evt"

process.WprimeAnalyzer.LooseZElectronType = cms.untracked.string("CiC2012LooseRelaxed")
process.WprimeAnalyzer.TightZElectronType = cms.untracked.string("CiC2012Loose")
process.WprimeAnalyzer.LooseWElectronType = cms.untracked.string("CiC2012MediumRelaxedPt20")
process.WprimeAnalyzer.TightWElectronType = cms.untracked.string("CiC2012MediumPt20")

process.WprimeAnalyzer.LooseZMuonType = cms.untracked.string("PFIsoHighPtBoostedZRelaxed")
process.WprimeAnalyzer.TightZMuonType = cms.untracked.string("PFIsoHighPtBoostedZLoose")
process.WprimeAnalyzer.LooseWMuonType = cms.untracked.string("PFIsoHighPtRelaxedPt20")
process.WprimeAnalyzer.TightWMuonType = cms.untracked.string("PFIsoHighPtTightPt20")

# +++++++++++++++++++Options
process.WprimeAnalyzer.minNLeptons = cms.untracked.uint32(2)
process.WprimeAnalyzer.preselect = False
process.WprimeAnalyzer.doMatrix = cms.untracked.bool(True)
