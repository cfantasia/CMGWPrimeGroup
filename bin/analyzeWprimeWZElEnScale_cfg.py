from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "WprimeWZElEnScale.root"
process.WprimeAnalyzer.logFile     = "WprimeWZElEnScale.dat"
process.WprimeAnalyzer.candEvtFile = "WprimeWZElEnScale.evt"

process.WprimeAnalyzer.Cuts = WprimeWZCuts

process.WprimeAnalyzer.LooseZElectronType = cms.untracked.string("CiC2012Loose")
process.WprimeAnalyzer.TightZElectronType = cms.untracked.string("CiC2012Loose")
process.WprimeAnalyzer.LooseWElectronType = cms.untracked.string("CiC2012MediumPt20")
process.WprimeAnalyzer.TightWElectronType = cms.untracked.string("CiC2012MediumPt20")

process.WprimeAnalyzer.LooseZMuonType = cms.untracked.string("EWKWZLoose")
process.WprimeAnalyzer.TightZMuonType = cms.untracked.string("EWKWZTight")
process.WprimeAnalyzer.LooseWMuonType = cms.untracked.string("EWKWZLoosePt20")
process.WprimeAnalyzer.TightWMuonType = cms.untracked.string("EWKWZTightPt20")

# +++++++++++++++++++Analysis Cuts
process.WprimeAnalyzer.minLt = cms.untracked.double(290.) #W'600 Cut
process.WprimeAnalyzer.minZpt =  cms.untracked.double(0.)
process.WprimeAnalyzer.minWpt = cms.untracked.double(0.)

# +++++++++++++++++++Options
process.WprimeAnalyzer.elScaleFactor = 1.02
