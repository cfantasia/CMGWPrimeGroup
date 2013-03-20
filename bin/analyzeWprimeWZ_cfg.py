from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "WprimeWZ.root"
process.WprimeAnalyzer.logFile     = "WprimeWZ.dat"
process.WprimeAnalyzer.candEvtFile = "WprimeWZ.evt"

process.WprimeAnalyzer.Cuts = WprimeWZCuts

process.WprimeAnalyzer.LooseZElectronType = cms.untracked.string("CiC2012Loose")
process.WprimeAnalyzer.TightZElectronType = cms.untracked.string("CiC2012Loose")
process.WprimeAnalyzer.LooseWElectronType = cms.untracked.string("CiC2012MediumPt20")
process.WprimeAnalyzer.TightWElectronType = cms.untracked.string("CiC2012MediumPt20")

process.WprimeAnalyzer.LooseZMuonType = cms.untracked.string("PFIsoHighPtBoostedZLoose")
process.WprimeAnalyzer.TightZMuonType = cms.untracked.string("PFIsoHighPtBoostedZTight")
process.WprimeAnalyzer.LooseWMuonType = cms.untracked.string("PFIsoHighPtLoosePt20")
process.WprimeAnalyzer.TightWMuonType = cms.untracked.string("PFIsoHighPtTightPt20")



# +++++++++++++++++++Analysis Cuts
process.WprimeAnalyzer.minLt = cms.untracked.double(500.) #W'1000 Cut
process.WprimeAnalyzer.minZpt =  cms.untracked.double(0.)
process.WprimeAnalyzer.minWpt = cms.untracked.double(0.)
