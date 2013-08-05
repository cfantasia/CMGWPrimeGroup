from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "ZJetsFakeRate-ZIso.root"
process.WprimeAnalyzer.logFile     = "ZJetsFakeRate-ZIso.dat"
process.WprimeAnalyzer.candEvtFile = "ZJetsFakeRate-ZIso.evt"
process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_WZDilepton.txt"

process.WprimeAnalyzer.Cuts = EWKWZCuts

process.WprimeAnalyzer.ExtraElectronType = cms.untracked.string("CiC2012Veto")
process.WprimeAnalyzer.ExtraMuonType = cms.untracked.string("PFIsoHighPtLoose")

process.WprimeAnalyzer.LooseZElectronType = cms.untracked.string("CiC2012Loose")
process.WprimeAnalyzer.TightZElectronType = cms.untracked.string("CiC2012Loose")
process.WprimeAnalyzer.LooseWElectronType = cms.untracked.string("CiC2012LooseRelaxed")
process.WprimeAnalyzer.TightWElectronType = cms.untracked.string("CiC2012Loose")

process.WprimeAnalyzer.LooseZMuonType = cms.untracked.string("PFIsoHighPtBoostedZTight")
process.WprimeAnalyzer.TightZMuonType = cms.untracked.string("PFIsoHighPtBoostedZTight")
process.WprimeAnalyzer.LooseWMuonType = cms.untracked.string("PFIsoHighPtBoostedZRelaxed")
process.WprimeAnalyzer.TightWMuonType = cms.untracked.string("PFIsoHighPtBoostedZTight")

process.WprimeAnalyzer.minNLeptons = 2

process.WprimeAnalyzer.minZmass =   71.188
process.WprimeAnalyzer.maxZmass =  111.188

process.WprimeAnalyzer.doMatrix = cms.untracked.bool(True)
process.WprimeAnalyzer.preselect = False
process.WprimeAnalyzer.useJSON = False
