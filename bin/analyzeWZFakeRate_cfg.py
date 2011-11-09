from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = 'WZFakeRateElec.root'
process.WprimeAnalyzer.logFile = "WZFakeRateElec.dat"
process.WprimeAnalyzer.candEvtFile = "WZFakeRateElec.evt"
process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_WZDilepton.txt"

######Measure Fake rate of Muons######################
#process.WprimeAnalyzer.outputFile  = 'WZFakeRateMuon.root'
#process.WprimeAnalyzer.logFile = "WZFakeRateMuon.dat"
#process.WprimeAnalyzer.candEvtFile = "WZFakeRateMuon.evt"
#process.WprimeAnalyzer.triggersToUse = SingleElecTriggers
#process.WprimeAnalyzer.Cuts = WZFakeMuonCuts
######Measure Fake rate of Electrons######################
process.WprimeAnalyzer.outputFile  = 'WZFakeRateElec.root'
process.WprimeAnalyzer.logFile = "WZFakeRateElec.dat"
process.WprimeAnalyzer.candEvtFile = "WZFakeRateElec.evt"
process.WprimeAnalyzer.triggersToUse = SingleMuonTriggers
process.WprimeAnalyzer.Cuts = WZFakeElecCuts

process.WprimeAnalyzer.LooseElectronType = "WZRelaxed"
process.WprimeAnalyzer.TightElectronType = "WZTight"
process.WprimeAnalyzer.LooseMuonType = "WZRelaxed"
process.WprimeAnalyzer.TightMuonType = "WZTight"

process.WprimeAnalyzer.minNLeptons = cms.untracked.uint32(2)
