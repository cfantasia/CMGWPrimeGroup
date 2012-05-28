from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = 'WZEffRate.root'
process.WprimeAnalyzer.logFile = "WZEffRate.dat"
process.WprimeAnalyzer.candEvtFile = "WZEffRate.evt"
process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_WZDilepton.txt"

process.WprimeAnalyzer.Cuts = WZEffCuts

process.WprimeAnalyzer.LooseElectronType = "WZRelaxed"
process.WprimeAnalyzer.TightElectronType = "WZTight"
process.WprimeAnalyzer.LooseMuonType = "EWKWZRelaxed"
process.WprimeAnalyzer.TightMuonType = "EWKWZTight"

process.WprimeAnalyzer.minNLeptons = 2
process.WprimeAnalyzer.doSystematics = cms.untracked.bool(True)
