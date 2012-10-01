from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = 'WZEffRate.root'
process.WprimeAnalyzer.logFile = "WZEffRate.dat"
process.WprimeAnalyzer.candEvtFile = "WZEffRate.evt"
process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_WZDilepton.txt"

process.WprimeAnalyzer.Cuts = WZEffCuts

process.WprimeAnalyzer.LooseZElectronType = "WZRelaxed80"
process.WprimeAnalyzer.TightZElectronType = "WZTight"

process.WprimeAnalyzer.LooseZMuonType = "EWKWZRelaxed"
process.WprimeAnalyzer.TightZMuonType = "EWKWZTight"

process.WprimeAnalyzer.minNLeptons = 2
process.WprimeAnalyzer.doSystematics = cms.untracked.bool(True)
