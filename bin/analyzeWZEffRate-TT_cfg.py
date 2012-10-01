from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = 'WZEffRate-TT.root'
process.WprimeAnalyzer.logFile = "WZEffRate-TT.dat"
process.WprimeAnalyzer.candEvtFile = "WZEffRate-TT.evt"
process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_WZDilepton.txt"

process.WprimeAnalyzer.Cuts = WZEffCuts

process.WprimeAnalyzer.LooseZElectronType = "WZRelaxed95"
process.WprimeAnalyzer.TightZElectronType = "WZLoose"

process.WprimeAnalyzer.LooseZMuonType = "EWKWZRelaxed"
process.WprimeAnalyzer.TightZMuonType = "EWKWZLoose"

process.WprimeAnalyzer.minNLeptons = 2
process.WprimeAnalyzer.doSystematics = cms.untracked.bool(True)
