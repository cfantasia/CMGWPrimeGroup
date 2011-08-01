from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = cms.string('WZEffRate_analysis.root')
process.WprimeAnalyzer.logFile = cms.string("WZEffRate_event_counts.txt")
process.WprimeAnalyzer.candEvtFile = cms.string("WZEffRate_CandEvts.txt")
process.WprimeAnalyzer.sample_cross_sections = cms.string("samples_cross_sections_WZDilepton.txt")

process.WprimeAnalyzer.Cuts = WZEffCuts

process.WprimeAnalyzer.LooseElectronType = "WZRelaxed"
process.WprimeAnalyzer.TightElectronType = "WZTight"
process.WprimeAnalyzer.LooseMuonType = "WZRelaxed"
process.WprimeAnalyzer.TightMuonType = "WZTight"

process.WprimeAnalyzer.minNLeptons = cms.uint32(2)
