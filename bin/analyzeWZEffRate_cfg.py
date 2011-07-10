from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = cms.string('WZEffRate_analysis.root')
process.WprimeAnalyzer.logFile = cms.string("WZEffRate_event_counts.txt")
process.WprimeAnalyzer.candEvtFile = cms.string("WZEffRate_CandEvts.txt")
process.WprimeAnalyzer.sample_cross_sections = cms.string("samples_cross_sections_WZDilepton.txt")

process.WprimeAnalyzer.Cuts = WZEffCuts

process.WprimeAnalyzer.LooseElecCuts = "WZRelaxed"
process.WprimeAnalyzer.TightElecCuts = "WZTight"
process.WprimeAnalyzer.LooseMuonCuts = "WZRelaxed"
process.WprimeAnalyzer.TightMuonCuts = "WZTight"

process.WprimeAnalyzer.minNLeptons = cms.uint32(2)
