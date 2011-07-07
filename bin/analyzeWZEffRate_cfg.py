from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = cms.string('WZEffRate_analysis.root')
process.WprimeAnalyzer.logFile = cms.string("WZEffRate_event_counts.txt")
process.WprimeAnalyzer.candEvtFile = cms.string("WZEffRate_CandEvts.txt")
process.WprimeAnalyzer.sample_cross_sections = cms.string("samples_cross_sections_WZDilepton.txt")

process.WprimeAnalyzer.Cuts = WZEffCuts

process.WprimeAnalyzer.LooseElecCuts = RelaxedElecCuts
process.WprimeAnalyzer.TightElecCuts = TightElecCuts
process.WprimeAnalyzer.LooseMuonCuts = RelaxedMuonCuts
process.WprimeAnalyzer.TightMuonCuts = TightMuonCuts

process.WprimeAnalyzer.minNLeptons = cms.uint32(2)
