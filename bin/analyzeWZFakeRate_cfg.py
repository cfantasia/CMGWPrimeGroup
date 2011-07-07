from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = cms.string('WZFakeRateMuon_analysis.root')## mandatory
process.WprimeAnalyzer.logFile = cms.string("WZFakeRateMuon_event_counts.txt")
process.WprimeAnalyzer.candEvtFile = cms.string("WZFakeRateMuon_CandEvts.txt")
process.WprimeAnalyzer.sample_cross_sections = cms.string("samples_cross_sections_WZDilepton.txt")

process.WprimeAnalyzer.triggersToUse = SingleElecTriggers
#process.WprimeAnalyzer.triggersToUse = SingleMuonTriggers

process.WprimeAnalyzer.Cuts = WZFakeElecCuts
#process.WprimeAnalyzer.Cuts = WZFakeMuonsCuts

process.WprimeAnalyzer.LooseElecCuts = RelaxedElecCuts
process.WprimeAnalyzer.TightElecCuts = TightElecCuts
process.WprimeAnalyzer.LooseMuonCuts = RelaxedMuonCuts
process.WprimeAnalyzer.TightMuonCuts = TightMuonCuts

process.WprimeAnalyzer.minNLeptons = cms.uint32(2)
