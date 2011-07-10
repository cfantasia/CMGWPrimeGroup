from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = cms.string('WZFakeRateMuon_analysis.root')## mandatory
process.WprimeAnalyzer.logFile = cms.string("WZFakeRateMuon_event_counts.txt")
process.WprimeAnalyzer.candEvtFile = cms.string("WZFakeRateMuon_CandEvts.txt")
process.WprimeAnalyzer.sample_cross_sections = cms.string("samples_cross_sections_WZDilepton.txt")

process.WprimeAnalyzer.triggersToUse = SingleElecTriggers
#process.WprimeAnalyzer.triggersToUse = SingleMuonTriggers

process.WprimeAnalyzer.Cuts = WZFakeElecCuts
#process.WprimeAnalyzer.Cuts = WZFakeMuonsCuts

process.WprimeAnalyzer.LooseElecCuts = "WZRelaxed"
process.WprimeAnalyzer.TightElecCuts = "WZTight"
process.WprimeAnalyzer.LooseMuonCuts = "WZRelaxed"
process.WprimeAnalyzer.TightMuonCuts = "WZTight"

process.WprimeAnalyzer.minNLeptons = cms.uint32(2)
