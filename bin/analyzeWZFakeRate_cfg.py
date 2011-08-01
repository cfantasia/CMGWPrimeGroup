from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = cms.string('WZFakeRateElec_analysis.root')## mandatory
process.WprimeAnalyzer.logFile = cms.string("WZFakeRateElec_event_counts.txt")
process.WprimeAnalyzer.candEvtFile = cms.string("WZFakeRateElec_CandEvts.txt")
process.WprimeAnalyzer.sample_cross_sections = cms.string("samples_cross_sections_WZDilepton.txt")

#process.WprimeAnalyzer.triggersToUse = SingleElecTriggers
#process.WprimeAnalyzer.Cuts = WZFakeMuonCuts
process.WprimeAnalyzer.triggersToUse = SingleMuonTriggers
process.WprimeAnalyzer.Cuts = WZFakeElecCuts

process.WprimeAnalyzer.LooseElectronType = "WZRelaxed"
process.WprimeAnalyzer.TightElectronType = "WZTight"
process.WprimeAnalyzer.LooseMuonType = "WZRelaxed"
process.WprimeAnalyzer.TightMuonType = "WZTight"

process.WprimeAnalyzer.minNLeptons = cms.uint32(2)
