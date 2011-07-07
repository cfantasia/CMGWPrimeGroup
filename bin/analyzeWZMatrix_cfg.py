from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = cms.string('WZMatrix_analysis.root')
process.WprimeAnalyzer.logFile = cms.string("WZMatrix_event_counts.txt")
process.WprimeAnalyzer.candEvtFile = cms.string("WZMatrix_CandEvts.txt")

process.WprimeAnalyzer.Cuts = EWKWZCuts

process.WprimeAnalyzer.LooseElecCuts = LooseElecCuts
process.WprimeAnalyzer.TightElecCuts = RelaxedElecCuts
process.WprimeAnalyzer.LooseMuonCuts = LooseMuonCuts
process.WprimeAnalyzer.TightMuonCuts = RelaxedMuonCuts

process.WprimeAnalyzer.minNLeptons = cms.uint32(2)
    
