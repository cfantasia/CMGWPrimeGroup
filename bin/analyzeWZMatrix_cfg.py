from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = cms.string('WZMatrix_analysis.root')
process.WprimeAnalyzer.logFile = cms.string("WZMatrix_event_counts.txt")
process.WprimeAnalyzer.candEvtFile = cms.string("WZMatrix_CandEvts.txt")

process.WprimeAnalyzer.Cuts = EWKWZCuts

process.WprimeAnalyzer.LooseElecType = "WZLoose"
process.WprimeAnalyzer.TightElecType = "WZRelaxed"
process.WprimeAnalyzer.LooseMuonType = "WZLoose"
process.WprimeAnalyzer.TightMuonType = "WZRelaxed"

process.WprimeAnalyzer.minNLeptons = cms.uint32(2)
    
