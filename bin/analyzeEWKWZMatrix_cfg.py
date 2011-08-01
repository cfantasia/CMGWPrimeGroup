from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = cms.string('EWKWZMatrix_analysis.root')
process.WprimeAnalyzer.logFile = cms.string("EWKWZMatrix_event_counts.txt")
process.WprimeAnalyzer.candEvtFile = cms.string("EWKWZMatrix_CandEvts.txt")

process.WprimeAnalyzer.Cuts = EWKWZCuts

process.WprimeAnalyzer.LooseElectronType = "WZLoose"
process.WprimeAnalyzer.TightElectronType = "WZRelaxed"
process.WprimeAnalyzer.LooseMuonType = "WZLoose"
process.WprimeAnalyzer.TightMuonType = "WZRelaxed"

process.WprimeAnalyzer.minNLeptons = cms.uint32(2)
    
