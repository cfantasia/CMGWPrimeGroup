from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = cms.string('WprimeWZMatrix_analysis.root')
process.WprimeAnalyzer.logFile     = cms.string("WprimeWZMatrix_event_counts.txt")
process.WprimeAnalyzer.candEvtFile = cms.string("WprimeWZMatrix_CandEvts.txt")

process.WprimeAnalyzer.Cuts = WprimeWZCuts

process.WprimeAnalyzer.LooseElectronType = "WZLoose"
process.WprimeAnalyzer.TightElectronType = "WZRelaxed"
process.WprimeAnalyzer.LooseMuonType = "WZLoose"
process.WprimeAnalyzer.TightMuonType = "WZRelaxed"

process.WprimeAnalyzer.minNLeptons = cms.uint32(2)
    
