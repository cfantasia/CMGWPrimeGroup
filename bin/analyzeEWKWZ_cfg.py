from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "EWKWZ.root"
process.WprimeAnalyzer.logFile     = "EWKWZ.dat"
process.WprimeAnalyzer.candEvtFile = "EWKWZ.evt"

process.WprimeAnalyzer.Cuts          = EWKWZCuts

process.WprimeAnalyzer.LooseZElectronType = "EWKWZTight"
process.WprimeAnalyzer.LooseZMuonType     = "EWKWZTight"
process.WprimeAnalyzer.TightZElectronType = "EWKWZTight"
process.WprimeAnalyzer.TightZMuonType     = "EWKWZTight"

process.WprimeAnalyzer.LooseWElectronType = "EWKWZTight"
process.WprimeAnalyzer.LooseWMuonType     = "EWKWZTight"
process.WprimeAnalyzer.TightWElectronType = "EWKWZTightPt20"
process.WprimeAnalyzer.TightWMuonType     = "EWKWZTightPt20"

process.WprimeAnalyzer.ExtraElectronType = cms.untracked.string("EWKWZTightPt20")
process.WprimeAnalyzer.ExtraMuonType     = cms.untracked.string("EWKWZTightPt20")

process.WprimeAnalyzer.minNLeptons = 2
process.WprimeAnalyzer.met = "patMETsPFType1"

process.WprimeAnalyzer.minZmass =   71.188
process.WprimeAnalyzer.maxZmass =  111.188
