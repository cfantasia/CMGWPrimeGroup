from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "EWKWZ.root"
process.WprimeAnalyzer.logFile     = "EWKWZ.dat"
process.WprimeAnalyzer.candEvtFile = "EWKWZ.evt"

process.WprimeAnalyzer.Cuts          = EWKWZCuts

process.WprimeAnalyzer.LooseElectronType = "WZLoose"
process.WprimeAnalyzer.LooseMuonType     = "EWKWZLoose"

process.WprimeAnalyzer.TightElectronType = "WZTight"
process.WprimeAnalyzer.TightMuonType     = "EWKWZTight"

process.WprimeAnalyzer.minNLeptons = 2
process.WprimeAnalyzer.met = "patMETsPFType1"
