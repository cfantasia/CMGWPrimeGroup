from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "EWKWZMatrix.root"
process.WprimeAnalyzer.logFile     = "EWKWZMatrix.dat"
process.WprimeAnalyzer.candEvtFile = "EWKWZMatrix.evt"

process.WprimeAnalyzer.Cuts = EWKWZCuts

process.WprimeAnalyzer.LooseElectronType = "WZLoose"
process.WprimeAnalyzer.TightElectronType = "WZRelaxed"
process.WprimeAnalyzer.LooseMuonType = "WZLoose"
process.WprimeAnalyzer.TightMuonType = "EWKWZRelaxed"

process.WprimeAnalyzer.minNLeptons = 2
process.WprimeAnalyzer.met = "patMETsPFType1"
    
