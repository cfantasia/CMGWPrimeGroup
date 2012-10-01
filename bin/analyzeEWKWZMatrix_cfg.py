from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "EWKWZMatrix.root"
process.WprimeAnalyzer.logFile     = "EWKWZMatrix.dat"
process.WprimeAnalyzer.candEvtFile = "EWKWZMatrix.evt"

process.WprimeAnalyzer.Cuts = EWKWZCuts

process.WprimeAnalyzer.LooseElectronType = "WZLoose"
process.WprimeAnalyzer.TightElectronType = "WZRelaxed"
process.WprimeAnalyzer.ExtraElectronType = cms.untracked.string("WZTight")
process.WprimeAnalyzer.LooseMuonType = "EWKWZLoose"
process.WprimeAnalyzer.TightMuonType = "EWKWZRelaxed"
process.WprimeAnalyzer.ExtraMuonType = cms.untracked.string("EWKWZTight")

process.WprimeAnalyzer.minNLeptons = 2
process.WprimeAnalyzer.met = "patMETsPFType1"

process.WprimeAnalyzer.doMatrix = cms.untracked.bool(True)

    
