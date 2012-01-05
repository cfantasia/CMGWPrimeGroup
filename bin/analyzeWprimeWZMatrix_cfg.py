from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = 'WprimeWZMatrix.root'
process.WprimeAnalyzer.logFile     = "WprimeWZMatrix.dat"
process.WprimeAnalyzer.candEvtFile = "WprimeWZMatrix.evt"

process.WprimeAnalyzer.Cuts = WprimeWZCuts

process.WprimeAnalyzer.minMET = 30.

process.WprimeAnalyzer.LooseElectronType = "WZLoose"
process.WprimeAnalyzer.TightElectronType = "WZRelaxed"
process.WprimeAnalyzer.LooseMuonType = "WZLoose"
process.WprimeAnalyzer.TightMuonType = "WZRelaxed"
process.WprimeAnalyzer.minNLeptons = cms.untracked.uint32(2)

# +++++++++++++++++++Analysis Cuts
process.WprimeAnalyzer.minHt = cms.untracked.double(270.) #W'600 Cut
process.WprimeAnalyzer.minZpt =  cms.untracked.double(0.)
process.WprimeAnalyzer.minWpt = cms.untracked.double(0.)
    
