from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_WZDilepton.txt"

process.WprimeAnalyzer.doSystematics = cms.untracked.bool(True)

######Measure Fake rate of Electrons######################
print 'Determining Fake Rate of Electrons'
process.WprimeAnalyzer.outputFile  = 'WZFakeRateElec.root'
process.WprimeAnalyzer.logFile = "WZFakeRateElec.dat"
process.WprimeAnalyzer.candEvtFile = "WZFakeRateElec.evt"
process.WprimeAnalyzer.triggersToUse = SingleMuonTriggers
process.WprimeAnalyzer.Cuts = WZFakeElecCuts

process.WprimeAnalyzer.VLElectronType = cms.untracked.string("WZLoose")
process.WprimeAnalyzer.LooseElectronType = "WZRelaxed"
process.WprimeAnalyzer.TightElectronType = "WZTight"
process.WprimeAnalyzer.VLMuonType = cms.untracked.string("EWKWZLoose")
process.WprimeAnalyzer.LooseMuonType = "EWKWZRelaxed"
process.WprimeAnalyzer.TightMuonType = "EWKWZTight"

process.WprimeAnalyzer.met = "patMETsPFType1"

###Cuts
process.WprimeAnalyzer.minWtransMass = 30.
process.WprimeAnalyzer.minMET = 30.
process.WprimeAnalyzer.maxNVLLeptons = cms.untracked.uint32(2) #This is for Very loose leptons
process.WprimeAnalyzer.minNLeptons = 2 #This is for loose leptonsOB
process.WprimeAnalyzer.minNTightLeptons = cms.untracked.uint32(1)
process.WprimeAnalyzer.maxNJets = cms.untracked.uint32(1)
