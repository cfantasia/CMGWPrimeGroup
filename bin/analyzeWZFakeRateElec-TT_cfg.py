from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_WZDilepton.txt"

process.WprimeAnalyzer.doSystematics = cms.untracked.bool(True)

######Measure Fake rate of Electrons######################
print 'Determining Fake Rate of Electrons'
process.WprimeAnalyzer.outputFile  = 'WZFakeRateElec-TT.root'
process.WprimeAnalyzer.logFile = "WZFakeRateElec-TT.dat"
process.WprimeAnalyzer.candEvtFile = "WZFakeRateElec-TT.evt"
process.WprimeAnalyzer.triggersToUse = SingleMuonTriggers
process.WprimeAnalyzer.Cuts = WZFakeElecCuts

process.WprimeAnalyzer.ExtraElectronType = cms.untracked.string("WZLoose")
process.WprimeAnalyzer.LooseZElectronType = "WZRelaxed95"
process.WprimeAnalyzer.TightZElectronType = "WZLoose"
process.WprimeAnalyzer.LooseWElectronType = "WZLoose"#Doesn't matter
process.WprimeAnalyzer.TightWElectronType = "WZTightPt20"

process.WprimeAnalyzer.ExtraMuonType = cms.untracked.string("EWKWZLoose")
process.WprimeAnalyzer.LooseZMuonType = "EWKWZRelaxed"
process.WprimeAnalyzer.TightZMuonType = "EWKWZLoose"
process.WprimeAnalyzer.LooseWMuonType = "EWKWZLoose"#Doesn't matter
process.WprimeAnalyzer.TightWMuonType = "EWKWZTightPt20"

process.WprimeAnalyzer.met = "patMETsPFType1"

###Cuts
process.WprimeAnalyzer.minWleptPt = 20.
process.WprimeAnalyzer.minWtransMass = 30.
process.WprimeAnalyzer.minMET = 30.
process.WprimeAnalyzer.maxNVLLeptons = cms.untracked.uint32(2) #This is for Very loose leptons
process.WprimeAnalyzer.minNLeptons = 2 #This is for loose leptonsOB
process.WprimeAnalyzer.minNTightLeptons = cms.untracked.uint32(1)
process.WprimeAnalyzer.maxNJets = cms.untracked.uint32(1)
