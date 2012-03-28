from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

#process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_WZ.txt"
process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_WZDilepton.txt"

process.WprimeAnalyzer.doSystematics = cms.untracked.bool(True)

######Measure Fake rate of Muons######################
if False:
#if True:
    print 'Determining Fake Rate of Electrons'
    process.WprimeAnalyzer.outputFile  = 'WZFakeRateElec.root'
    process.WprimeAnalyzer.logFile = "WZFakeRateElec.dat"
    process.WprimeAnalyzer.candEvtFile = "WZFakeRateElec.evt"
    process.WprimeAnalyzer.triggersToUse = SingleMuonTriggers
    process.WprimeAnalyzer.Cuts = WZFakeElecCuts
else :
    print 'Determining Fake Rate of Muons'
    process.WprimeAnalyzer.outputFile  = 'WZFakeRateMuon.root'
    process.WprimeAnalyzer.logFile = "WZFakeRateMuon.dat"
    process.WprimeAnalyzer.candEvtFile = "WZFakeRateMuon.evt"
    process.WprimeAnalyzer.triggersToUse = SingleElecTriggers
    process.WprimeAnalyzer.Cuts = WZFakeMuonCuts

process.WprimeAnalyzer.LooseElectronType = "WZRelaxed"
process.WprimeAnalyzer.TightElectronType = "WZTight"
process.WprimeAnalyzer.LooseMuonType = "WZRelaxed"
process.WprimeAnalyzer.TightMuonType = "WZTight"

###Cuts
process.WprimeAnalyzer.minWtransMass = 30.
process.WprimeAnalyzer.minMET = 30.
process.WprimeAnalyzer.minMET = 30.
process.WprimeAnalyzer.minNLeptons = 2
process.WprimeAnalyzer.minNTightLeptons = cms.untracked.uint32(1)
