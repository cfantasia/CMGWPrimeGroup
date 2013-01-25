from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "WprimeWZPUSys.root"
process.WprimeAnalyzer.logFile     = "WprimeWZPUSys.dat"
process.WprimeAnalyzer.candEvtFile = "WprimeWZPUSys.evt"

process.WprimeAnalyzer.Cuts = WprimeWZCuts

if False: #2012 pfIso
    process.WprimeAnalyzer.LooseZElectronType = cms.untracked.string("exotica")
    process.WprimeAnalyzer.TightZElectronType = cms.untracked.string("exotica")
    process.WprimeAnalyzer.LooseWElectronType = cms.untracked.string("exoticaPt20")
    process.WprimeAnalyzer.TightWElectronType = cms.untracked.string("exoticaPt20")

    process.WprimeAnalyzer.LooseZMuonType = cms.untracked.string("EWKWZLoose")
    process.WprimeAnalyzer.TightZMuonType = cms.untracked.string("EWKWZTight")
    process.WprimeAnalyzer.LooseWMuonType = cms.untracked.string("EWKWZLoose")
    process.WprimeAnalyzer.TightWMuonType = cms.untracked.string("EWKWZTight")

elif False: #2011
    process.WprimeAnalyzer.LooseZElectronType = cms.untracked.string("WZLoose")
    process.WprimeAnalyzer.TightZElectronType = cms.untracked.string("WZTight")
    process.WprimeAnalyzer.LooseWElectronType = cms.untracked.string("WZLoose")
    process.WprimeAnalyzer.TightWElectronType = cms.untracked.string("WZTight")

    process.WprimeAnalyzer.LooseZMuonType = cms.untracked.string("WZLoose")
    process.WprimeAnalyzer.TightZMuonType = cms.untracked.string("WZTight")
    process.WprimeAnalyzer.LooseWMuonType = cms.untracked.string("WZLoose")
    process.WprimeAnalyzer.TightWMuonType = cms.untracked.string("WZTight")

else:
    process.WprimeAnalyzer.LooseZElectronType = cms.untracked.string("CiC2012Loose")
    process.WprimeAnalyzer.TightZElectronType = cms.untracked.string("CiC2012Loose")
    process.WprimeAnalyzer.LooseWElectronType = cms.untracked.string("CiC2012MediumPt20")
    process.WprimeAnalyzer.TightWElectronType = cms.untracked.string("CiC2012MediumPt20")

    process.WprimeAnalyzer.LooseZMuonType = cms.untracked.string("EWKWZLoose")
    process.WprimeAnalyzer.TightZMuonType = cms.untracked.string("EWKWZTight")
    process.WprimeAnalyzer.LooseWMuonType = cms.untracked.string("EWKWZLoosePt20")
    process.WprimeAnalyzer.TightWMuonType = cms.untracked.string("EWKWZTightPt20")



# +++++++++++++++++++Analysis Cuts
process.WprimeAnalyzer.minLt = cms.untracked.double(290.) #W'600 Cut
process.WprimeAnalyzer.minZpt =  cms.untracked.double(0.)
process.WprimeAnalyzer.minWpt = cms.untracked.double(0.)

# +++++++++++++++++++Options
process.WprimeAnalyzer.DataPUDistFile = 'UserCode/CMGWPrimeGroup/pileup/DataPileupHistogram-73p5mb.root'
