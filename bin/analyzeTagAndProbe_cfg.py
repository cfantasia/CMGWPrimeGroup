from UserCode.CMGWPrimeGroup.commonWZ_cfg import *

process.WprimeAnalyzer.runWZAnalysis    = False
process.WprimeAnalyzer.runTagAndProbe   = True

#process.WprimeAnalyzer.maxEvents   = 100
#process.WprimeAnalyzer.debug = True

process.WprimeAnalyzer.outputFile  = "TagAndProbe.root"
process.WprimeAnalyzer.logFile     = "TagAndProbe.dat"
process.WprimeAnalyzer.candEvtFile = "TagAndProbe.evt"

process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_WZDilepton.txt"

process.WprimeAnalyzer.Cuts = WZEffCuts

process.WprimeAnalyzer.CountElectronType = cms.untracked.string("CiC2012Loose")
process.WprimeAnalyzer.TagElectronType   = cms.untracked.string("CiC2012Medium")
process.WprimeAnalyzer.LooseProbeElectronType = cms.untracked.string("CiC2012LooseRelaxed")
process.WprimeAnalyzer.TightProbeElectronType = cms.untracked.string("CiC2012Loose")

process.WprimeAnalyzer.CountMuonType = cms.untracked.string("EWKWZLoose")
process.WprimeAnalyzer.TagMuonType   = cms.untracked.string("EWKWZTight")
process.WprimeAnalyzer.LooseProbeMuonType = cms.untracked.string("EWKWZRelaxed")
process.WprimeAnalyzer.TightProbeMuonType = cms.untracked.string("EWKWZLoose")

# +++++++++++++++++++General Cut values
process.WprimeAnalyzer.maxNumZs = cms.uint32(1)
process.WprimeAnalyzer.minNLeptons = 1

# +++++++++++++++++++Z Cuts
process.WprimeAnalyzer.minZeePt1 =  cms.double(35.)
process.WprimeAnalyzer.minZeePt2 =  cms.double(35.)
process.WprimeAnalyzer.minZmmPt1 =  cms.double(25.)
process.WprimeAnalyzer.minZmmPt2 =  cms.double(10.)

process.WprimeAnalyzer.minZmass =  cms.untracked.double( 71.188)
process.WprimeAnalyzer.maxZmass =  cms.untracked.double(111.188)


# +++++++++++++++++++Options
process.WprimeAnalyzer.adjustMETPhi = False
