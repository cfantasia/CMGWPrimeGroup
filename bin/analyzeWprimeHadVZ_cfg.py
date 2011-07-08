from UserCode.CMGWPrimeGroup.commonWprime_cfg import *

process.WprimeAnalyzer.outputFile  = cms.string('test.root')## mandatory
process.WprimeAnalyzer.maxEvents   = cms.int32(-1)                      ## optional
process.WprimeAnalyzer.reportAfter = cms.uint32(1)                     ## optional
process.WprimeAnalyzer.useJSON = cms.bool(True)
process.WprimeAnalyzer.logFile     = cms.string("test.log")
process.WprimeAnalyzer.candEvtFile = cms.string("testCandEvtFile.txt")
process.WprimeAnalyzer.sample_cross_sections = cms.string("samples_cross_sections_HadVZ.txt")

## enable analysis in individual channels
process.WprimeAnalyzer.runHadVZAnalysis = cms.bool(True)

## input specific for this analyzer
process.WprimeAnalyzer.muons = cms.string('selectedPatMuons')
process.WprimeAnalyzer.muonAlgo = cms.uint32(0)
process.WprimeAnalyzer.jets = cms.string('selectedPatJets')
process.WprimeAnalyzer.genParticles = cms.InputTag('prunedGenParticles')
process.WprimeAnalyzer.hltEventTag = cms.string('patTriggerEvent')
process.WprimeAnalyzer.pileupTag  = cms.string('addPileupInfo')
#
process.WprimeAnalyzer.maxNumZs = cms.uint32(1)
process.WprimeAnalyzer.minNLeptons =cms.uint32(2)
process.WprimeAnalyzer.minNumJets = cms.uint32(0) # Larger than
process.WprimeAnalyzer.maxNumJets = cms.uint32(9999) # Smaller than 
process.WprimeAnalyzer.minLeadPt = cms.double(20.0)
process.WprimeAnalyzer.maxAngleBetweenJets = cms.double(9999.9)
#
process.WprimeAnalyzer.minZpt = cms.double(0.0) # All units in GeV
process.WprimeAnalyzer.minZmass = cms.double(80.0)
process.WprimeAnalyzer.maxZmass = cms.double(100.0)
process.WprimeAnalyzer.minHadVPt = cms.double(0.0)
process.WprimeAnalyzer.minHadVmass = cms.double(0.0)
process.WprimeAnalyzer.maxHadVmass = cms.double(9999.9)
# +++++++++++++++++++Muon General Cuts
process.WprimeAnalyzer.maxMuonEta = cms.double(2.5)
process.WprimeAnalyzer.minMuonLoosePt = cms.double(10.)
process.WprimeAnalyzer.minMuonTightPt = cms.double(20.)
#VBTF Recommended Cuts
process.WprimeAnalyzer.maxMuonDxy = cms.double(0.2)
process.WprimeAnalyzer.maxMuonNormChi2 = cms.double(10.)
process.WprimeAnalyzer.minMuonNPixHit = cms.int32(0)
process.WprimeAnalyzer.minMuonNTrkHit = cms.int32(10)
process.WprimeAnalyzer.minMuonStations = cms.int32(1)
process.WprimeAnalyzer.minMuonHitsUsed = cms.int32(0)
# +++++++++++++++++++Jet General Cuts
process.WprimeAnalyzer.minJetPt = cms.double(30.0)
process.WprimeAnalyzer.maxJetEta = cms.double(2.4)
process.WprimeAnalyzer.maxJetNHF = cms.double(0.99)
process.WprimeAnalyzer.maxJetNEF = cms.double(0.99)
process.WprimeAnalyzer.minJetnumConst = cms.uint32(1)
process.WprimeAnalyzer.minJetCHF = cms.double(0.0)
process.WprimeAnalyzer.minJetcMult = cms.uint32(0)
process.WprimeAnalyzer.maxJetCEF = cms.double(0.99)

