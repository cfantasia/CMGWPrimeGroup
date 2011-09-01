from UserCode.CMGWPrimeGroup.commonWprime_cfg import *

process.WprimeAnalyzer.outputFile  = "HadVZAnalyzer.root"
process.WprimeAnalyzer.logFile     = "HadVZAnalyzer.log"
process.WprimeAnalyzer.candEvtFile = cms.string("HadVZ_CandEvtFile.txt")
process.WprimeAnalyzer.sample_cross_sections = "samples_cross_sections_HadVZ.txt"

process.WprimeAnalyzer.maxEvents   = -1
process.WprimeAnalyzer.reportAfter = 10000
process.WprimeAnalyzer.useJSON = False
process.WprimeAnalyzer.debugme = cms.bool(False)
process.WprimeAnalyzer.preselect = cms.bool(False)

## enable analysis in individual channels
process.WprimeAnalyzer.runHadVZAnalysis = True
process.WprimeAnalyzer.triggersToUse = cms.vstring()

## input specific for this analyzer
process.WprimeAnalyzer.muons = 'selectedPatMuons'
process.WprimeAnalyzer.muonReconstructor = 3
process.WprimeAnalyzer.jets = 'selectedPatJets'
#
process.WprimeAnalyzer.minNLeptons =cms.untracked.uint32(2)
process.WprimeAnalyzer.minNJets = cms.untracked.uint32(1) # Larger than
process.WprimeAnalyzer.maxNJets = cms.uint32(99999) # Smaller than
process.WprimeAnalyzer.maxAngleBetweenJets = cms.double(9999.9)
#
process.WprimeAnalyzer.minZpt = cms.untracked.double(50.0) # All units in GeV
process.WprimeAnalyzer.minZmass = cms.untracked.double(70.0)
process.WprimeAnalyzer.maxZmass = cms.untracked.double(110.0)
process.WprimeAnalyzer.minVPt = cms.untracked.double(50.0)
process.WprimeAnalyzer.minVmass = cms.untracked.double(60.0)
process.WprimeAnalyzer.maxVmass = cms.untracked.double(110.0)
# +++++++++++++++++++Muon General Cuts
process.WprimeAnalyzer.maxMuonEta = cms.double(2.5)
process.WprimeAnalyzer.minMuonLoosePt = cms.double(10.)
process.WprimeAnalyzer.minMuonTightPt = cms.double(35.)
#VBTF Recommended Cuts
process.WprimeAnalyzer.maxMuonDxy = cms.double(0.2)
process.WprimeAnalyzer.maxMuonNormChi2 = cms.double(10.)
process.WprimeAnalyzer.minMuonNPixHit = cms.uint32(1)
process.WprimeAnalyzer.minMuonNTrkHit = cms.uint32(10)
process.WprimeAnalyzer.minMuonStations = cms.uint32(1)
process.WprimeAnalyzer.minMuonHitsUsed = cms.uint32(0)
# +++++++++++++++++++Jet General Cuts
process.WprimeAnalyzer.minJetPt = cms.double(150.)
process.WprimeAnalyzer.maxJetEta = cms.double(2.4)
process.WprimeAnalyzer.maxJetNHF = cms.double(0.99)
process.WprimeAnalyzer.maxJetNEF = cms.double(0.99)
process.WprimeAnalyzer.minJetnumConst = cms.uint32(1)
process.WprimeAnalyzer.minJetCHF = cms.double(0.0)
process.WprimeAnalyzer.minJetcMult = cms.uint32(0)
process.WprimeAnalyzer.maxJetCEF = cms.double(0.99)

process.WprimeAnalyzer.Cuts = cms.vstring(
        "NoCuts",

        "MinNLeptons",
        "MinNJets",

        "ValidV",
        "VMass",
        "VPt",
        
        "ValidZ",
        "ZMass",
        "Zpt",
        "ValidVZ",
              
        #    "HLT",
        "AllCuts")
