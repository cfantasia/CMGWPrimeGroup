from UserCode.CMGWPrimeGroup.commonWprime_cfg import *

process.WprimeAnalyzer.outputFile  = cms.string('HadVZAnalyzer.root')## mandatory
process.WprimeAnalyzer.maxEvents   = cms.int32(-1)                      ## optional
process.WprimeAnalyzer.reportAfter = cms.uint32(1000)                     ## optional
process.WprimeAnalyzer.useJSON = cms.bool(True)
process.WprimeAnalyzer.logFile     = cms.string("HadVZAnalyzer.log")
process.WprimeAnalyzer.candEvtFile = cms.string("HadVZ_CandEvtFile.txt")
process.WprimeAnalyzer.sample_cross_sections = cms.string("samples_cross_sections_HadVZ.txt")
process.WprimeAnalyzer.debugme = cms.bool(False)
process.WprimeAnalyzer.preselect = cms.bool(False)
## enable analysis in individual channels
process.WprimeAnalyzer.runHadVZAnalysis = cms.bool(True)
process.WprimeAnalyzer.triggersToUse = cms.vstring()
## input specific for this analyzer
process.WprimeAnalyzer.muons = cms.string('selectedPatMuons')
process.WprimeAnalyzer.muonAlgo = cms.uint32(3)
process.WprimeAnalyzer.jets = cms.string('selectedPatJets')
process.WprimeAnalyzer.genParticles = cms.InputTag('prunedGenParticles')
process.WprimeAnalyzer.hltEventTag = cms.string('patTriggerEvent')
process.WprimeAnalyzer.pileupTag  = cms.string('addPileupInfo')
#
process.WprimeAnalyzer.maxNumZs = cms.uint32(2)
process.WprimeAnalyzer.minNLeptons =cms.uint32(1)
process.WprimeAnalyzer.minNJets = cms.uint32(0) # Larger than
process.WprimeAnalyzer.maxNJets = cms.uint32(9999) # Smaller than 
process.WprimeAnalyzer.minLeadPt = cms.double(20.0)
process.WprimeAnalyzer.maxAngleBetweenJets = cms.double(9999.9)
#
process.WprimeAnalyzer.minZpt = cms.double(0.0) # All units in GeV
process.WprimeAnalyzer.minZmass = cms.double(70.0)
process.WprimeAnalyzer.maxZmass = cms.double(110.0)
process.WprimeAnalyzer.minHadVPt = cms.double(0.0)
process.WprimeAnalyzer.minHadVmass = cms.double(60.0)
process.WprimeAnalyzer.maxHadVmass = cms.double(110.0)
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
        
        "MinNTightLeptons",
        "ValidZTight",
        "ZTightMass",
        "ZTightpt",
        "ValidVZTight",
        
        #    "HLT",
        "AllCuts")
