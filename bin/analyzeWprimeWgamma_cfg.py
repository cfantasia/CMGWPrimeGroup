from UserCode.CMGWPrimeGroup.commonWprime_cfg import *

process.WprimeAnalyzer.outputFile  = cms.string('Wgamma_analysis.root')## mandatory
process.WprimeAnalyzer.sample_cross_sections = cms.string("samples_cross_sections_MuMET.txt")
process.WprimeAnalyzer.logFile = cms.string("Wprime_event_counts.txt")

## enable analysis in individual channels
process.WprimeAnalyzer.runWgammaAnalysis = cms.bool(True)

## input specific for this analyzer
process.WprimeAnalyzer.muons = cms.InputTag('selectedPatMuons')
process.WprimeAnalyzer.mets   = cms.InputTag('patMETsPF')
process.WprimeAnalyzer.particleFlow = cms.InputTag('selectedPatPFParticles')
process.WprimeAnalyzer.photons = cms.InputTag('selectedPatPhotons')
process.WprimeAnalyzer.genParticles = cms.InputTag('prunedGenParticles')
#
process.WprimeAnalyzer.muonReconstructor = cms.int32(3) ## see TeVMuon_tracking.h
# do not consider muons below this pt-threshold
process.WprimeAnalyzer.muonPtThreshold   = cms.double(10) ## in GeV
process.WprimeAnalyzer.oneMuPtTrackCut   = cms.double(25) ## in GeV
process.WprimeAnalyzer.chi2Cut           = cms.double(10)
process.WprimeAnalyzer.muonEtaCut        = cms.double(2.1)
process.WprimeAnalyzer.combRelCut        = cms.double(0.15)
process.WprimeAnalyzer.highestPtMuonOnly = cms.bool(True)
process.WprimeAnalyzer.highestPtPhotonOnly = cms.bool(False)
process.WprimeAnalyzer.dumpHighPtMuons   = cms.bool(True)
process.WprimeAnalyzer.dumpHighPtMuonThreshold = cms.double(200)
process.WprimeAnalyzer.dumpHighPtPhotons = cms.bool(True)
process.WprimeAnalyzer.dumpHighPtPhotonThreshold = cms.double(100)
process.WprimeAnalyzer.BarrelJurrasicECALIsoConst = cms.double(4.2)
process.WprimeAnalyzer.BarrelJurrasicECALIsoSlope = cms.double(0.006)
process.WprimeAnalyzer.BarrelTowerHCALIsoConst = cms.double(2.2)
process.WprimeAnalyzer.BarrelTowerHCALIsoSlope = cms.double(0.0025)
process.WprimeAnalyzer.BarrelMaxHadronicOverEm = cms.double(0.05)
process.WprimeAnalyzer.BarrelHollowConeTrkIsoConst = cms.double(3.5)
process.WprimeAnalyzer.BarrelHollowConeTrkIsoSlope = cms.double(0.001)
process.WprimeAnalyzer.BarrelMaxSigmaIetaIeta = cms.double(99999.9)
process.WprimeAnalyzer.EndcapJurrasicECALIsoConst = cms.double(4.2)
process.WprimeAnalyzer.EndcapJurrasicECALIsoSlope = cms.double(0.006)
process.WprimeAnalyzer.EndcapTowerHCALIsoConst = cms.double(2.2)
process.WprimeAnalyzer.EndcapTowerHCALIsoSlope = cms.double(0.0025)
process.WprimeAnalyzer.EndcapMaxHadronicOverEm = cms.double(0.05)
process.WprimeAnalyzer.EndcapHollowConeTrkIsoConst = cms.double(3.5)
process.WprimeAnalyzer.EndcapHollowConeTrkIsoSlope = cms.double(0.001)
process.WprimeAnalyzer.EndcapMaxSigmaIetaIeta = cms.double(99999.9)
process.WprimeAnalyzer.ApplyTrackVeto = cms.bool(True)
process.WprimeAnalyzer.minPt = cms.double(10)
process.WprimeAnalyzer.maxEta = cms.double(5)
