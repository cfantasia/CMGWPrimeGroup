from UserCode.CMGWPrimeGroup.commonWprime_cfg import *

process.WprimeAnalyzer.outputFile  = cms.string('Wprime_analysis_MuMET.root')## mandatory
process.WprimeAnalyzer.reportAfter = cms.uint32(15000)                     ## optional
process.WprimeAnalyzer.useJSON = True
process.WprimeAnalyzer.doRecoilCorrectionForW = cms.bool(False)
process.WprimeAnalyzer.sample_cross_sections = cms.string("samples_cross_sections_MuMET.txt")
process.WprimeAnalyzer.logFile = cms.string("event_counts_MuMET.txt")
## enable analysis in individual channels
process.WprimeAnalyzer.runMuMETAnalysis = cms.bool(True)
process.WprimeAnalyzer.maxEvents = cms.int32(-1)

## input specific for this analyzer
process.WprimeAnalyzer.muons = cms.InputTag('selectedPatMuons')
process.WprimeAnalyzer.met   = cms.InputTag('patMETsPF')
process.WprimeAnalyzer.particleFlow = cms.InputTag('selectedPatPFParticles')
process.WprimeAnalyzer.genParticles = cms.InputTag('prunedGenParticles')
#
process.WprimeAnalyzer.muonReconstructor = cms.int32(3) ## see TeVMuon_tracking.h
# do not consider muons below this pt-threshold
process.WprimeAnalyzer.muonPtThreshold   = cms.double(90) ## in GeV
process.WprimeAnalyzer.oneMuPtTrackCut   = cms.double(25) ## in GeV
process.WprimeAnalyzer.chi2Cut           = cms.double(10)
process.WprimeAnalyzer.muonEtaCut        = cms.double(2.1)
process.WprimeAnalyzer.relIsoCut        = cms.double(0.10)
process.WprimeAnalyzer.highestPtMuonOnly = cms.bool(False)
process.WprimeAnalyzer.dumpHighPtMuons   = cms.bool(True)
process.WprimeAnalyzer.dumpHighPtMuonThreshold = cms.double(300)
process.WprimeAnalyzer.dumpHighMtMuonThreshold = cms.double(800)
process.WprimeAnalyzer.useAdjustedMET = cms.bool(True)

process.inputs.lumisToProcess.extend(MuonPhysJSONList)
