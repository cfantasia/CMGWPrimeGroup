from UserCode.CMGWPrimeGroup.commonWprime_cfg import *

process.WprimeAnalyzer.outputFile  = cms.string('Wprime_analysis_ElMET.root')## mandatory
process.WprimeAnalyzer.reportAfter = cms.uint32(15000)                     ## optional
process.WprimeAnalyzer.doRecoilCorrectionForW = cms.bool(False)
process.WprimeAnalyzer.sample_cross_sections = cms.string("samples_cross_sections_ElMET.txt")
process.WprimeAnalyzer.logFile = cms.string("event_counts_ElMET.txt")

## enable analysis in individual channels
process.WprimeAnalyzer.runElMETAnalysis = cms.bool(True)

## input specific for this analyzer
process.WprimeAnalyzer.electrons = cms.InputTag('selectedPatElectrons')
process.WprimeAnalyzer.met   = cms.InputTag('patMETsPF')
process.WprimeAnalyzer.genParticles = cms.InputTag('prunedGenParticles')
#
process.WprimeAnalyzer.oneEleEtCut   = cms.double(25) ## in GeV
process.WprimeAnalyzer.highestEtElectronOnly = cms.bool(False)
process.WprimeAnalyzer.dumpHighEtElectrons   = cms.bool(True)
process.WprimeAnalyzer.dumpHighEtElectronThreshold = cms.double(200)
process.WprimeAnalyzer.barrelCuts = heepBarrelCuts
process.WprimeAnalyzer.endcapCuts = heepEndcapCuts

