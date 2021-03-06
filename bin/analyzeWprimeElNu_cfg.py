from UserCode.CMGWPrimeGroup.commonWprime_cfg import *

process.WprimeAnalyzer.outputFile  = cms.string('Wprime_analysis_ElMET.root')## mandatory
process.WprimeAnalyzer.reportAfter = 15000                     ## optional
process.WprimeAnalyzer.useJSON = True
process.WprimeAnalyzer.doRecoilCorrectionForW = cms.bool(False)
process.WprimeAnalyzer.sample_cross_sections = cms.string("samples_cross_sections_ElMET.txt")
process.WprimeAnalyzer.logFile = cms.string("event_counts_ElMET.txt")

## enable analysis in individual channels
process.WprimeAnalyzer.runElMETAnalysis = cms.bool(True)
process.WprimeAnalyzer.maxEvents = cms.int32(-1)

## input specific for this analyzer
process.WprimeAnalyzer.electrons = cms.InputTag('selectedPatElectrons')
process.WprimeAnalyzer.met   = cms.InputTag('patMETsPF')
process.WprimeAnalyzer.genParticles = cms.InputTag('prunedGenParticles')
#
# do not consider electrons below this pt-threshold
process.WprimeAnalyzer.electronPtThreshold   = cms.double(90) ## in GeV
process.WprimeAnalyzer.oneEleEtCut   = cms.double(25) ## in GeV
process.WprimeAnalyzer.highestEtElectronOnly = cms.bool(False)
process.WprimeAnalyzer.dumpHighEtElectrons   = cms.bool(True)
process.WprimeAnalyzer.dumpHighEtElectronThreshold = cms.double(300)
process.WprimeAnalyzer.dumpHighMtElectronThreshold = cms.double(800)

# do or not make Ntuple
process.WprimeAnalyzer.mkTuple = cms.bool(False)

# do or not E/p
process.WprimeAnalyzer.doEoP = cms.bool(False)
process.WprimeAnalyzer.electronEoverPthreshold = cms.double(30)

#HEEP Selection
process.WprimeAnalyzer.barrelCuts = heepBarrelCuts
process.WprimeAnalyzer.endcapCuts = heepEndcapCuts

#WP80 Selection
#process.WprimeAnalyzer.barrelCuts = wp80BarrelCuts
#process.WprimeAnalyzer.endcapCuts = wp80EndcapCuts

#for trigger
# Summer11 MC :
# "HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2","HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1","HLT_Ele45_CaloIdVT_TrkIdT_v2"
# 
process.WprimeAnalyzer.triggerResults = cms.InputTag("TriggerResults::HLT")
process.WprimeAnalyzer.hltPaths = cms.vstring(
    'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1','HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2','HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3',
    'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1','HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2','HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3',
    'HLT_Ele45_CaloIdVT_TrkIdT_v1','HLT_Ele45_CaloIdVT_TrkIdT_v2','HLT_Ele45_CaloIdVT_TrkIdT_v3',
    'HLT_Ele52_CaloIdVT_TrkIdT_v1','HLT_Ele52_CaloIdVT_TrkIdT_v2','HLT_Ele52_CaloIdVT_TrkIdT_v3',
    'HLT_Ele65_CaloIdVT_TrkIdT_v1','HLT_Ele65_CaloIdVT_TrkIdT_v2','HLT_Ele65_CaloIdVT_TrkIdT_v3','HLT_Ele65_CaloIdVT_TrkIdT_v4',
    'HLT_Ele25_WP80_PFMT40_v1',
    'HLT_Ele27_WP80_PFMT50_v1',
    'HLT_Ele80_CaloIdVT_TrkIdT_v2','HLT_Ele80_CaloIdVT_TrkIdT_v3'
    )

process.WprimeAnalyzer.Cuts = cms.vstring("hlt","qual","1el","iso", "met")
