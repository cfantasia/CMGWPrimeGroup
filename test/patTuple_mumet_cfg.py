## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *

## ------------------------------------------------------
#  NOTE: you can use a bunch of core tools of PAT to
#  taylor your PAT configuration; for a few examples
#  uncomment the lines below
## ------------------------------------------------------
from PhysicsTools.PatAlgos.tools.coreTools import *

## remove MC matching from the default sequence
removeMCMatching(process, ['All'])

from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')

from PhysicsTools.PatAlgos.tools.pfTools import *
addPFCandidates(process, 'particleFlow')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

### Prune the GEN particle collection ###
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.prunedGenParticles = cms.EDProducer("GenParticlePruner",
#                                            src = cms.InputTag("genParticles"),
#                                            select = cms.vstring(
#    "drop  *",
    #keeps all particles from the hard matrix element
#    "keep status = 3",
    #keeps all stable muons (13) and electrons (11) + neutrinos (14, 12)
#    # + W (24) + W' (34)  and their (direct) mothers.
#    "+keep (abs(pdgId) = 11 | abs(pdgId) = 13 | abs(pdgId) = 12 | abs(pdgId) = 14 | abs(pdgId) = 24 | abs(pdgId) = 34) & status = 1"
#    )
#)

process.selectedPatMuons.cut = "pt > 25. & abs(eta) < 2.5"
process.selectedPatPFParticles.cut = "particleId == 3"

process.highPtMuonFilter = cms.EDFilter("CandViewCountFilter",
                                src = cms.InputTag("selectedPatMuons"),
                                minNumber = cms.uint32(1)
                                )

## let it run
process.p = cms.Path(
    process.patDefaultSequence *
#    process.prunedGenParticles * 
    process.highPtMuonFilter
)

## ------------------------------------------------------
#  In addition you usually want to change the following
#  parameters:
## ------------------------------------------------------
#
#   process.GlobalTag.globaltag =  ...    ##  (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)
#                                         ##
process.source.fileNames = [          ##
    '/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/161/312/F8AEC745-DF57-E011-8D23-001D09F290BF.root'
#    '/store/relval/CMSSW_3_8_6/RelValTTbar/GEN-SIM-RECO/START38_V13-v1/0068/98EA8C65-25E8-DF11-B0A0-0018F3D095F8.root'
    ]                                     ##  (e.g. 'file:AOD.root')
#                                         ##
process.maxEvents.input = -1        ##  (e.g. -1 to run on all events)
#                                         ##
process.out.outputCommands = [
    # RECO
    'keep *_selectedPatMuons*_*_*',
    'keep *_patMETs*_*_*',
    'keep *_selectedPatPFParticles*_*_*',
    # GEN
    'keep *_prunedGenParticles_*_*',
#    'keep recoGenParticles_genParticles*_*_*',
    'keep GenEventInfoProduct_*_*_*',
    'keep GenRunInfoProduct_*_*_*',
    # TRIGGER
    'keep edmTriggerResults_TriggerResults*_*_*',
    'keep *_hltTriggerSummaryAOD_*_*'
     ]

process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )

#                                         ##
#   process.out.fileName = ...            ##  (e.g. 'myTuple.root')
#                                         ##
process.options.wantSummary = True        ##  (to suppress the long output at the end of the job)    

