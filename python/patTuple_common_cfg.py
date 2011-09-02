from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
from PhysicsTools.PatAlgos.tools.metTools import addPfMET
from PhysicsTools.PatAlgos.tools.pfTools import addPFCandidates

def common_config(process, reportEveryNum=100, maxEvents=-1) :
    addPfMET(process, 'PF')

    switchOnTrigger(process)
    process.patTrigger.addL1Algos = cms.bool(True)
    switchOnTrigger(process) # to fix event content
                
    # this is needed so we can correct the pfMET by adjusting the e/mu-pt
    # when switching to one of the dedicated Heep/TeV muon reconstructors
    addPFCandidates(process, 'particleFlow')
    process.selectedPatPFParticles.cut = "abs(pdgId())==11 || abs(pdgId())==13"

    process.load("FWCore.MessageLogger.MessageLogger_cfi")
    process.MessageLogger.cerr.FwkReport.reportEvery = reportEveryNum

    process.maxEvents.input = maxEvents    ##  (e.g. -1 to run on all events)
    process.GlobalTag.globaltag = cms.string('GR_R_42_V18::All')
    #                                         ##
    process.out.outputCommands = [
    # GEN
        'keep *_prunedGenParticles_*_*',
        'keep GenEventInfoProduct_*_*_*',
        'keep GenRunInfoProduct_*_*_*',
    # TRIGGER
        'keep edmTriggerResults_TriggerResults*_*_*',
        'keep *_hltTriggerSummaryAOD_*_*',
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*',
        'keep *_userPat*_*_*',
    # PILEUP
        'keep *_addPileupInfo_*_*',     
    # PF CANDS
        'keep *_selectedPatPFParticles*_*_*'
        ]
    
##  (to suppress the long output at the end of the job)    
    process.options.wantSummary = True        

