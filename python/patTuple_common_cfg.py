from PhysicsTools.PatAlgos.tools.metTools import addPfMET

def common_config(process, reportEveryNum=100, maxEvents=-1) :
    addPfMET(process, 'PF')
    process.load("FWCore.MessageLogger.MessageLogger_cfi")
    process.MessageLogger.cerr.FwkReport.reportEvery = reportEveryNum

    process.maxEvents.input = maxEvents    ##  (e.g. -1 to run on all events)
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
        'keep *_addPileupInfo_*_*',     ]

##  (to suppress the long output at the end of the job)    
    process.options.wantSummary = True        

