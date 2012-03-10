from PhysicsTools.PatAlgos.tools.coreTools import *

singleMuonPaths = [
    'HLT_Mu15_v*',
    'HLT_Mu17_v*',
    'HLT_Mu20_v*',
    'HLT_Mu24_v*',
    'HLT_Mu30_v*',
    'HLT_Mu40_v*',
    'HLT_Mu40_eta2p1_v*',
    'HLT_Mu50_eta2p1_v*',
    ]

doubleMuonPaths = [
    'HLT_DoubleMu7_v*',
    'HLT_Mu13_Mu8_v*', #1e33 unprescaled
    'HLT_Mu17_Mu8_v*' #3e33 unprescaled
    'HLT_Mu17_TkMu8_v*', #2011B
    ]

singleElectronPaths = [
    ]

doubleElectronPaths = [
    'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*',
    'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*',#early version, MC
    'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*'
    ]


def addHLTFilter(process, hltProcess='HLT', mode='mc'):
    "Add HLT filter used to keep datasets orthogonal."
    if mode == 'mc' or mode == 'allmueg':
        ## Give hltFilter a dummy producer (so it always passes).
        ## This will be removed anyway in configureFilters.
        process.hltFilter = cms.EDProducer("EventCountProducer")
        return
    from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
    hltHighLevel.andOr = True # True means 'OR'; False means 'AND'
    hltHighLevel.throw = False # Don't die on unknown path names
    hltHighLevel.TriggerResultsTag = cms.InputTag("TriggerResults", "", hltProcess)
    process.doubleMuonFilter = hltHighLevel.clone()
    process.doubleMuonFilter.HLTPaths = doubleMuonPaths
    process.singleMuonFilter = hltHighLevel.clone()
    process.singleMuonFilter.HLTPaths = singleMuonPaths
    process.doubleElectronFilter = hltHighLevel.clone()
    process.doubleElectronFilter.HLTPaths = doubleElectronPaths
    process.singleElectronFilter = hltHighLevel.clone()
    process.singleElectronFilter.HLTPaths = singleElectronPaths

    if mode == 'doublemu':
        process.hltFilter = cms.Sequence(process.doubleMuonFilter)
    if mode == 'doubleelectron':
        process.hltFilter = cms.Sequence(process.doubleElectronFilter *
                                         ~process.doubleMuonFilter)
    if mode == 'singlemu':
        process.hltFilter = cms.Sequence(process.singleMuonFilter)
        
    if mode == 'singleelectron':
        process.hltFilter = cms.Sequence(process.singleElectronFilter *
                                         ~process.singleMuonFilter)
        
        
