from PhysicsTools.PatAlgos.tools.coreTools import *

doubleMuonPaths = [
    'HLT_DoubleMu5_v*',
    'HLT_DoubleMu7_v*',
    'HLT_TripleMu5_v*',
    'HLT_Mu13_Mu8_v*', #1e33 unprescaled
    'HLT_Mu17_Mu8_v*' #3e33 unprescaled
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
    
    if mode == 'mu':
        process.hltFilter = cms.Sequence(process.doubleMuonFilter)
    if mode == 'electron':
        process.hltFilter = cms.Sequence(~process.doubleMuonFilter)
