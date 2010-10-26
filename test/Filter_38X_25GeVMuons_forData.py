import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

    'rfio:/castor/cern.ch/cms/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/146/511/E0E62F3F-8FC7-DF11-B66F-0019B9F709A4.root'


)
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
# Vertex Filter can be applied if needed
process.primaryVertexFilter = cms.EDFilter("VertexSelector",
   src = cms.InputTag("offlinePrimaryVertices"),
   cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
   filter = cms.bool(True),
)


process.TheHighPtTeVDefaultMuFilter = cms.EDFilter("TrackRefSelector",
                                                   src = cms.InputTag("tevMuons","default"),
                                                   cut = cms.string('pt >25'),
                                                   filter = cms.bool(True),
                                                   minN    = cms.int32(1)
)

process.TheHighPtTeVPickyMuFilter = cms.EDFilter("TrackRefSelector",
                                                 src = cms.InputTag("tevMuons","picky"),
                                                 cut = cms.string('pt >25'),
                                                 filter = cms.bool(True),
                                                 minN    = cms.int32(1)
)

process.TheHighPtTeVFirstHitMuFilter = cms.EDFilter("TrackRefSelector",
                                                    src = cms.InputTag("tevMuons","firstHit"),
                                                    cut = cms.string('pt >25'),
                                                    filter = cms.bool(True),
                                                    minN    = cms.int32(1)
)
process.TheHighPtTeVDYTMuFilter = cms.EDFilter("TrackRefSelector",
                                               src = cms.InputTag("tevMuons","dyt"),
                                               cut = cms.string('pt >25'),
                                               filter = cms.bool(True),
                                               minN    = cms.int32(1)
)

process.TheHighPtTrackerMuFilter = cms.EDFilter("MuonRefSelector",
                                                    src = cms.InputTag("muons"),
                                                    cut = cms.string('isGlobalMuon=1 && innerTrack().pt()>25'),
                                                                                                filter = cms.bool(True),
                                                                                                minN    = cms.int32(1)
                                                                                                )
process.TheHighPtGlbMuFilter = cms.EDFilter("MuonRefSelector",
                                                    src = cms.InputTag("muons"),
                                                    cut = cms.string('isGlobalMuon=1 && globalTrack().pt()>25'),
                                                                                                filter = cms.bool(True),
                                                                                                minN    = cms.int32(1)
                                                                                                )
process.TheHighPtDefaultMuFilter = cms.EDFilter("MuonRefSelector",
                                                    src = cms.InputTag("muons"),
                                                    cut = cms.string('isGlobalMuon=1 && pt>25'),
                                                                                                filter = cms.bool(True),
                                                                                                minN    = cms.int32(1)
                                                                                                )

process.Output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('WprimeMuSkim_25GeVSkim.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p1', 
                                   'p2', 
                                   'p3', 
                                   'p4',
                                   'p5',
                                   'p6',
                                   'p7'
      )
    ),
    basketSize = cms.untracked.int32(4096),
    outputCommands = cms.untracked.vstring('keep *')
)

#process.p1 = cms.Path(process.primaryVertexFilter + process.TheHighPtGlbMuFilter)
#process.p2 = cms.Path(process.primaryVertexFilter + process.TheHighPtTeVDefaultMuFilter)
#process.p3 = cms.Path(process.primaryVertexFilter + process.TheHighPtTeVPickyMuFilter)
#process.p4 = cms.Path(process.primaryVertexFilter + process.TheHighPtTeVFirstHitMuFilter)
#process.p5 = cms.Path(process.primaryVertexFilter + process.TheHighPtTeVDYTMuFilter)
#process.p6 = cms.Path(process.primaryVertexFilter + process.TheHighPtTrackerMuFilter)
process.p1 = cms.Path(process.TheHighPtGlbMuFilter)
process.p2 = cms.Path(process.TheHighPtTeVDefaultMuFilter)
process.p3 = cms.Path(process.TheHighPtTeVPickyMuFilter)
process.p4 = cms.Path(process.TheHighPtTeVFirstHitMuFilter)
process.p5 = cms.Path(process.TheHighPtTeVDYTMuFilter)
process.p6 = cms.Path(process.TheHighPtTrackerMuFilter)
process.p7 = cms.Path(process.TheHighPtDefaultMuFilter)

process.saveout = cms.EndPath(process.Output)

