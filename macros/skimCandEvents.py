import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

evtListFile = open("candEvts.txt","r").readlines()
evtList = []
for k in evtListFile:
    evtList.append(k.rstrip())

fileListFile = open("fileList.txt","r").readlines()
fileList = []
for k in fileListFile:
    fileList.append(k.rstrip())

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(fileList),

    eventsToProcess = cms.untracked.VEventRange(evtList)                          
)

process.load("Configuration.EventContent.EventContent_cff")
process.hltPoolOutput = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:candEvts.root')
)
process.HLTOutput = cms.EndPath( process.hltPoolOutput)

