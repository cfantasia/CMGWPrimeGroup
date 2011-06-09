from PhysicsTools.PatAlgos.patTemplate_cfg import *
#process.load('CommonTools.ParticleFlow.PF2PAT_cff')
#from PhysicsTools.PFCandProducer.PF2PAT_cff import *
process.load('PhysicsTools.PFCandProducer.PF2PAT_cff')
from PhysicsTools.PatAlgos.tools.coreTools import *
from UserCode.CMGWPrimeGroup.patTuple_jet_cfg import *

## remove MC matching from the default sequence when running on data
removeMCMatching(process, ['All'])

# 2nd argument: message-logger frequency
# 3rd argument: # of events to process
jet_config(process, 1000, 5000)

# keep all events with jet-pt above 30 GeV, |eta| < 2.4
process.selectedPatJets.cut = "pt > 30. & abs(eta) < 2.4"
process.lowPtJetFilter = cms.EDFilter("CandViewCountFilter",
                                src = cms.InputTag("selectedPatJets"),
                                minNumber = cms.uint32(1)
                                )

# Setup so that my patJets are PFJets and not calojets
process.patJets.jetSource='pfJets'
# turn off all the corrections and extra factors for now
process.patJets.addJetCorrFactors = False
process.patJets.addBTagInfo = False
process.patJets.addDiscriminators = False
process.patJets.addAssociatedTracks = False
process.patJets.addJetCharge = False
process.patJets.addJetID = False
process.patJets.addGenPartonMatch = False
process.patJets.embedGenPartonMatch = False
process.patJets.genPartonMatch = ''
process.patJets.addGenJetMatch = False
process.patJets.embedGenJetMatch = False
process.patJets.getJetMCFlavour = False
process.patJets.JetPartonMapSource = '' 

# here I define which sequence I want to be made from PF2PAT - will made the process up to jets only, not care about it does later
process.PF2PATmod = cms.Sequence(
    process.pfNoPileUpSequence +
    process.pfAllNeutralHadrons+
    process.pfAllChargedHadrons+
    process.pfAllPhotons+
    process.pfMuonSequence +
    process.pfNoMuon +
    process.pfElectronSequence +
    process.pfNoElectron +
    process.pfJetSequence 
    )

## let it run
process.p = cms.Path(
    process.PF2PATmod *
    (process.patJets +
     process.selectedPatJets +
     process.lowPtJetFilter)
)

#                                         ##
process.source.fileNames = [          ##
    '/store/data/Run2011A/SingleElectron/AOD/PromptReco-v1/000/161/312/90646AF9-F957-E011-B0DB-003048F118C4.root'
    ] 
#                                         ##

process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )

#process.out.fileName = '/tmp/cleonido/patTuple.root'
