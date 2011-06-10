from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *
process.load('PhysicsTools.PFCandProducer.PF2PAT_cff')

from UserCode.CMGWPrimeGroup.patTuple_jet_cfg import *
from UserCode.CMGWPrimeGroup.patTuple_mumet_cfg import *
from UserCode.CMGWPrimeGroup.patTuple_mc_cfg import *

## remove MC matching from the default sequence when running on data
#removeMCMatching(process, ['All'])

# 2nd argument: message-logger frequency
# 3rd argument: # of events to process
mumet_config(process, 1, 5000)
jetExtra_config(process, 1, 5000)
mc_config(process, cms)

# keep all events with jet-pt above 30 GeV, |eta| < 2.4
process.selectedPatJets.cut = "pt > 30. & abs(eta) < 2.4"
process.lowPtJetFilter = cms.EDFilter("CandViewCountFilter",
                                src = cms.InputTag("selectedPatJets"),
                                minNumber = cms.uint32(1)
                                )

# Setup so that my patJets are PFJets and not calojets
process.patJets.jetSource='pfJets'
# Corrections
process.patJetCorrFactors.src = 'pfJets'
process.patJetCorrFactors.levels = cms.vstring('L2Relative',
                                               'L3Absolute')
process.patJetCorrFactors.payload = cms.string('AK5PF')

# Track association
process.jetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
                                                     jets = cms.InputTag("pfJets"),
                                                     tracks = cms.InputTag("generalTracks"),
                                                     coneSize = cms.double(0.5)
                                                     )

# turn off all the corrections and extra factors for now
process.patJets.addJetCorrFactors = True
process.patJets.addBTagInfo = False
process.patJets.addDiscriminators = False
process.patJets.addJetCharge = False
process.patJets.addJetID = False
process.patJets.addGenPartonMatch = False
process.patJets.embedGenPartonMatch = False
process.patJets.genPartonMatch = ''
process.patJets.addGenJetMatch = False
process.patJets.embedGenJetMatch = False
process.patJets.getJetMCFlavour = False
process.patJets.JetPartonMapSource = '' 
process.patJets.embedPFCandidates = True
process.patJets.addAssociatedTracks = True
process.patJets.trackAssociationSource = "jetTracksAssociatorAtVertex"
    

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
    process.pfJetSequence +
    process.jetTracksAssociatorAtVertex
    )

## let it run
process.p = cms.Path(
    process.muonMatch + 
    process.patMuons + 
    process.selectedPatMuons +
    process.PF2PATmod *
    (process.patJetCorrFactors + 
     process.patJets +
     process.selectedPatJets +
     process.lowPtJetFilter)
)

#                                         ##
process.source.fileNames = [          ##
    'file:/afs/cern.ch/user/t/tomei/public/PYTHIA6_EXOTICA_RSGravZZ_kMpl005_M1000_7TeV_mumujj_cff_py_GEN_FASTSIM_HLT.root',
    ] 
#                                         ##

process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )

#process.out.fileName = '/tmp/cleonido/patTuple.root'
