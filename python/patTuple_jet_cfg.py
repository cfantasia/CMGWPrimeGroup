import FWCore.ParameterSet.Config as cms
from UserCode.CMGWPrimeGroup.patTuple_common_cfg import *
#from PhysicsTools.PFCandProducer.PF2PAT_cff import *
from CommonTools.ParticleFlow.PF2PAT_cff import *
from PhysicsTools.PatAlgos.tools.pfTools import *


def CMGWPswitchToPFJets(process) :
    ## Trying to change Jet algorithm
    process.pfJets = jetAlgo('AK7')


    # IMPORTANT: must have patTemplate loaded
    # Setup so that my patJets are PFJets and not calojets
    process.patJets.jetSource='pfJets'
    
    # Corrections
    process.patJetCorrFactors.src = 'pfJets'
    process.patJetCorrFactors.levels = cms.vstring('L2Relative',
                                                   'L3Absolute')
    process.patJetCorrFactors.payload = cms.string('AK7PF')
    
    # Turn off other extra factors
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

    process.selectedPatPFParticles.cut = ""


    from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
    process.goodOfflinePrimaryVertices = cms.EDFilter(
        "PrimaryVertexObjectFilter",
        filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
        src=cms.InputTag('offlinePrimaryVertices')
        )

    
def jetFull_config(process, reportEveryNum=100, maxEvents=-1) :
    common_config(process, reportEveryNum, maxEvents)
    jet_config(process, reportEveryNum, maxEvents)
       
def jet_config(process, reportEveryNum=100, maxEvents=-1) :
    CMGWPswitchToPFJets(process)
    # RECO
    process.out.outputCommands.append('keep *_selectedPatJets_*_*')
    process.out.outputCommands.append('keep *_patJets_*_*')
    
# Modules and sequences

# Track association
jetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
                                             jets = cms.InputTag("pfJets"),
                                             tracks = cms.InputTag("generalTracks"),
                                             coneSize = cms.double(0.5)
                                             )

# Adaptation of PF2PAT
PF2PATmod = cms.Sequence(pfNoPileUpSequence +
                         pfAllNeutralHadrons +  
                         pfAllChargedHadrons +
                         pfAllPhotons +
                         pfMuonSequence +
                         pfNoMuon +
                         pfElectronSequence +
                         pfNoElectron +
                         pfJetSequence +
                         jetTracksAssociatorAtVertex
                         )
