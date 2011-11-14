import FWCore.ParameterSet.Config as cms
from UserCode.CMGWPrimeGroup.patTuple_common_cfg import *
#from PhysicsTools.PFCandProducer.PF2PAT_cff import *
#from CommonTools.ParticleFlow.PF2PAT_cff import *
from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *
process.load("PhysicsTools.PatAlgos.patSequences_cff")



def CMGWPswitchToPFJets(process) :

    addJetCollection(process,cms.InputTag('ak7PFJets'),'AK7','PF',
                     doJTA        = True,
                     doBTagging   = True,
#                     addDiscriminators    = True,   ## addition btag discriminators                    
#                     discriminatorSources = cms.VInputTag(
#            cms.InputTag("combinedSecondaryVertexBJetTags"),
#            cms.InputTag("combinedSecondaryVertexMVABJetTags"),
#            cms.InputTag("jetBProbabilityBJetTags"),
#            cms.InputTag("jetProbabilityBJetTags"),
#            cms.InputTag("simpleSecondaryVertexHighEffBJetTags"),
#            cms.InputTag("simpleSecondaryVertexHighPurBJetTags"),
#            cms.InputTag("softElectronByPtBJetTags"),
#            cms.InputTag("softElectronByIP3dBJetTags"),
#            cms.InputTag("softMuonBJetTags"),
#            cms.InputTag("softMuonByPtBJetTags"),
#            cms.InputTag("softMuonByIP3dBJetTags"),
#            cms.InputTag("trackCountingHighEffBJetTags"),
#            cms.InputTag("trackCountingHighPurBJetTags"),
#            ),
                     jetCorrLabel = ('AK7PF', cms.vstring(['L1Offset', 'L2Relative', 'L3Absolute'])),
                     doType1MET   = True,
                     genJetCollection=cms.InputTag("ak7GenJets"),
                     doJetID      = False
                     )
    


    addJetCollection(process,cms.InputTag('ak5PFJets'),'AK5','PF',
                                          doJTA        = True,
                                          doBTagging   = True,
                     #                     addDiscriminators    = True,   ## addition btag discriminators
                     #                     discriminatorSources = cms.VInputTag(
                     #            cms.InputTag("combinedSecondaryVertexBJetTags"),
                     #            cms.InputTag("combinedSecondaryVertexMVABJetTags"),
                     #            cms.InputTag("jetBProbabilityBJetTags"),
                     #            cms.InputTag("jetProbabilityBJetTags"),
                     #            cms.InputTag("simpleSecondaryVertexHighEffBJetTags"),
                     #            cms.InputTag("simpleSecondaryVertexHighPurBJetTags"),
                     #            cms.InputTag("softElectronByPtBJetTags"),
                     #            cms.InputTag("softElectronByIP3dBJetTags"),
                     #            cms.InputTag("softMuonBJetTags"),
                     #            cms.InputTag("softMuonByPtBJetTags"),
                     #            cms.InputTag("softMuonByIP3dBJetTags"),
                     #            cms.InputTag("trackCountingHighEffBJetTags"),
                     #            cms.InputTag("trackCountingHighPurBJetTags"),
                     #            ),
                                          jetCorrLabel = ('AK5PF', cms.vstring(['L1Offset', 'L2Relative', 'L3Absolute'])),
                                          doType1MET   = True,
                                          genJetCollection=cms.InputTag("ak5GenJets"),
                                          doJetID      = False
                                          )

    

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
    # keep all events with jet-pt above 30 GeV, |eta| < 2.4
    process.selectedPatJets.cut = "pt > 30. & abs(eta) < 2.4"

    # RECO
    process.out.outputCommands.append('keep *_selectedPatJetsAK7PF_*_*')
    process.out.outputCommands.append('keep *_ak7GenJets_*_*')


    process.out.outputCommands.append('keep *_selectedPatJetsAK5PF_*_*')
    process.out.outputCommands.append('keep *_ak5GenJets_*_*')
        

# Modules and sequences

# Track association
jetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
                                             jets = cms.InputTag("pfJets"),
                                             tracks = cms.InputTag("generalTracks"),
                                             coneSize = cms.double(0.5)
                                             )




PF2PATmod = cms.Sequence(process.patDefaultSequence)

# Adaptation of PF2PAT
#PF2PATmod = cms.Sequence(pfNoPileUpSequence +
#                         pfAllNeutralHadrons +  
#                         pfAllChargedHadrons +
#                         pfAllPhotons +
#                         pfMuonSequence +
#                         pfNoMuon +
#                         pfElectronSequence +
#                         pfNoElectron +
#                         pfJetSequence +
#                         jetTracksAssociatorAtVertex
#                         )
