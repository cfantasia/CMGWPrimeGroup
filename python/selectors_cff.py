import FWCore.ParameterSet.Config as cms

#####################
######  Muons  ######
#####################

#Note: These cuts are <= and >= 
muonSelectors = cms.PSet(
    VBTF = cms.PSet(
       minPt = cms.untracked.double(10.),
       maxEta = cms.untracked.double(2.4),
       minIsGlobal = cms.untracked.int32(1),
       maxDxy = cms.untracked.double(0.2),
       maxNormalizedChi2 = cms.untracked.double(10.0),
       minNTrackerHits = cms.untracked.int32(11),
       minNPixelHits = cms.untracked.int32(1),
       minNMuonHits = cms.untracked.int32(1),
       minNMatches = cms.untracked.int32(2),
    ),
    PFIso = cms.PSet(
       minPt = cms.untracked.double(10.),
       maxEta = cms.untracked.double(2.4),
       minIsGlobal = cms.untracked.int32(1),
       maxDxy = cms.untracked.double(0.2),
       maxNormalizedChi2 = cms.untracked.double(10.0),
       minNPixelHits = cms.untracked.int32(1),
       minNMuonHits = cms.untracked.int32(1),
       minNTrackerLayers = cms.untracked.int32(9),
       minNMatches = cms.untracked.int32(2),
       #maxPFIso = cms.untracked.double(???),#Added below in clones
    ),
    WZLoose = cms.PSet(),
    WZRelaxed = cms.PSet(),
    WZTight = cms.PSet(),
    HadVZLoose = cms.PSet(),
    HadVZTight = cms.PSet(),
    exotica = cms.PSet(
       maxEta = cms.untracked.double(2.4),
       minIsGlobal = cms.untracked.int32(1),
       minIsTracker = cms.untracked.int32(1),
       maxDxy = cms.untracked.double(0.2),
#       maxNormalizedChi2 = cms.untracked.double(10.), not needed for TeV muons
#       minNTrackerHits = cms.untracked.int32(11), deprecated
       minNPixelHits = cms.untracked.int32(1),
       minNMuonHits = cms.untracked.int32(1),
       minNMatches = cms.untracked.int32(2),
       minNTrackerLayers = cms.untracked.int32(9),
#       minTrackerValidFrac = cms.untracked.double(),
    ),
    )
###WprimeWZ
muonSelectors.WZLoose = muonSelectors.VBTF.clone(
    minPt = 10.,
    maxIso03 = cms.untracked.double(0.15)
    )
muonSelectors.WZRelaxed = muonSelectors.VBTF.clone(
    minPt = 20.
    )
muonSelectors.WZTight = muonSelectors.VBTF.clone(
    minPt = 20.,
    maxIso03 = cms.untracked.double(0.1)
    )
###EWKWZ
muonSelectors.EWKWZLoose = muonSelectors.PFIso.clone(
    minPt = 10.,
    maxPFIso = cms.untracked.double(0.2)
    )
muonSelectors.EWKWZRelaxed = muonSelectors.PFIso.clone(
    minPt = 20.
    )
muonSelectors.EWKWZTight = muonSelectors.PFIso.clone(
    minPt = 20.,
    maxPFIso = cms.untracked.double(0.12)
    )
###HadVZ
muonSelectors.HadVZLoose = muonSelectors.exotica.clone(
    minPt = cms.untracked.double(20.),
#    maxIso03 = cms.untracked.double(0.1)
    
    )
muonSelectors.HadVZTight = muonSelectors.exotica.clone(
    minPt = cms.untracked.double(35.),
#    maxIso03 = cms.untracked.double(0.1)
    
    )
#muonSelectors.HadVZTight.remove(maxNormalizedChi2)

#####################
####  Electrons  ####
#####################

cutsMissingHits = [0, 0, 0, 0, 0, 0]
cutsConvDist = [0., 0., 0.02, 0.02, 0.02, 0.02]
cutsConvDcot = [0., 0., 0.02, 0.02, 0.02, 0.02]
cutsEBCombRelIso = [0.150, 0.085, 0.053, 0.070, 0.030, 0.016]#2010
#cutsEBCombRelIso = [0.150, 0.085, 0.053, 0.040, 0.030, 0.016]#2011
cutsEBSigmaIEtaIEta = [0.012, 0.01, 0.01, 0.01, 0.01, 0.01]
cutsEBDeltaPhi = [0.800, 0.071, 0.039, 0.027, 0.020, 0.020]
cutsEBDeltaEta = [0.007, 0.007, 0.005, 0.005, 0.004, 0.004]
#cutsEECombRelIso = [0.100, 0.051, 0.042, 0.033, 0.016, 0.008]#2011
cutsEECombRelIso = [0.100, 0.051, 0.042, 0.060, 0.016, 0.008]#2010
cutsEESigmaIEtaIEta = [0.031, 0.031, 0.031, 0.031, 0.031, 0.031]
cutsEEDeltaPhi = [0.7, 0.047, 0.028, 0.021, 0.021, 0.021]
cutsEEDeltaEta = [0.011, 0.011, 0.007, 0.006, 0.005, 0.004]

electronSelectors = cms.PSet(
    WZLoose = cms.PSet(),
    WZRelaxed = cms.PSet(),
    WZTight = cms.PSet(),
    exotica = cms.PSet(),
    )

for i, s in enumerate(["wp95", "wp90", "wp85", "wp80", "wp70", "wp60"]):
    pset = cms.PSet(
        barrel = cms.PSet(
           maxMissingHits = cms.untracked.int32(cutsMissingHits[i]),
           minConv = cms.untracked.double(cutsConvDist[i]),#Hack bc we need an OR of these two cuts below
           maxCombRelIso = cms.untracked.double(cutsEBCombRelIso[i]),
           maxSigmaIEtaIEta = cms.untracked.double(cutsEBSigmaIEtaIEta[i]),
           maxDeltaPhi = cms.untracked.double(cutsEBDeltaPhi[i]),
           maxDeltaEta = cms.untracked.double(cutsEBDeltaEta[i]),
           ),
        endcap = cms.PSet(
           maxMissingHits = cms.untracked.int32(cutsMissingHits[i]),
           minConv = cms.untracked.double(cutsConvDist[i]),#Hack bc we need an or of these two cuts below
           maxCombRelIso = cms.untracked.double(cutsEECombRelIso[i]),
           maxSigmaIEtaIEta = cms.untracked.double(cutsEESigmaIEtaIEta[i]),
           maxDeltaPhi = cms.untracked.double(cutsEEDeltaPhi[i]),
           maxDeltaEta = cms.untracked.double(cutsEEDeltaEta[i]),
           ),
        )
    setattr(electronSelectors, s, pset)
    
electronSelectors.WZLoose = electronSelectors.wp95.clone()
electronSelectors.WZLoose.barrel.minPt = cms.untracked.double(10.)
electronSelectors.WZLoose.endcap.minPt = cms.untracked.double(10.)

electronSelectors.WZRelaxed = electronSelectors.wp80.clone()
electronSelectors.WZRelaxed.barrel.minPt = cms.untracked.double(20.)
electronSelectors.WZRelaxed.endcap.minPt = cms.untracked.double(20.)
electronSelectors.WZRelaxed.barrel.maxCombRelIso = cms.untracked.double(999999999.)
electronSelectors.WZRelaxed.endcap.maxCombRelIso = cms.untracked.double(999999999.)
#electronSelectors.WZRelaxed.barrel.remove(maxCombRelIso)
#electronSelectors.WZRelaxed.endcap.remove(maxCombRelIso)

electronSelectors.WZTight = electronSelectors.wp80.clone()
electronSelectors.WZTight.barrel.minPt = cms.untracked.double(20.)
electronSelectors.WZTight.endcap.minPt = cms.untracked.double(20.)


####################
#####  Jets  #######
####################
jetSelectors = cms.PSet(
    Base = cms.PSet(
       minPt = cms.untracked.double(30.),
       maxEta = cms.untracked.double(2.4),
       maxNHF = cms.untracked.double(0.99),
       maxNEF = cms.untracked.double(0.99),
       minNDaughters = cms.untracked.int32(1),
       minCHF = cms.untracked.double(0.0),
       maxCEF = cms.untracked.double(0.99),
       minCMult = cms.untracked.int32(0),
       ),
    Pat = cms.PSet(
       minPt = cms.untracked.double(30.),
       maxEta = cms.untracked.double(3.0),
       ),
    )

####################
#### Photons #######
####################
photonSelectors = cms.PSet(
    Base = cms.PSet(#Made up values
       minPt = cms.untracked.double(10.),
       maxEta = cms.untracked.double(5),
       maxECalIso = cms.untracked.double(0.1),
       maxHCalIso = cms.untracked.double(0.1),
       maxTrkIso = cms.untracked.double(0.1),
       maxHoE = cms.untracked.double(0.05),
       maxSigmaee = cms.untracked.double(9e9),
       minHasSeed = cms.untracked.bool(True),       
       ),
    )

