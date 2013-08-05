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
    PFIsoHighPtBoostedZRelaxed = cms.PSet(
       minPt = cms.untracked.double(10.),
       maxEta = cms.untracked.double(2.4),
       minIsGblOrTrk = cms.untracked.int32(1),
       maxDxy = cms.untracked.double(0.2),
       maxDz = cms.untracked.double(0.5),
       maxRelPtErr = cms.untracked.double(0.3),
       minNPixelHits = cms.untracked.int32(1),
       minNTrackerLayers = cms.untracked.int32(6),
       minNMatches = cms.untracked.int32(2),
       #maxPFIso = cms.untracked.double(???),#Added below in clones
    ),
    PFIsoHighPtRelaxed = cms.PSet(
       minPt = cms.untracked.double(10.),
       maxEta = cms.untracked.double(2.4),
       minIsGlobal = cms.untracked.int32(1),
       maxDxy = cms.untracked.double(0.2),
       maxDz = cms.untracked.double(0.5),
       maxRelPtErr = cms.untracked.double(0.3),
       minNPixelHits = cms.untracked.int32(1),
       minNMuonHits = cms.untracked.int32(1), 
       minNTrackerLayers = cms.untracked.int32(6),
       minNMatches = cms.untracked.int32(2),
       #maxPFIso = cms.untracked.double(???),#Added below in clones
    ),
    PFIsoRelaxed = cms.PSet(
       minPt = cms.untracked.double(10.),
       maxEta = cms.untracked.double(2.4),
       minIsGlobal = cms.untracked.int32(1),
       minIsPF = cms.untracked.int32(1),
       maxDxy = cms.untracked.double(0.2),
       maxDz = cms.untracked.double(0.5),
       maxNormalizedChi2 = cms.untracked.double(10.0), 
       maxRelPtErr = cms.untracked.double(0.3),
       minNPixelHits = cms.untracked.int32(1),
       minNMuonHits = cms.untracked.int32(1), 
       minNTrackerLayers = cms.untracked.int32(6),
       minNMatches = cms.untracked.int32(2),
       #maxPFIso = cms.untracked.double(???),#Added below in clones
    ),
    Reco = cms.PSet(
       minPt = cms.untracked.double(10.),
    ),
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
muonSelectors.WZRelaxed = muonSelectors.VBTF.clone(
    )
muonSelectors.WZRelaxedPt20 = muonSelectors.WZRelaxed.clone(
    minPt = 20.
    )
muonSelectors.WZLoose = muonSelectors.WZRelaxed.clone(
    maxIso03 = cms.untracked.double(0.15)
    )
muonSelectors.WZTight = muonSelectors.WZRelaxed.clone(
    maxIso03 = cms.untracked.double(0.1)
    )
muonSelectors.WZTightPt20 = muonSelectors.WZTight.clone(
    minPt = 20.,
    )
###EWKWZ
muonSelectors.EWKWZRelaxed = muonSelectors.PFIsoRelaxed.clone(
    )
muonSelectors.EWKWZRelaxedPt20 = muonSelectors.EWKWZRelaxed.clone(
    minPt = 20.
    )
muonSelectors.EWKWZLoose = muonSelectors.EWKWZRelaxed.clone(
    maxPFIso = cms.untracked.double(0.2)
    )
muonSelectors.EWKWZLoosePt20 = muonSelectors.EWKWZLoose.clone(
    minPt = 20.,
    )
muonSelectors.EWKWZTight = muonSelectors.EWKWZRelaxed.clone(
    maxPFIso = cms.untracked.double(0.12)
    )
muonSelectors.EWKWZTightPt20 = muonSelectors.EWKWZTight.clone(
    minPt = 20.,
    )

###High Pt PfIso Boosted Z
muonSelectors.PFIsoHighPtBoostedZRelaxedPt20 = muonSelectors.PFIsoHighPtBoostedZRelaxed.clone(
    minPt = 20.
    )
muonSelectors.PFIsoHighPtBoostedZLoose = muonSelectors.PFIsoHighPtBoostedZRelaxed.clone(
    maxPFIso = cms.untracked.double(0.2)
    )
muonSelectors.PFIsoHighPtBoostedZLoosePt20 = muonSelectors.PFIsoHighPtBoostedZLoose.clone(
    minPt = 20.,
    )
muonSelectors.PFIsoHighPtBoostedZTight = muonSelectors.PFIsoHighPtBoostedZRelaxed.clone(
    maxPFIso = cms.untracked.double(0.12)
    )
muonSelectors.PFIsoHighPtBoostedZTightPt20 = muonSelectors.PFIsoHighPtBoostedZTight.clone(
    minPt = 20.,
    )

#####High Pt PfIso
muonSelectors.PFIsoHighPtRelaxedPt20 = muonSelectors.PFIsoHighPtRelaxed.clone(
    minPt = 20.
    )
muonSelectors.PFIsoHighPtLoose = muonSelectors.PFIsoHighPtRelaxed.clone(
    maxPFIso = cms.untracked.double(0.2)
    )
muonSelectors.PFIsoHighPtLoosePt20 = muonSelectors.PFIsoHighPtLoose.clone(
    minPt = 20.,
    )
muonSelectors.PFIsoHighPtTight = muonSelectors.PFIsoHighPtRelaxed.clone(
    maxPFIso = cms.untracked.double(0.12)
    )
muonSelectors.PFIsoHighPtTightPt20 = muonSelectors.PFIsoHighPtTight.clone(
    minPt = 20.,
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
    WZRelaxed95 = cms.PSet(),
    WZRelaxed80 = cms.PSet(),
    WZTight = cms.PSet(),
    exotica = cms.PSet(
       joint = cms.PSet(
          minPt = cms.untracked.double(10.),
          minIsEcalDriven  = cms.untracked.int32(1),
          maxDeltaPhi = cms.untracked.double(0.06),
          maxHoverE = cms.untracked.double(0.05),
          #maxTrackIso = cms.untracked.double(5.),
          maxPFIso = cms.untracked.double(0.15),          
          maxMissingHits = cms.untracked.int32(1),
       ),
       barrel = cms.PSet(
          maxDeltaEta = cms.untracked.double(0.005),
          #maxSigmaIEtaIEta = cms.untracked.double(),
          minPassEX5overE55 = cms.untracked.int32(1),
          #minPassEMHadDepth1Iso = cms.untracked.int32(1),
          maxd0 = cms.untracked.double(0.02),
       ),
       endcap = cms.PSet(
          maxDeltaEta = cms.untracked.double(0.007),
          maxSigmaIEtaIEta = cms.untracked.double(0.03),
          minPassEX5overE55 = cms.untracked.int32(1),
          #minPassEMHadDepth1Iso = cms.untracked.int32(1),
          maxd0 = cms.untracked.double(0.05),
       ),
    ),
    mva = cms.PSet(
       joint = cms.PSet(
          minPt = cms.untracked.double(10.),
          minpassMVAPresel = cms.untracked.int32(1),
          minpassMVATrig = cms.untracked.int32(1),
          #maxPFIso = cms.untracked.double(???),#Added below in clones
       ),
       barrel = cms.PSet(
       ),
       endcap = cms.PSet(
       ),
    ),
    EWKWZLoose = cms.PSet(),
    EWKWZTight = cms.PSet(),

    CiC2012VetoRelaxed = cms.PSet(
       joint = cms.PSet(
          minPt = cms.untracked.double(10.),
          maxIsGap = cms.untracked.int32(0),
          minPassConvVeto = cms.untracked.int32(1),
          maxIsSCGap = cms.untracked.int32(0),
       ),
       barrel = cms.PSet(
           #maxMissingHits = cms.untracked.int32(1),
           #minVtxFitProb = cms.untracked.double(1e-6),
           maxSigmaIEtaIEta = cms.untracked.double(0.01),
           maxDeltaPhi = cms.untracked.double(0.8),
           maxDeltaEta = cms.untracked.double(0.007),
           maxHoverE = cms.untracked.double(0.15),
           maxd0 = cms.untracked.double(0.04),
           maxdz = cms.untracked.double(0.2),
           #maxfabsdiffEp = cms.untracked.double(0.05),
       ),
       endcap = cms.PSet(
           #maxMissingHits = cms.untracked.int32(1),
           #minVtxFitProb = cms.untracked.double(1e-6),
           maxSigmaIEtaIEta = cms.untracked.double(0.03),
           maxDeltaPhi = cms.untracked.double(0.7),
           maxDeltaEta = cms.untracked.double(0.01),
           #maxHoverE = cms.untracked.double(0.10),
           maxd0 = cms.untracked.double(0.04),
           maxdz = cms.untracked.double(0.2),
           #maxfabsdiffEp = cms.untracked.double(0.05),
       ),
    ),
    CiC2012LooseRelaxed = cms.PSet(
       joint = cms.PSet(
          minPt = cms.untracked.double(10.),
          maxIsGap = cms.untracked.int32(0),
          minPassConvVeto = cms.untracked.int32(1),
          maxIsSCGap = cms.untracked.int32(0),
       ),
       barrel = cms.PSet(
           maxMissingHits = cms.untracked.int32(1),
           minVtxFitProb = cms.untracked.double(1e-6),
           maxSigmaIEtaIEta = cms.untracked.double(0.01),
           maxDeltaPhi = cms.untracked.double(0.15),
           maxDeltaEta = cms.untracked.double(0.007),
           maxHoverE = cms.untracked.double(0.12),
           maxd0 = cms.untracked.double(0.02),
           maxdz = cms.untracked.double(0.2),
           maxfabsdiffEp = cms.untracked.double(0.05),
       ),
       endcap = cms.PSet(
           maxMissingHits = cms.untracked.int32(1),
           minVtxFitProb = cms.untracked.double(1e-6),
           maxSigmaIEtaIEta = cms.untracked.double(0.03),
           maxDeltaPhi = cms.untracked.double(0.10),
           maxDeltaEta = cms.untracked.double(0.009),
           maxHoverE = cms.untracked.double(0.10),
           maxd0 = cms.untracked.double(0.02),
           maxdz = cms.untracked.double(0.2),
           maxfabsdiffEp = cms.untracked.double(0.05),
       ),
    ),
    Reco = cms.PSet(
       joint = cms.PSet(
          minPt = cms.untracked.double(10.),
          maxIsGap = cms.untracked.int32(0),
          minPassConvVeto = cms.untracked.int32(1),
          maxIsSCGap = cms.untracked.int32(0),
       ),
       barrel = cms.PSet(),
       endcap = cms.PSet(),
    ),
    CiC2012Loose = cms.PSet(),
    CiC2012LoosePt20 = cms.PSet(),
    CiC2012MediumRelaxed = cms.PSet(
       joint = cms.PSet(
          minPt = cms.untracked.double(10.),
          maxIsGap = cms.untracked.int32(0),
          minPassConvVeto = cms.untracked.int32(1),
          maxIsSCGap = cms.untracked.int32(0),
       ),
       barrel = cms.PSet(
           maxMissingHits = cms.untracked.int32(1),
           minVtxFitProb = cms.untracked.double(1e-6),
           maxSigmaIEtaIEta = cms.untracked.double(0.01),
           maxDeltaPhi = cms.untracked.double(0.06),
           maxDeltaEta = cms.untracked.double(0.004),
           maxHoverE = cms.untracked.double(0.12),
           maxd0 = cms.untracked.double(0.02),
           maxdz = cms.untracked.double(0.1),
           maxfabsdiffEp = cms.untracked.double(0.05),
       ),
       endcap = cms.PSet(
           maxMissingHits = cms.untracked.int32(1),
           minVtxFitProb = cms.untracked.double(1e-6),
           maxSigmaIEtaIEta = cms.untracked.double(0.03),
           maxDeltaPhi = cms.untracked.double(0.03),
           maxDeltaEta = cms.untracked.double(0.007),
           maxHoverE = cms.untracked.double(0.10),
           maxd0 = cms.untracked.double(0.02),
           maxdz = cms.untracked.double(0.1),
           maxfabsdiffEp = cms.untracked.double(0.05),
       ),
    ),
    CiC2012Medium = cms.PSet(),
    CiC2012MediumPt20 = cms.PSet(),
    )

for i, s in enumerate(["wp95", "wp90", "wp85", "wp80", "wp70", "wp60"]):
    pset = cms.PSet(
        joint = cms.PSet(
           minPt = cms.untracked.double(10.),
           maxIsSCGap = cms.untracked.int32(0),
        ),
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

electronSelectors.WZRelaxed95 = electronSelectors.wp95.clone()
electronSelectors.WZRelaxed95.barrel.maxCombRelIso = cms.untracked.double(999999999.)
electronSelectors.WZRelaxed95.endcap.maxCombRelIso = cms.untracked.double(999999999.)

electronSelectors.WZRelaxed80 = electronSelectors.wp80.clone()
electronSelectors.WZRelaxed80.barrel.maxCombRelIso = cms.untracked.double(999999999.)
electronSelectors.WZRelaxed80.endcap.maxCombRelIso = cms.untracked.double(999999999.)

electronSelectors.WZRelaxed80Pt20 = electronSelectors.wp80.clone()
electronSelectors.WZRelaxed80Pt20.joint.minPt = cms.untracked.double(20.)
electronSelectors.WZRelaxed80Pt20.barrel.maxCombRelIso = cms.untracked.double(999999999.)
electronSelectors.WZRelaxed80Pt20.endcap.maxCombRelIso = cms.untracked.double(999999999.)

electronSelectors.WZTight = electronSelectors.wp80.clone()

electronSelectors.WZTightPt20 = electronSelectors.wp80.clone()
electronSelectors.WZTightPt20.joint.minPt = cms.untracked.double(20.)

##EWKWZ
electronSelectors.EWKWZLoose = electronSelectors.mva.clone()
electronSelectors.EWKWZLoose.joint.maxPFIso = cms.untracked.double(0.1)

electronSelectors.EWKWZRelaxed = electronSelectors.mva.clone()

electronSelectors.EWKWZRelaxedPt20 = electronSelectors.mva.clone()
electronSelectors.EWKWZRelaxedPt20.joint.minPt = 20.

electronSelectors.EWKWZTight = electronSelectors.mva.clone()
electronSelectors.EWKWZTight.joint.maxPFIso = cms.untracked.double(0.1)

electronSelectors.EWKWZTightPt20 = electronSelectors.mva.clone()
electronSelectors.EWKWZTightPt20.joint.minPt = 20.
electronSelectors.EWKWZTightPt20.joint.maxPFIso = cms.untracked.double(0.1)

#CiC2012
electronSelectors.CiC2012Veto = electronSelectors.CiC2012VetoRelaxed.clone()
electronSelectors.CiC2012Veto.joint.maxPFIso = cms.untracked.double(0.15)

electronSelectors.CiC2012VetoPt20 = electronSelectors.CiC2012Veto.clone()
electronSelectors.CiC2012VetoPt20.joint.minPt = 20.

electronSelectors.CiC2012Loose = electronSelectors.CiC2012LooseRelaxed.clone()
electronSelectors.CiC2012Loose.joint.maxPFIso = cms.untracked.double(0.15)

electronSelectors.CiC2012LoosePt20 = electronSelectors.CiC2012Loose.clone()
electronSelectors.CiC2012LoosePt20.joint.minPt = 20.
    
electronSelectors.CiC2012MediumRelaxedPt20 = electronSelectors.CiC2012MediumRelaxed.clone()
electronSelectors.CiC2012MediumRelaxedPt20.joint.minPt = 20.

electronSelectors.CiC2012Medium = electronSelectors.CiC2012MediumRelaxed.clone()
electronSelectors.CiC2012Medium.joint.maxPFIso = cms.untracked.double(0.15)
#electronSelectors.CiC2012Medium.joint.maxTrackIso = cms.untracked.double(2)

electronSelectors.CiC2012MediumPt20 = electronSelectors.CiC2012Medium.clone()
electronSelectors.CiC2012MediumPt20.joint.minPt = 20.

##Exotica
electronSelectors.exoticaPt20 = electronSelectors.exotica.clone()
electronSelectors.exoticaPt20.joint.minPt = cms.untracked.double(20.)

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

#CiC2012Relaxed
