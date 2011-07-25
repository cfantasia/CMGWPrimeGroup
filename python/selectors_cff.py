import FWCore.ParameterSet.Config as cms
#Note: These cuts are <= and >= 
muonSelectors = cms.PSet(
    WZBase = cms.PSet(
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
    WZLoose = cms.PSet(),
    WZRelaxed = cms.PSet(),
    WZTight = cms.PSet(),
    exotica = cms.PSet(
    ),
    )
muonSelectors.WZLoose = muonSelectors.WZBase.clone()
muonSelectors.WZLoose.minPt = cms.untracked.double(10.)
muonSelectors.WZLoose.maxIso03 = cms.untracked.double(0.15)

muonSelectors.WZRelaxed = muonSelectors.WZBase.clone()
muonSelectors.WZRelaxed.minPt = cms.untracked.double(20.)

muonSelectors.WZTight = muonSelectors.WZBase.clone()
muonSelectors.WZTight.minPt = cms.untracked.double(20.)
muonSelectors.WZTight.maxIso03 = cms.untracked.double(0.1)

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

electronSelectors.WZTight = electronSelectors.wp80.clone()
electronSelectors.WZTight.barrel.minPt = cms.untracked.double(20.)
electronSelectors.WZTight.endcap.minPt = cms.untracked.double(20.)
