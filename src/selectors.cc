#include "UserCode/CMGWPrimeGroup/interface/selectors.h"

ElectronSelectorBase::ElectronSelectorBase() {}
ElectronSelectorBase::ElectronSelectorBase(Pset const params) {
  // set the last parameter to false to turn off the cut
  loadFromPset<double>(params, "minPt", true);
  loadFromPset<double>(params, "minConv", true);
  loadFromPset<double>(params, "maxSigmaIEtaIEta", true);
  loadFromPset<double>(params, "maxDeltaEta", true);
  loadFromPset<double>(params, "maxDeltaPhi", true);
  loadFromPset<double>(params, "maxCombRelIso", true);
  loadFromPset<int>(params, "maxMissingHits", true);
  loadFromPset<int>(params, "minpassMVAPresel", true);
  loadFromPset<int>(params, "minpassMVATrig", true);
  loadFromPset<int>(params, "minIsEcalDriven", true);
  loadFromPset<int>(params, "minPassEX5overE55", true);
  loadFromPset<int>(params, "minPassEMHadDepth1Iso", true);
  loadFromPset<double>(params, "maxPFIso", true);
  loadFromPset<double>(params, "maxHoverE", true);
  loadFromPset<double>(params, "maxd0", true);
  loadFromPset<double>(params, "maxdz", true);
  loadFromPset<double>(params, "maxfabsdiffEp", true);
  loadFromPset<double>(params, "maxTrackIso", true);

  bitmask = getBitTemplate();
}

bool ElectronSelectorBase::operator()(const pat::Electron & p, pat::strbitset & ret) {
  bool val = (*this)(p,0.);
  ret = bitmask;
  return val;
}
bool ElectronSelectorBase::operator()(const pat::Electron & p, const float pu, const reco::Vertex & vtx) {
  bitmask.set(false);
  setpassCut("minPt", p.pt(), bitmask);
  if(ignoreCut("minConv") || 
     fabs(p.convDist()) >= cut("minConv", double()) || fabs(p.convDcot()) >= cut("minConv", double()))
    passCut(bitmask, "minConv");
  setpassCut("maxSigmaIEtaIEta", p.sigmaIetaIeta(), bitmask);
  setpassCut("maxDeltaEta", fabs(p.deltaEtaSuperClusterTrackAtVtx()), bitmask);
  setpassCut("maxDeltaPhi", fabs(p.deltaPhiSuperClusterTrackAtVtx()), bitmask);   
  setpassCut("maxCombRelIso",calcCombRelIso(p, pu), bitmask);
  setpassCut("maxMissingHits", 
             p.gsfTrack()->trackerExpectedHitsInner().numberOfHits(), bitmask);
  if(ignoreCut("minpassMVAPresel") || p.userInt("gsf_passPreselMVA2011") >= cut("minpassMVAPresel", int())) passCut(bitmask, "minpassMVAPresel");
  if(ignoreCut("minpassMVATrig") || p.userInt("gsf_BDT_MVA_HWW2011_Trig_pass") >= cut("minpassMVATrig", int())) passCut(bitmask, "minpassMVATrig");

  if(ignoreCut("maxPFIso") || (pfIso(p,pu) <= cut("maxPFIso", double()))) passCut(bitmask, "maxPFIso");
  setpassCut("maxHoverE", p.hadronicOverEm(), bitmask);
  setpassCut("maxd0", fabs(p.dB()), bitmask);
  //setpassCut("maxd0", fabs(p.gsfTrack()->dxy(vtx.position())), bitmask);
  setpassCut("maxdz", fabs(p.gsfTrack()->dz (vtx.position())), bitmask);
  setpassCut("maxfabsdiffEp", fabs(1.0/p.ecalEnergy() - p.eSuperClusterOverP()/p.ecalEnergy()), bitmask);

  setpassCut("minIsEcalDriven", p.ecalDrivenSeed(), bitmask);
  setpassCut("maxTrackIso", p.dr03TkSumPt(), bitmask);
  if(ignoreCut("minPassEX5overE55") || (p.e2x5Max()  > p.e5x5()*0.94 || p.e1x5() > p.e5x5()*0.83)) passCut(bitmask, "minPassEX5overE55"); 
  double scet = p.caloEnergy()*sin(p.p4().theta());
  double EMHadDepthCut = p.isEB() ? (2 +0.03*scet + 0.28*pu) : scet < 50 ? (2.5 + 0.28*pu) : (2.5 + 0.03*(scet-50)+0.28*pu);
  if(ignoreCut("minPassEMHadDepth1Iso") || (p.dr03EcalRecHitSumEt() + p.dr03HcalDepth1TowerSumEt() < EMHadDepthCut) ) passCut(bitmask, "minPassEMHadDepth1Iso"); 
  
  setIgnored(bitmask); 
  return (bool) bitmask;
}
bool ElectronSelectorBase::operator()(const heep::Ele & p, pat::strbitset & ret) {
  bool val = (*this)(p,0.);
  ret = bitmask;
  return val;
}
bool ElectronSelectorBase::operator()(const heep::Ele & p, const float pu) {
  bitmask.set(false);
  const pat::Electron & e = p.patEle();
  setpassCut("minPt", p.et(), bitmask);
  if(ignoreCut("minConv") || 
     fabs(e.convDist()) >= cut("minConv", double()) || fabs(e.convDcot()) >= cut("minConv", double()))
    passCut(bitmask, "minConv");
  setpassCut("maxSigmaIEtaIEta", e.sigmaIetaIeta(), bitmask);
  setpassCut("maxDeltaEta", fabs(e.deltaEtaSuperClusterTrackAtVtx()), bitmask);
  setpassCut("maxDeltaPhi", fabs(e.deltaPhiSuperClusterTrackAtVtx()), bitmask);   
  setpassCut("maxCombRelIso",calcCombRelIso(e, pu), bitmask);
  setpassCut("maxMissingHits", 
             e.gsfTrack()->trackerExpectedHitsInner().numberOfHits(), bitmask);

  setIgnored(bitmask);
  return (bool) bitmask;
}
  
float ElectronSelectorBase::pfIso(const pat::Electron & p, const float & pu) const{
  return (p.chargedHadronIso() + std::max(0., p.neutralHadronIso() + p.photonIso()-0.5*p.puChargedHadronIso())) / p.pt();
}

/// A wrapper to handle the barrel/endcap split for electrons
ElectronSelector::ElectronSelector() {}

ElectronSelector::ElectronSelector(Pset pset, std::string selectorName) {
  Pset const params = pset.getParameter<Pset>(selectorName);
  //I want to combine psets, does this work
  ElectronSelectorBase jointSelector = ElectronSelectorBase(params.getParameter<Pset>("joint"));
  Pset jointPSet = params.getParameter<Pset>("joint");
  Pset barrelPSet = params.getParameter<Pset>("barrel");
  Pset endcapPSet = params.getParameter<Pset>("endcap");
  std::vector<std::string> jointParameterNames = jointPSet.getParameterNames();
  for ( std::vector<std::string>::iterator iter = jointParameterNames.begin();
        iter != jointParameterNames.end(); iter ++ ) {
    barrelPSet.copyFrom(jointPSet, (*iter));
    endcapPSet.copyFrom(jointPSet, (*iter));
  }
  
  barrelSelector_ = ElectronSelectorBase(barrelPSet);
  endcapSelector_ = ElectronSelectorBase(endcapPSet);
}

bool ElectronSelector::operator()(const pat::Electron & p, const float pu, const reco::Vertex & vtx) {
  if     (p.isEB()) return barrelSelector_(p, pu);
  else if(p.isEE()) return endcapSelector_(p, pu);
  return false;
}
bool ElectronSelector::operator()(const heep::Ele & p, const float pu) {
  if     (p.isEB()) return barrelSelector_(p, pu);
  else if(p.isEE()) return endcapSelector_(p, pu);
  return false;
}
pat::strbitset ElectronSelector::getBitTemplate() { return barrelSelector_.getBitTemplate(); }


/// Selector for muons based on params
MuonSelector::MuonSelector() {}
MuonSelector::MuonSelector(Pset pset, std::string selectorName) {
  Pset const params = pset.getParameter<Pset>(selectorName);
  // set the last parameter to false to turn off the cut
  loadFromPset<double>(params, "minPt", true);
  loadFromPset<double>(params, "maxEta", true);
  loadFromPset<double>(params, "maxDxy", true);
  loadFromPset<double>(params, "maxDz", true);
  loadFromPset<double>(params, "maxIso", true);
  loadFromPset<double>(params, "maxIso03", true);
  loadFromPset<double>(params, "maxTrkCorIso", true);
  loadFromPset<double>(params, "maxPFIso", true);
  loadFromPset<int>(params, "minIsGlobal", true);
  loadFromPset<int>(params, "minIsTracker", true);
  loadFromPset<int>(params, "minIsGblOrTrk", true);
  loadFromPset<int>(params, "minIsPF", true);
  loadFromPset<int>(params, "minNMatches", true);
  loadFromPset<int>(params, "minNMatchedStations", true);
  loadFromPset<double>(params, "maxNormalizedChi2", true);
  loadFromPset<int>(params, "minNTrackerHits", true);
  loadFromPset<int>(params, "minNPixelHits", true);
  loadFromPset<int>(params, "minNMuonHits", true);
  loadFromPset<int>(params, "minNTrackerLayers", true);
  loadFromPset<double>(params, "minTrackerValidFrac", true);

  bitmask = getBitTemplate();
}
bool MuonSelector::operator()(const TeVMuon & p, pat::strbitset & ret) {
  bool val = (*this)(p,0.);
  ret = bitmask;
  return val;
}
bool MuonSelector::operator()(const TeVMuon & p, const float pu, const reco::Vertex & vtx) {
  bitmask.set(false);
  setpassCut("minPt", p.pt(), bitmask);
  setpassCut("maxEta", fabs(p.eta()), bitmask);
  setpassCut("maxDxy", fabs(p.dB()), bitmask);
  setpassCut("maxIso", p.combRelIsolation(), bitmask);
  setpassCut("maxIso03", p.combRelIsolation03(pu), bitmask);
  setpassCut("maxTrkCorIso", p.trkCorRelIsolation(pu), bitmask);
  setpassCut("maxPFIso", p.combRelPFIsolation(), bitmask);
  setpassCut("minIsGlobal", p.isGlobalMuon(), bitmask);
  setpassCut("minIsTracker", p.isTrackerMuon(), bitmask);
  setpassCut("minIsGblOrTrk", p.isGlobalMuon() || p.isTrackerMuon(), bitmask);
  setpassCut("minIsPF", p.isPFMuon(), bitmask);
  setpassCut("minNMatches", p.numberOfMatches(), bitmask);
  setpassCut("minNMatchedStations", p.numberOfMatchedStations(), bitmask);
    
    reco::TrackRef global = p.globalTrack();
    if(!global.isNull()){
      const reco::HitPattern& gtHP = global->hitPattern();
      setpassCut("maxNormalizedChi2", global->normalizedChi2(), bitmask);
      setpassCut("minNMuonHits", gtHP.numberOfValidMuonHits(), bitmask);
      setpassCut("minTrackerValidFrac", global->validFraction(), bitmask);
      //delete below
      //setpassCut("minNTrackerHits", gtHP.numberOfValidTrackerHits(), bitmask);
      //setpassCut("minNPixelHits", gtHP.numberOfValidPixelHits(), bitmask);
      //setpassCut("minNTrackerLayers", gtHP.trackerLayersWithMeasurement(), bitmask);
      setpassCut("maxDz",  fabs(p.muonBestTrack()->dz(vtx.position())), bitmask);
    }

    reco::TrackRef inner = p.innerTrack();
    if(!inner.isNull()){
      const reco::HitPattern& inHP = inner->hitPattern();
      setpassCut("minNTrackerHits", inHP.numberOfValidTrackerHits(), bitmask);
      setpassCut("minNPixelHits", inHP.numberOfValidPixelHits(), bitmask);
      setpassCut("minNTrackerLayers", inHP.trackerLayersWithMeasurement(), bitmask);
      //setpassCut("maxDz",  fabs(p.innerTrack()->dxy(vtx.position())), bitmask);
    }


  setIgnored(bitmask);
  return (bool) bitmask;
}




/// Selector for jets based on params
JetSelector::JetSelector() {}
JetSelector::JetSelector(Pset pset, std::string selectorName) {
  Pset const params = pset.getParameter<Pset>(selectorName);
  // set the last parameter to false to turn off the cut
  loadFromPset<double>(params, "minPt", true);
  loadFromPset<double>(params, "maxEta", true);
  loadFromPset<double>(params, "maxNHF", true);
  loadFromPset<double>(params, "maxNEF", true);
  loadFromPset<int>(params, "minNDaughters", true);
  loadFromPset<double>(params, "minCHF", true);
  loadFromPset<double>(params, "maxCEF", true);
  loadFromPset<int>(params, "minCMult", true);

  bitmask = getBitTemplate();
}
bool JetSelector::operator()(const pat::Jet & p, pat::strbitset & ret) {
  bool val = (*this)(p);
  ret = bitmask;
  return val;
}
bool JetSelector::operator()(const pat::Jet & p) {
  bitmask.set(false);
  bool inTracking = fabs(p.eta()) < 2.4;
  if(ignoreCut("minPt")  || p.pt() > cut("minPt", double())) passCut(bitmask, "minPt");
  if(ignoreCut("maxEta") || fabs(p.eta()) < cut("maxEta", double())) passCut(bitmask, "maxEta");
  if(ignoreCut("maxNHF") || p.neutralHadronEnergyFraction() < cut("maxNHF", double())) passCut(bitmask, "maxNHF");
  if(ignoreCut("maxNEF") || p.neutralEmEnergyFraction()     < cut("maxNEF", double())) passCut(bitmask, "maxNEF");
  if(ignoreCut("minNDaughters") || (int)p.numberOfDaughters() > cut("minNDaughters", int())) passCut(bitmask, "minNDaughters");
  //Below are only used for fabs(eta) < 2.4 b/c tracking needed
  if(ignoreCut("minCHF")   || !inTracking || p.chargedHadronEnergyFraction() > cut("minCHF",   double())) passCut(bitmask, "minCHF");
  if(ignoreCut("maxCEF")   || !inTracking || p.chargedEmEnergyFraction()     < cut("maxCEF",   double())) passCut(bitmask, "maxCEF");
  if(ignoreCut("minCMult") || !inTracking || (int)p.chargedMultiplicity()    > cut("minCMult",    int())) passCut(bitmask, "minCMult");
  setIgnored(bitmask);
  return (bool) bitmask;
}
  
PhotonSelector::PhotonSelector(){}
PhotonSelector::PhotonSelector(Pset pset, std::string selectorName){
  Pset const params = pset.getParameter<Pset>(selectorName);
  // set the last parameter to false to turn off the cut
  loadFromPset<double>(params, "minPt", true);
  loadFromPset<double>(params, "maxEta", true);
  loadFromPset<double>(params, "maxECalIso", true);
  loadFromPset<double>(params, "maxHCalIso", true);
  loadFromPset<double>(params, "maxTrkIso", true);
  loadFromPset<double>(params, "maxHoE", true);
  loadFromPset<double>(params, "maxSigmaee", true);
  loadFromPset<bool>  (params, "minHasSeed", true);

  bitmask = getBitTemplate();
}
  
bool PhotonSelector::operator()(const pat::Photon & p, pat::strbitset & ret) {
  bool val = (*this)(p);
  ret = bitmask;
  return val;
}
bool PhotonSelector::operator()(const pat::Photon & p) {
  bitmask.set(false);
  if(ignoreCut("minPt")       || p.pt() > cut("minPt", double())) passCut(bitmask, "minPt");
  if(ignoreCut("maxEta")      || fabs(p.superCluster()->eta()) < cut("maxEta", double())) passCut(bitmask, "maxEta");
  if(ignoreCut("maxECalIso")  || p.ecalRecHitSumEtConeDR04() < cut("maxECalIso", double())) passCut(bitmask, "maxECalIso");
  if(ignoreCut("maxHCalIso")  || p.hcalTowerSumEtConeDR04() < cut("maxHCalIso", double())) passCut(bitmask, "maxHCalIso");
  if(ignoreCut("maxTrkIso")   || p.trkSumPtHollowConeDR04() < cut("maxTrkIso", double())) passCut(bitmask, "maxTrkIso");
  if(ignoreCut("maxHoE")      || p.hadronicOverEm() < cut("maxHoE", double())) passCut(bitmask, "maxHoE");
  if(ignoreCut("maxSigmaee")  || p.sigmaIetaIeta() < cut("maxSigmaee", double())) passCut(bitmask, "maxSigmaee");
  if(ignoreCut("minHasSeed")  || p.hasPixelSeed() >= cut("minHasSeed", bool())) passCut(bitmask, "minHasSeed");
  setIgnored(bitmask);
  return (bool) bitmask;
}
 



