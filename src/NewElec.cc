#include "UserCode/CMGWPrimeGroup/interface/NewElec.h"

float
CalcTrkIso(const pat::Electron& elec){
  return elec.dr03TkSumPt();
}

float
CalcECalIso(const pat::Electron& elec){
  return elec.isEB() ? 
    std::max(0., elec.dr03EcalRecHitSumEt() - 1.)
    : elec.dr03EcalRecHitSumEt();
}

float
CalcHCalIso(const pat::Electron& elec){
  return  elec.dr03HcalTowerSumEt();
}

float
CalcHCalIsoCor(const pat::Electron& elec){
  return elec.hadronicOverEm() * elec.energy() * fabs(sin(elec.superCluster()->position().theta()));
}

float
CalcAdjHCalIso(const pat::Electron& elec){
  return CalcHCalIso(elec) + CalcHCalIsoCor(elec);
}

float
CalcCombRelIso(const pat::Electron& elec, const float pu){
  return (CalcTrkIso(elec) + CalcECalIso(elec) + CalcAdjHCalIso(elec) - pu) /
    elec.p4().Pt();
}

inline bool PassEtCut(const pat::Electron& elec, const float cut){
  return elec.et() > cut;
}//--- PassEtCut

inline bool PassPtCut(const pat::Electron& elec, const float cut){
  return elec.pt() > cut;
}//--- PassPtCut

inline bool PassWPCut(const pat::Electron& elec, const std::string type, const int mask){
  return ((int)elec.electronID(type) & mask) == mask;
}//--- PassLooseWPCut

inline bool PassEtaCut(const pat::Electron& elec){
  return elec.isEB() || elec.isEE();
}//--- PassEtaCut

inline bool PassNMissingHitsCut(const pat::Electron& elec, const float cut){
  return elec.gsfTrack().get()->trackerExpectedHitsInner().numberOfHits() <= cut;
}

inline bool PassDistDCotCut(const pat::Electron& elec, const float cut){
  return PassDistCut(elec, cut) || PassDeltaCotThetaCut(elec, cut);
}

inline bool PassDistCut(const pat::Electron& elec, const float cut){
  return fabs(elec.convDist()) >= cut;
}

inline bool PassDeltaCotThetaCut(const pat::Electron& elec, const float cut){
  return fabs(elec.convDcot()) >= cut;
}

inline bool PassSigmaIEtaIEtaCut(const pat::Electron& elec, const float cut){
  return elec.sigmaIetaIeta() < cut;
}//--- PassSigmaEtaEtaCut

inline bool PassDeltaPhiCut(const pat::Electron& elec, const float cut){
  return fabs(elec.deltaPhiSuperClusterTrackAtVtx()) < cut;
}//--- PassDeltaPhiCut

inline bool PassDeltaEtaCut(const pat::Electron& elec, const float cut){
  return fabs(elec.deltaEtaSuperClusterTrackAtVtx()) < cut;
}//--- PassDeltaEtaCut

inline bool PassHOverECut(const pat::Electron& elec, const float cut){
  return elec.hadronicOverEm() < cut;
}//--- PassHOverECut

inline bool PassCombRelIsoCut(const pat::Electron& elec, const float cut, const float pu){
  return CalcCombRelIso(elec, pu) < cut;
}//--- PassTrkRelIsoCut


