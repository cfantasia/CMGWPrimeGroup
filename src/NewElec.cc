#include "UserCode/CMGWPrimeGroup/interface/NewElec.h"

float
calcTrkIso(const pat::Electron& elec){
  return elec.dr03TkSumPt();
}

float
calcECalIso(const pat::Electron& elec){
  return elec.isEB() ? 
    std::max(0., elec.dr03EcalRecHitSumEt() - 1.)
    : elec.dr03EcalRecHitSumEt();
}

float
calcHCalIso(const pat::Electron& elec){
  return  elec.dr03HcalTowerSumEt();
}

float
calcHCalIsoCor(const pat::Electron& elec){
  return elec.hadronicOverEm() * elec.energy() * fabs(sin(elec.superCluster()->position().theta()));
}

float
calcAdjHCalIso(const pat::Electron& elec){
  return calcHCalIso(elec) + calcHCalIsoCor(elec);
}

float
calcCombRelIso(const pat::Electron& elec, const float pu){
  return (calcTrkIso(elec) + calcECalIso(elec) + calcAdjHCalIso(elec) - pu) /
    elec.p4().Pt();
}

inline bool passEtCut(const pat::Electron& elec, const float cut){
  return elec.et() > cut;
}//--- passEtCut

inline bool passPtCut(const pat::Electron& elec, const float cut){
  return elec.pt() > cut;
}//--- passPtCut

inline bool passWPCut(const pat::Electron& elec, const std::string type, const int mask){
  return ((int)elec.electronID(type) & mask) == mask;
}//--- passLooseWPCut

inline bool passEtaCut(const pat::Electron& elec){
  return elec.isEB() || elec.isEE();
}//--- passEtaCut

inline bool passNMissingHitsCut(const pat::Electron& elec, const float cut){
  return elec.gsfTrack().get()->trackerExpectedHitsInner().numberOfHits() <= cut;
}

inline bool passDistDCotCut(const pat::Electron& elec, const float cut){
  return passDistCut(elec, cut) || passDeltaCotThetaCut(elec, cut);
}

inline bool passDistCut(const pat::Electron& elec, const float cut){
  return fabs(elec.convDist()) >= cut;
}

inline bool passDeltaCotThetaCut(const pat::Electron& elec, const float cut){
  return fabs(elec.convDcot()) >= cut;
}

inline bool passSigmaIEtaIEtaCut(const pat::Electron& elec, const float cut){
  return elec.sigmaIetaIeta() < cut;
}//--- passSigmaEtaEtaCut

inline bool passDeltaPhiCut(const pat::Electron& elec, const float cut){
  return fabs(elec.deltaPhiSuperClusterTrackAtVtx()) < cut;
}//--- passDeltaPhiCut

inline bool passDeltaEtaCut(const pat::Electron& elec, const float cut){
  return fabs(elec.deltaEtaSuperClusterTrackAtVtx()) < cut;
}//--- passDeltaEtaCut

inline bool passHOverECut(const pat::Electron& elec, const float cut){
  return elec.hadronicOverEm() < cut;
}//--- passHOverECut

inline bool passCombRelIsoCut(const pat::Electron& elec, const float cut, const float pu){
  return calcCombRelIso(elec, pu) < cut;
}//--- passTrkRelIsoCut


