#ifndef _NewElec_h_
#define _NewElec_h_

#include "DataFormats/PatCandidates/interface/Electron.h"

float calcTrkIso(const pat::Electron& elec);

float calcECalIso(const pat::Electron& elec);
float calcHCalIso(const pat::Electron& elec);
float calcHCalIsoCor(const pat::Electron& elec);
float calcAdjHCalIso(const pat::Electron& elec);
float calcCombRelIso(const pat::Electron& elec, const float pu);
bool passEtCut(const pat::Electron& elec, const float cut);
bool passPtCut(const pat::Electron& elec, const float cut);
bool passWPCut(const pat::Electron& elec, const std::string type, const int mask);
bool passEtaCut(const pat::Electron& elec);
bool passNMissingHitsCut(const pat::Electron& elec, const float cut);
bool passDistDCotCut(const pat::Electron& elec, const float cut);
bool passDistCut(const pat::Electron& elec, const float cut);
bool passDeltaCotThetaCut(const pat::Electron& elec, const float cut);
bool passSigmaIEtaIEtaCut(const pat::Electron& elec, const float cut);
bool passDeltaPhiCut(const pat::Electron& elec, const float cut);
bool passDeltaEtaCut(const pat::Electron& elec, const float cut);
bool passHOverECut(const pat::Electron& elec, const float cut);
bool passCombRelIsoCut(const pat::Electron& elec, const float cut, const float pu);
#endif // #define _NewElec_h_
