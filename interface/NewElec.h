#ifndef _NewElec_h_
#define _NewElec_h_

#include "DataFormats/PatCandidates/interface/Electron.h"

float CalcTrkIso(const pat::Electron& elec);

float CalcECalIso(const pat::Electron& elec);
float CalcHCalIso(const pat::Electron& elec);
float CalcHCalIsoCor(const pat::Electron& elec);
float CalcAdjHCalIso(const pat::Electron& elec);
float CalcCombRelIso(const pat::Electron& elec, const float pu);
bool PassEtCut(const pat::Electron& elec, const float cut);
bool PassPtCut(const pat::Electron& elec, const float cut);
bool PassWPCut(const pat::Electron& elec, const std::string type, const int mask);
bool PassEtaCut(const pat::Electron& elec);
bool PassNMissingHitsCut(const pat::Electron& elec, const float cut);
bool PassDistDCotCut(const pat::Electron& elec, const float cut);
bool PassDistCut(const pat::Electron& elec, const float cut);
bool PassDeltaCotThetaCut(const pat::Electron& elec, const float cut);
bool PassSigmaIEtaIEtaCut(const pat::Electron& elec, const float cut);
bool PassDeltaPhiCut(const pat::Electron& elec, const float cut);
bool PassDeltaEtaCut(const pat::Electron& elec, const float cut);
bool PassHOverECut(const pat::Electron& elec, const float cut);
bool PassCombRelIsoCut(const pat::Electron& elec, const float cut, const float pu);
#endif // #define _NewElec_h_
