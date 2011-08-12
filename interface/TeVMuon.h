#ifndef _TeVMuon_h_
#define _TeVMuon_h_

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TLorentzVector.h"

#include "UserCode/CMGWPrimeGroup/interface/TeVMuon_tracking.h"

class TeVMuon : public pat::Muon{
 public:
  TeVMuon(const reco::Muon & muon, const unsigned muReconstructor=kTEV, bool isValid=1) :
    pat::Muon(muon), muReconstructor_(muReconstructor),isValid_(isValid),patP4_(muon.px(),muon.py(),muon.pz(),muon.energy()){
    setMuLorentzVector(muReconstructor);
  }
  TeVMuon(const pat::Muon & muon, const unsigned muReconstructor=kTEV, bool isValid=1) :
    pat::Muon(muon), muReconstructor_(muReconstructor),isValid_(isValid),patP4_(muon.px(),muon.py(),muon.pz(),muon.energy()){
    
    setMuLorentzVector(muReconstructor);
  }
  ~TeVMuon(){}
  
  //////////////////////////////
  ///TeV Accessor Functions/////
  //////////////////////////////
  const unsigned getmuReconstructor() const;
  const bool isValid() const;
  const bool isValid(const unsigned muReconstructor);
  const TLorentzVector getPatP4() const;
  const reco::TrackRef getTrack(const unsigned muReconstructor) const;
  const TLorentzVector getTrkLorentzVector(const reco::TrackRef trk) const;

  //const TLorentzVector P4() const;
  const TLorentzVector p4(const unsigned muReconstructor) const;
  double pt() const ;
  const double pt(const unsigned muReconstructor) const ;

  //////////////////////////////
  ///TeV Helper Functions///////
  //////////////////////////////

  TVector2 getPtDiff() const;

  // computes the tracker-based relative isolation value
  float trkRelIsolation() const;

  //computes the combined rel isolation value
  float combRelIsolation() const;
  float combRelIsolation03(const float offset) const;
  
//////////////////////////////
///TeV Cut Functions//////////
//////////////////////////////

  bool goodQualityMuon(float chi2Cut, float muonEtaCut) const;
  bool PassEtaCut(const float cut);
  bool PassPtCut(const float cut);
  bool PassIsGlobalCut();
  bool PassIsTrackerCut();
  bool PassDxyCut(const float cut);
  bool PassNPixHitCut(const float cut);
  bool PassNTrkHitCut(const float cut);
  bool PassNormChi2Cut(const float cut);
  bool PassHitsUsedCut(const float cut);
  bool PassStationsCut(const float cut);
  bool PassCombRelIsoCut(const float cut);
  bool PassCombRelIso03Cut(const float cut, const float puOffset=0.);
  
 private:
  unsigned muReconstructor_;
  bool isValid_;
  TLorentzVector patP4_;

  //////////////////////////////
  ///TeV Modifier Functions/////
  //////////////////////////////
  void setMuLorentzVector(const TLorentzVector& p4);
  void setMuLorentzVector(const reco::TrackRef & trk);
  void setMuLorentzVector(const unsigned muReconstructor);
  void setisValid(const bool val);
};

//////////////////////////////
///TeV Accessor Functions/////
//////////////////////////////

inline const unsigned TeVMuon::getmuReconstructor() const {
  return muReconstructor_;
}

inline const bool TeVMuon::isValid() const {
  return isValid_;
}

inline const bool
TeVMuon::isValid(const unsigned muReconstructor){
  isValid_ = !(getTrack(muReconstructor).isNull());
  return isValid_;
}

inline const TLorentzVector TeVMuon::getPatP4() const{
  return patP4_;
}

/*
inline const TLorentzVector TeVMuon::P4() const{
  return TLorentzVector(px(), py(), pz(), energy());
}
*/
inline const TLorentzVector TeVMuon::p4(const unsigned muReconstructor) const{
  return muReconstructor==kPAT ? getPatP4() :
    getTrkLorentzVector(getTrack(muReconstructor));
}

inline double TeVMuon::pt() const{
  return pat::Muon::pt();
}

inline const double TeVMuon::pt(const unsigned muReconstructor) const {
  return p4(muReconstructor).Pt();
}

//////////////////////////////
///TeV Helper Functions///////
//////////////////////////////

// computes the tracker-based relative isolation value
inline float TeVMuon::trkRelIsolation() const
{
  return trackIso()/pt();
}

//computes the combined rel isolation value
inline float TeVMuon::combRelIsolation() const
{
  return ( ecalIso() + hcalIso() + trackIso() ) 
    / pt();
}

inline float TeVMuon::combRelIsolation03(const float offset) const{
  return (isolationR03().emEt + isolationR03().hadEt + isolationR03().sumPt - offset)
    / pt();
}

//////////////////////////////
///TeV Cut Functions//////////
//////////////////////////////

inline bool TeVMuon::PassIsGlobalCut(){
  return isGlobalMuon(); 
}//--- PassIsGlobalCut

inline bool TeVMuon::PassIsTrackerCut(){
  return isTrackerMuon(); 
}//--- PassIsTrackerCut

inline bool TeVMuon::PassPtCut(const float cut){
  return pt() > cut;
}

inline bool TeVMuon::PassEtaCut(const float cut){
  return fabs(eta()) < cut;
}//--- PassEta Cut

inline bool TeVMuon::PassDxyCut(const float cut){
  return (fabs(dB()) < cut);
}//--- PassDxyCut

inline bool TeVMuon::PassNPixHitCut(const float cut){
  return (globalTrack()->hitPattern().numberOfValidPixelHits() > cut);
}//--- PassNpixhitCut

inline bool TeVMuon::PassNTrkHitCut(const float cut){
  return (globalTrack()->hitPattern().numberOfValidTrackerHits() > cut);
}//--- PassNtrkhitCut

inline bool TeVMuon::PassNormChi2Cut(const float cut){
  return (globalTrack()->normalizedChi2() < cut);
}//--- PassChi2Cut

inline bool TeVMuon::PassHitsUsedCut(const float cut){
  return (globalTrack()->hitPattern().numberOfValidMuonHits() > cut);
}//--- PassHits Used Cut

inline bool TeVMuon::PassStationsCut(const float cut){
  return numberOfMatches() > cut;
}//--- PassStationsCut

inline bool TeVMuon::PassCombRelIso03Cut(const float cut, const float puOffset){
  return combRelIsolation03(puOffset) < cut;
}//--- PassCombRelIsoCut

//////////////////////////////
///TeV Modifier Functions/////
//////////////////////////////

inline void TeVMuon::setMuLorentzVector(const TLorentzVector& p4)
{
  setP4(LorentzVector(p4.Px(), p4.Py(), p4.Pz(), p4.E()));
}

inline void TeVMuon::setMuLorentzVector(const reco::TrackRef & trk)
{
  setMuLorentzVector(getTrkLorentzVector(trk));
}

inline void TeVMuon::setMuLorentzVector(const unsigned muReconstructor)
{
  if(muReconstructor!=kPAT && isValid(muReconstructor))
    setMuLorentzVector(getTrack(muReconstructor));
}

inline void TeVMuon::setisValid(const bool val){
  isValid_ = val;
}

#endif // #define _TeVMuon_h_
