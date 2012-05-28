#ifndef _TeVMuon_h_
#define _TeVMuon_h_

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TLorentzVector.h"

#include "UserCode/CMGWPrimeGroup/interface/TeVMuon_tracking.h"

class TeVMuon : public pat::Muon{
 public:
  TeVMuon(const reco::Muon & muon, const unsigned muReconstructor=kCOCKTAIL) :
    pat::Muon(muon), muReconstructor_(muReconstructor),patP4_(muon.px(),muon.py(),muon.pz(),muon.energy()){
    setMuLorentzVector(muReconstructor);
  }
  TeVMuon(const pat::Muon & muon, const unsigned muReconstructor=kCOCKTAIL) :
    pat::Muon(muon), muReconstructor_(muReconstructor),patP4_(muon.px(),muon.py(),muon.pz(),muon.energy()){
    
    setMuLorentzVector(muReconstructor);
  }
  ~TeVMuon(){}
  
  //////////////////////////////
  ///TeV Accessor Functions/////
  //////////////////////////////
  const unsigned getmuReconstructor() const;
  const bool isValid() const;
  const bool isValid(const unsigned muReconstructor) const;
  const TLorentzVector getPatP4() const;
  const reco::TrackRef getTrack(const unsigned muReconstructor) const;
  const TLorentzVector getTrkLorentzVector(const reco::TrackRef trk) const;

  //const TLorentzVector P4() const;
  const TLorentzVector P4(const unsigned muReconstructor) const;
  const TLorentzVector P4() const;
  double pt() const ;
  const double pt(const unsigned muReconstructor) const ;


  //////////////////////////////
  ///TeV Helper Functions///////
  //////////////////////////////

  // computes the tracker-based relative isolation value
  float trkRelIsolation() const;

  //computes the combined rel isolation value
  float combRelIsolation() const;
  float combRelIsolation03(const float offset) const;
  float combRelPFIsolation() const;
  
//////////////////////////////
///TeV Cut Functions//////////
//////////////////////////////

  bool goodQualityMuon(float chi2Cut, float muonEtaCut) const;
  bool passEtaCut(const float cut);
  bool passPtCut(const float cut);
  bool passIsGlobalCut();
  bool passIsTrackerCut();
  bool passDxyCut(const float cut);
  bool passNPixHitCut(const float cut);
  bool passNTrkHitCut(const float cut);
  bool passNormChi2Cut(const float cut);
  bool passHitsUsedCut(const float cut);
  bool passStationsCut(const float cut);
  bool passCombRelIsoCut(const float cut);
  bool passCombRelIso03Cut(const float cut, const float puOffset=0.);
  
  void printTrackerInfo() const;
  void printPtInfo(unsigned reconstructor) const;

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
TeVMuon::isValid(const unsigned muReconstructor) const{
  return !(getTrack(muReconstructor).isNull());
}

inline const TLorentzVector TeVMuon::getPatP4() const{
  return patP4_;
}

/*
inline const TLorentzVector TeVMuon::P4() const{
  return TLorentzVector(px(), py(), pz(), energy());
}
*/
inline const TLorentzVector TeVMuon::P4(const unsigned muReconstructor) const{
  return muReconstructor==kPAT ? getPatP4() :
    getTrkLorentzVector(getTrack(muReconstructor));
}

inline const TLorentzVector TeVMuon::P4() const{
  return P4(muReconstructor_);
}

inline double TeVMuon::pt() const{
  return pat::Muon::pt();
}

inline const double TeVMuon::pt(const unsigned muReconstructor) const {
  return P4(muReconstructor).Pt();
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

//computes the combined rel isolation value
inline float TeVMuon::combRelPFIsolation() const
{
  return ( chargedHadronIso() + neutralHadronIso() + photonIso() ) 
    / pt();
}

//////////////////////////////
///TeV Cut Functions//////////
//////////////////////////////

inline bool TeVMuon::passIsGlobalCut(){
  return isGlobalMuon(); 
}//--- passIsGlobalCut

inline bool TeVMuon::passIsTrackerCut(){
  return isTrackerMuon(); 
}//--- passIsTrackerCut

inline bool TeVMuon::passPtCut(const float cut){
  return pt() > cut;
}

inline bool TeVMuon::passEtaCut(const float cut){
  return fabs(eta()) < cut;
}//--- passEta Cut

inline bool TeVMuon::passDxyCut(const float cut){
  return (fabs(dB()) < cut);
}//--- passDxyCut

inline bool TeVMuon::passNPixHitCut(const float cut){
  return (globalTrack()->hitPattern().numberOfValidPixelHits() > cut);
}//--- passNpixhitCut

inline bool TeVMuon::passNTrkHitCut(const float cut){
  return (globalTrack()->hitPattern().numberOfValidTrackerHits() > cut);
}//--- passNtrkhitCut

inline bool TeVMuon::passNormChi2Cut(const float cut){
  return (globalTrack()->normalizedChi2() < cut);
}//--- passChi2Cut

inline bool TeVMuon::passHitsUsedCut(const float cut){
  return (globalTrack()->hitPattern().numberOfValidMuonHits() > cut);
}//--- passHits Used Cut

inline bool TeVMuon::passStationsCut(const float cut){
  return numberOfMatches() > cut;
}//--- passStationsCut

inline bool TeVMuon::passCombRelIso03Cut(const float cut, const float puOffset){
  return combRelIsolation03(puOffset) < cut;
}//--- passCombRelIsoCut

//////////////////////////////
///TeV Modifier Functions/////
//////////////////////////////

inline void TeVMuon::setMuLorentzVector(const TLorentzVector& p4)
{
  setP4(LorentzVector(p4.Px(), p4.Py(), p4.Pz(), p4.E()));
}

inline void TeVMuon::setMuLorentzVector(const reco::TrackRef & trk)
{
  setisValid(!trk.isNull());
  setMuLorentzVector(getTrkLorentzVector(trk));
}

inline void TeVMuon::setMuLorentzVector(const unsigned muReconstructor)
{
  if(muReconstructor==kPAT)
    {
      setisValid(true);
      return;
    }
  setMuLorentzVector(getTrack(muReconstructor));
}

inline void TeVMuon::setisValid(const bool val){
  isValid_ = val;
}

#endif // #define _TeVMuon_h_
