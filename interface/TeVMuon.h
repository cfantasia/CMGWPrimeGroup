#ifndef _TeVMuon_h_
#define _TeVMuon_h_

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TLorentzVector.h"

#include "UserCode/CMGWPrimeGroup/interface/TeVMuon_tracking.h"

class TeVMuon : public pat::Muon{
 public:
  TeVMuon(const pat::Muon & muon, unsigned int muReconstructor=kTEV, bool isValid=1) :
    pat::Muon(muon), muReconstructor_(muReconstructor) {p4_ = p4(isValid_);}
  ~TeVMuon(){}
  
  // get muon 4-d momentum according to muonReconstructor_ value
  const TLorentzVector & p4(unsigned muReconstructor, 
                            bool & isInvalidMuon);

  // isInvalidMuon: set to false if muon track (for given reconstructor) is NULL
  const TLorentzVector & p4(bool & isInvalidMuon);

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

  double pt() const ;
  double pt(unsigned muReconstructor) const ;
  // computes the tracker-based relative isolation value
  float trkRelIsolation() const;
  //computes the combined rel isolation value
  float combRelIsolation() const;
  float combRelIsolation03(const float offset) const;
  
  TVector2 getPtDiff() const;
  const reco::TrackRef GetTrack(unsigned muReconstructor) const;
  const TLorentzVector getTrkLorentzVector(unsigned muReconstructor) const;
  const TLorentzVector getTrkLorentzVector(const reco::TrackRef & trk, bool & isInvalidMuon) const;

 private:
  bool isValid_;
  unsigned int muReconstructor_;
  TLorentzVector p4_;

  void setMuLorentzVector(const reco::TrackRef & trk, bool & isInvalidMuon);

};

#endif // #define _TeVMuon_h_
