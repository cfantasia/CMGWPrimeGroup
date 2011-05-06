#ifndef _TeVMuon_h_
#define _TeVMuon_h_

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TLorentzVector.h"

#include "UserCode/CMGWPrimeGroup/interface/TeVMuon_tracking.h"

class TeVMuon

{
 public:
  TeVMuon(const pat::Muon & muon) :
    patMuon_(&muon) {}
  ~TeVMuon(){}
  
  const pat::Muon * patMuon() const {return patMuon_;}
  
  reco::TrackRef innerTrack() const {return patMuon_->innerTrack();}

  // get muon 4-d momentum according to muonReconstructor_ value
  const TLorentzVector & p4(unsigned theMu, unsigned muReconstructor, 
			    bool & isInvalidMuon);

  bool goodQualityMuon(float chi2Cut, float muonEtaCut) const;


  //computes the combined rel isolation value
  float combRelIsolation() const;

 private:
  const pat::Muon* patMuon_;
  TLorentzVector p4_;

  void setMuLorentzVector(const reco::TrackRef & trk, bool & isInvalidMuon);

};

#endif // #define _TeVMuon_h_
