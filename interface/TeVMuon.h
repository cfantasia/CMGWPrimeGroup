#ifndef _TeVMuon_h_
#define _TeVMuon_h_

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TLorentzVector.h"

#include "UserCode/CMGWPrimeGroup/interface/TeVMuon_tracking.h"

class TeVMuon : public pat::Muon{
 public:
  TeVMuon(const pat::Muon & muon, unsigned int algo=0, bool isValid=1) :
    pat::Muon(muon), algo_(algo) {p4_ = p4(0,algo_,isValid_);}
  ~TeVMuon(){}
  
  // get muon 4-d momentum according to muonReconstructor_ value
  const TLorentzVector & p4(unsigned theMu, unsigned muReconstructor, 
			    bool & isInvalidMuon);

  bool goodQualityMuon(float chi2Cut, float muonEtaCut) const;

  double pt() const ;
  //computes the combined rel isolation value
  float combRelIsolation() const;

 private:
  bool isValid_;
  unsigned int algo_;
  TLorentzVector p4_;

  void setMuLorentzVector(const reco::TrackRef & trk, bool & isInvalidMuon);

};

#endif // #define _TeVMuon_h_
