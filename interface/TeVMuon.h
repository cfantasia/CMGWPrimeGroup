#ifndef _TeVMuon_h_
#define _TeVMuon_h_

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TLorentzVector.h"

#include "UserCode/CMGWPrimeGroup/interface/TeVMuon_tracking.h"

class TeVMuon : public pat::Muon{
 public:
  TeVMuon(const pat::Muon & muon, unsigned int muReconstructor=0, bool isValid=1) :
    pat::Muon(muon), muReconstructor_(muReconstructor) {p4_ = p4(isValid_);}
  ~TeVMuon(){}
  
  // get muon 4-d momentum according to muonReconstructor_ value;
  // isInvalidMuon: set to false if muon track (for given reconstructor) is NULL
  const TLorentzVector & p4(bool & isInvalidMuon);

  bool goodQualityMuon(float chi2Cut, float muonEtaCut) const;

  double pt() const ;
  //computes the combined rel isolation value
  float combRelIsolation() const;

 private:
  bool isValid_;
  unsigned int muReconstructor_;
  TLorentzVector p4_;

  void setMuLorentzVector(const reco::TrackRef & trk, bool & isInvalidMuon);

};

#endif // #define _TeVMuon_h_
