#include "UserCode/CMGWPrimeGroup/interface/TeVMuon.h"
#include "UserCode/CMGWPrimeGroup/interface/util.h"



// get muon 4-d momentum according to muonReconstructor_ value
// see interface/TeVMuon.h
const TLorentzVector & TeVMuon::p4(unsigned muReconstructor,
                                   bool & isInvalidMuon){
  setMuLorentzVector(GetTrack(muReconstructor), isInvalidMuon);
 
  return p4_;
}

const TLorentzVector & TeVMuon::p4(bool & isInvalidMuon){
  return p4(muReconstructor_, isInvalidMuon);
}

const reco::TrackRef
TeVMuon::GetTrack(unsigned muReconstructor) const{
  switch(muReconstructor)
  {
  case kGLOBAL:
    return globalTrack();
  case kINNER:
    return innerTrack();
  case kTPFMS:
    return tpfmsMuon();
  case kCOCKTAIL:
    return cocktailMuon();
  case kPICKY:
    return pickyMuon();
  case kTEV:
    return defaultTeVMuon();
  case kDYT:
    return dytMuon();
  }
  std::cout<<"Failed to find a track requested\n";
  return defaultTeVMuon();
}

double TeVMuon::pt() const {
  return p4_.Pt();
}

double TeVMuon::pt(unsigned muReconstructor) const {
  return getTrkLorentzVector(muReconstructor).Pt();
}

TVector2 TeVMuon::getPtDiff() const {
  TVector2  chosenAlgo( p4_.Px(), p4_.Py() );
  bool isInvalid = false;
  TLorentzVector p4Def = getTrkLorentzVector(0);//Global Muon Pt
  if(isInvalid){
    std::cout<<"Failed getting muon\n";
    return TVector2();
  }
  TVector2 defaultAlgo( p4Def.Px(), p4Def.Py() );
  return chosenAlgo - defaultAlgo;
}

void TeVMuon::setMuLorentzVector(const reco::TrackRef & trk, bool & isInvalidMuon)
{
  p4_ = getTrkLorentzVector(trk, isInvalidMuon);
}

const TLorentzVector TeVMuon::getTrkLorentzVector(unsigned muReconstructor) const{
  bool isInvalidMuon = false;
  return getTrkLorentzVector(GetTrack(muReconstructor), isInvalidMuon);
}

const TLorentzVector TeVMuon::getTrkLorentzVector(const reco::TrackRef & trk, bool & isInvalidMuon) const
{
  TLorentzVector trkP4;
  if(trk.isNull()){
    isInvalidMuon = true;
  }else{
    TVector3 p3(trk->px(), trk->py(), trk->pz());
    trkP4 = TLorentzVector(p3, wprime::MUON_MASS);
  } 
  return trkP4;
}

//computes the combined rel isolation value
float TeVMuon::combRelIsolation() const
{
  return ( ecalIso() + hcalIso() + trackIso() ) 
    / pt();
}

bool TeVMuon::goodQualityMuon(float chi2Cut, float muonEtaCut) const
{
  // See twiki: https://twiki.cern.ch/twiki/bin/view/CMS/ExoticaWprime
  // for the latest quality cuts
  
  bool muonID = isGood("AllGlobalMuons") && isGood("AllTrackerMuons");
  
  reco::TrackRef glb = globalTrack();
  if(glb.isNull())
    return false;
  
  bool muon_hits = glb->hitPattern().numberOfValidTrackerHits() > 10
    && glb->hitPattern().numberOfValidMuonHits() > 0
    && glb->hitPattern().numberOfValidPixelHits() > 0
    && numberOfMatches() > 1;
  
  TVector3 p3(glb->px(), glb->py(), glb->pz());
  
  bool checkqual = (glb->chi2()/glb->ndof() / chi2Cut)
    && TMath::Abs(p3.Eta()) < muonEtaCut
    // is this d0 wrt to the beamspot???
    && TMath::Abs(dB()) < 0.02;
  
  if(!muonID || !muon_hits || !checkqual)
    return false;
  else
    return true;

}
