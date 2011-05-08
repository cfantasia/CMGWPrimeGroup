#include "UserCode/CMGWPrimeGroup/interface/TeVMuon.h"
#include "UserCode/CMGWPrimeGroup/interface/util.h"

// get muon 4-d momentum according to muonReconstructor_ value
// see interface/TeVMuon.h
const TLorentzVector & TeVMuon::p4(unsigned theMu, unsigned muReconstructor,
				   bool & isInvalidMuon)
{
  switch(muReconstructor)
    {
    case 0:
      setMuLorentzVector(patMuon_->globalTrack(), isInvalidMuon);
      break;
    case 1:
      setMuLorentzVector(patMuon_->innerTrack(), isInvalidMuon);
      break;
    case 2:
      setMuLorentzVector(patMuon_->tpfmsMuon(), isInvalidMuon);
      break;
    case 3:
      setMuLorentzVector(patMuon_->cocktailMuon(), isInvalidMuon);
      break;
    case 4:
      setMuLorentzVector(patMuon_->pickyMuon(), isInvalidMuon);
      break;
    case 5:
      setMuLorentzVector(patMuon_->defaultTeVMuon(), isInvalidMuon);
      break;
    case 6:
      setMuLorentzVector(patMuon_->dytMuon(), isInvalidMuon);
      break;
    }

  return p4_;
}

void TeVMuon::setMuLorentzVector(const reco::TrackRef & trk, bool & isInvalidMuon)
{
  if(trk.isNull())
    {
      isInvalidMuon = true;
      return;
    }
  TVector3 p3(trk->px(), trk->py(), trk->pz());
  p4_.SetVectM(p3, wprime::MUON_MASS);
}

//computes the combined rel isolation value
float TeVMuon::combRelIsolation() const
{
  return ( patMuon_->ecalIso() + patMuon_->hcalIso() + 
	   patMuon_->trackIso() ) / p4_.Pt();
}

bool TeVMuon::goodQualityMuon(float chi2Cut, float muonEtaCut) const
{
  // See twiki: https://twiki.cern.ch/twiki/bin/view/CMS/ExoticaWprime
  // for the latest quality cuts
  
  bool muonID = patMuon_->isGood("AllGlobalMuons") && patMuon_->isGood("AllTrackerMuons");
  
  reco::TrackRef glb = patMuon_->globalTrack();
  if(glb.isNull())
    return false;
  
  bool muon_hits = glb->hitPattern().numberOfValidTrackerHits() > 10
    && glb->hitPattern().numberOfValidMuonHits() > 0
    && glb->hitPattern().numberOfValidPixelHits() > 0
    && patMuon_->numberOfMatches() > 1;
  
  TVector3 p3(glb->px(), glb->py(), glb->pz());
  
  bool checkqual = (glb->chi2()/glb->ndof() / chi2Cut)
    && TMath::Abs(p3.Eta()) < muonEtaCut
    // is this d0 wrt to the beamspot???
    && TMath::Abs(patMuon_->dB()) < 0.02;
  
  if(!muonID || !muon_hits || !checkqual)
    return false;
  else
    return true;

}