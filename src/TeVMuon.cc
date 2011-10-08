#include "UserCode/CMGWPrimeGroup/interface/TeVMuon.h"
#include "UserCode/CMGWPrimeGroup/interface/util.h"

//////////////////////////////
///TeV Modifier Functions/////
//////////////////////////////

//////////////////////////////
///TeV Accessor Functions/////
//////////////////////////////

const reco::TrackRef
TeVMuon::getTrack(const unsigned muReconstructor) const{
  switch(muReconstructor)
  {// see interface/TeVMuon.h
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
  case kSTANDALONE:
    return standAloneMuon();
  }
  std::cout<<"Failed to find a track requested\n";
  return defaultTeVMuon();
}

const TLorentzVector TeVMuon::getTrkLorentzVector(const reco::TrackRef trk) const
{
  const float dummy = -9999;
  TLorentzVector trkP4;
  if(trk.isNull())
    trkP4.SetXYZM(dummy, dummy, dummy, wprime::MUON_MASS);
  else
    trkP4.SetXYZM(trk->px(), trk->py(), trk->pz(), wprime::MUON_MASS);

  return trkP4;
}

//////////////////////////////
///TeV Helper Functions///////
//////////////////////////////

//////////////////////////////
///TeV Cut Functions//////////
//////////////////////////////

bool TeVMuon::goodQualityMuon(float chi2Cut, float muonEtaCut) const
{
  // See twiki: https://twiki.cern.ch/twiki/bin/view/CMS/ExoticaWprime
  // for the latest quality cuts
  
  bool muonID = isGood("AllGlobalMuons") && isGood("AllTrackerMuons");
  
  if(!muonID) return false;

  reco::TrackRef glb = globalTrack();
  if(glb.isNull())
    return false;

  TVector3 p3(glb->px(), glb->py(), glb->pz());
  int num_layers_needed = -1;
  float eta = TMath::Abs(p3.Eta());
  if(eta > 0.9 && eta < 1.5)
     num_layers_needed = 8;
   else
     num_layers_needed = 9;
 
  bool muon_hits = glb->hitPattern().trackerLayersWithMeasurement() > num_layers_needed
    // glb->hitPattern().numberOfValidTrackerHits() > 10
    && glb->hitPattern().numberOfValidMuonHits() > 0
    && glb->hitPattern().numberOfValidPixelHits() > 0
    && numberOfMatches() > 1;
  
  if(!muon_hits)
    return false;

  bool checkqual = (glb->chi2()/glb->ndof() / chi2Cut)
    && TMath::Abs(p3.Eta()) < muonEtaCut
    // is this d0 wrt to the beamspot???
    && TMath::Abs(dB()) < 0.02;

  if(!checkqual)return false;
  
  return true;
}


