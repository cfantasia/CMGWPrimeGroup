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
  case kPICKY:
    return pickyMuon();
  case kSTANDALONE:
    return standAloneMuon();
  case kCOCKTAIL:
    return (muon::tevOptimized(*this, 200, 17., 40., 0.25)).first;
  }
  std::cout<<"Failed to find a track requested ("<<algo_desc_long[muReconstructor]<<")\n";
  return globalTrack();
}

reco::TrackRef TeVMuon::muonBestTrack() const{
  return muReconstructor_ == kCOCKTAIL ? getTrack(kCOCKTAIL) : Muon::muonBestTrack();
}

void TeVMuon::printTrackerInfo() const
{
  std::cout << " trk_hits = " << innerTrack()->hitPattern().numberOfValidHits()
       << ", trk_layers = " << innerTrack()->hitPattern().trackerLayersWithMeasurement()
	 << ", trk_validFraction = " << innerTrack()->validFraction()
	    << std::endl;
    
}

void TeVMuon::printPtInfo(unsigned reconstructor) const
{
  if(!isValid(reconstructor))return;

  float rpt = 0;
  float rpt2 = 0;
  if(getTrack(reconstructor)->ptError()!=0) {
    rpt = getTrack(reconstructor)->ptError()/getTrack(reconstructor)->pt();
    rpt2 = 1000*rpt/getTrack(reconstructor)->pt();
  }
  std::cout << "\n " << algo_desc_long[reconstructor] << " pt = "
	    << getTrack(reconstructor)->pt() << " +- " 
	    << getTrack(reconstructor)->ptError()
	    << " GeV, charge = " << getTrack(reconstructor)->charge() 
	    << "\n ptError/pt = " << rpt
	    << ", ptError/pt^2 = " << rpt2 << std::endl;  
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

  // this should really be a config parameter, but this is only temporary
  const int num_layers_needed = 9;
  bool muon_hits = glb->hitPattern().trackerLayersWithMeasurement() >= num_layers_needed
    // glb->hitPattern().numberOfValidTrackerHits() > 10
    && glb->hitPattern().numberOfValidMuonHits() > 0
    && glb->hitPattern().numberOfValidPixelHits() > 0
    && numberOfMatchedStations() > 1;
  
  if(!muon_hits)
    return false;

  bool chi2_qual = true;
  if (muReconstructor_ == kGLOBAL ||
      muReconstructor_ == kINNER ||
      muReconstructor_ == kPAT ||
      muReconstructor_ == kSTANDALONE)
    chi2_qual = (glb->chi2()/glb->ndof() < chi2Cut);
  
  bool checkqual = chi2_qual
    && TMath::Abs(p3.Eta()) < muonEtaCut
    // is this d0 wrt to the beamspot???
    && TMath::Abs(dB()) < 0.02;

  if(!checkqual)return false;
  
  return true;
}


