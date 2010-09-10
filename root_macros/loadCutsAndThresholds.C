// wprime stuff

#include "UserCode/CMGWPrimeGroup/root_macros/loadCutsAndThresholds.h"


// Calculate efficiencies
//------------------------------------------------------------------------
void getEff(float & eff, float & deff,float Num,float Denom)
{
  //------------------------------------------------------------------------
  eff = Num/Denom;
  deff = TMath::Sqrt(eff * (1-eff)/Denom);
}//---------------getEff



// Get the hardest muon (based on tracker-pt) [iMu is the index in ev->mu array]
//------------------------------------------------------------------------
void GetTheHardestMuon(const wprime::Event * ev, int & iMu)
{
//------------------------------------------------------------------------

    int nmuons = ev->mu->GetLast() + 1;
    float temp_muPT = -999;
    for(int j = 0; j != nmuons; ++j)
      {
	wprime::Muon * mu = (wprime::Muon *) ev->mu->At(j);
	float current_muPT = mu->tracker.p.Pt();
	if (current_muPT > temp_muPT) 
	  {
	    temp_muPT = current_muPT;
	    iMu = j;
	  }
      }
}//-----GetTheHardestMuon


// returns SumPtIsolation
//-------------------------------------------------------------------
float SumPtIsolation(const wprime::Muon* the_mu, 
                     unsigned detR_iso_index)
{
//-------------------------------------------------------------------

  return the_mu->SumPtIso[detR_iso_index];

}//------ SumPtIsolation



//computes the combined rel isolation value
//-------------------------------------------------------------------
float CombRelIsolation(const wprime::Muon* the_mu,
                       unsigned detR_iso_index)
{
//-------------------------------------------------------------------    

  return (
      (the_mu->SumPtIso[detR_iso_index]+
       the_mu->ECALIso[detR_iso_index]+
       the_mu->HCALIso[detR_iso_index])/the_mu->tracker.p.Pt()
      );
  
      
}//-------- CombRelIsolation



//returns Abs(DeltaPhi) between an object(TLorentzVector)
// and the leading jet in the event
//------------------------------------------------------------------------
float XJetDPhi(const TLorentzVector& lv, const wprime::Event * ev)
{
//------------------------------------------------------------------------

 float dphi = -1;
  
  int njets = ev->jet->GetLast() + 1;
  if (njets > 0){ 
      TLorentzVector * jet = (TLorentzVector *) ev->jet->At(0);
      dphi = jet->DeltaPhi(lv);
    }
  
  
  return TMath::Abs(dphi);

}//-------------------XJetDPhi



//transverse mass with a given MET object (TVector2)
//------------------------------------------------------------------------
float TMass(const TLorentzVector& lv, const TVector2& themet)
{
//------------------------------------------------------------------------

    float tmass = 0;
    float cdphi = TMath::Cos(lv.Phi()-themet.Phi());
    float tmass_sqr = 2*lv.Pt()*themet.Mod()*(1-cdphi);
    tmass = (tmass_sqr>0) ? sqrt(tmass_sqr) : 0;
  
    return tmass;

}//-------------------TMass



// returns # of (global) muons with tracker-pt above <tracker_muon_pt>
//------------------------------------------------------------------------
unsigned NmuAboveThresh(float tracker_muon_pt, const wprime::Event * ev)
{
//------------------------------------------------------------------------
  unsigned N = 0;

  int nmuons = ev->mu->GetLast() + 1;
  for(int j = 0; j != nmuons; ++j)
    { // loop over muons
      wprime::Muon * mu = (wprime::Muon *) ev->mu->At(j);
      if(mu->tracker.p.Pt() > tracker_muon_pt)
	++N;
            
    } // loop over muons

  return N;
}//------------NumAboveThresh



// returns # of jets with Et above threshold
//------------------------------------------------------------------------
unsigned NjetAboveThresh(float threshold, const wprime::Event * ev)
{
  unsigned N = 0;
  
  int njets = ev->jet->GetLast() + 1;
  for(int j = 0; j != njets; ++j)
    { // loop over jets
      TLorentzVector * jet = (TLorentzVector *) ev->jet->At(j);
      if(jet->Et() > threshold)
	++N;
    } // loop over jets
  
  return N;
}//--------------NjetAboveThresh



// returns # of jets with Et above threshold 
// with angle greater than delta_phi from muon
//------------------------------------------------------------------------
unsigned NjetAboveThresh(float threshold, float delta_phi, 
			 const wprime::Muon * mu, const wprime::Event * ev)
{
  //------------------------------------------------------------------------
  unsigned N = 0;
  
  int njets = ev->jet->GetLast() + 1;
  for(int j = 0; j != njets; ++j)
    { // loop over jets
      TLorentzVector * jet = (TLorentzVector *) ev->jet->At(j);
      if(jet->Et() >threshold && 
	 TMath::Abs(jet->DeltaPhi(mu->tracker.p)) > delta_phi)
	++N;
    } // loop over jets
  
  return N;

}//-----------------NjetAboveThresh



// true if HLT conditions are met
//-------------------------------------------------------------------
bool PassedHLT(const wprime::Event* ev)
{
  //-------------------------------------------------------------------
#if debugmemore
  cout << " Executing PassedHLT()" << endl;
#endif
  
  // here the triggers that are to be used
  bool HLT = (ev->HLT_Mu9 == 1);
  return HLT;
}//-------PassedHLT()



// check if muon is in pt-range for the different algorithms, fill isThere
//-------------------------------------------------------------------
void CheckMuonPtInRange(const wprime::Muon* mu, bool isThere[],
			float min_MuPt, float max_MuPt)
{
//-------------------------------------------------------------------
#if debugmemore
  cout << " Processing IsMuonPtInRange() " << endl;
#endif
  
  // Make sure there is a muon
  if(mu == 0) 
    {
      for(int i = 0; i != Num_trkAlgos; ++i)
	isThere[i] = false;

      return;
    }

  isThere[0] = mu->global.p.Pt() >=min_MuPt 
    && mu->global.p.Pt() <=max_MuPt;

  isThere[1] = mu->tracker.p.Pt() >=min_MuPt 
    && mu->tracker.p.Pt() <=max_MuPt;

  isThere[2] = mu->tpfms.p.Pt() >=min_MuPt 
    && mu->tpfms.p.Pt() <=max_MuPt;

  isThere[3] = mu->cocktail.p.Pt() >=min_MuPt 
    && mu->cocktail.p.Pt() <=max_MuPt;

  isThere[4] = mu->picky.p.Pt() >=min_MuPt 
    && mu->picky.p.Pt() <=max_MuPt;

  isThere[5] = mu->tmr.p.Pt() >=min_MuPt 
    && mu->tmr.p.Pt() <=max_MuPt;

}//-------IsMuonPtInRange



// true if only one muon with track pT > the threshold
//-------------------------------------------------------------------
bool OnlyOneHighTrackPtMuon(const wprime::Event* ev, float one_mu_pt_trkcut)
{
//-------------------------------------------------------------------
#if debugmemore
 cout << " Processing OnlyOneHighTrackPtMuon() " << endl;
#endif

 return (NmuAboveThresh(one_mu_pt_trkcut, ev) == 1);
}//----------ThereIsOneMuonOnly



// true if isolation requirements satisfied for muon
//-------------------------------------------------------------------
bool Isolation(const wprime::Muon* the_mu, unsigned detR_iso_index,
               float iso_cut)
{
//-------------------------------------------------------------------
#if debugmemore
  cout << " Processing Isolation() " << endl;
#endif

  return (CombRelIsolation(the_mu,detR_iso_index) <= iso_cut);
  //return (SumPtIsolation(theMu,detR_iso_index) <= iso_cut);

}//--------Isolation




// true if energetic Jet(s) found back to back with muon 
//-------------------------------------------------------------------
bool ExceedMaxNumJetsOpposedToMu(unsigned max_jets_aboveThresh,
				 float et_jet_cut,  
				 float delta_phi,
				 const wprime::Muon* the_mu,
				 const wprime::Event* ev)
{
//-------------------------------------------------------------------
#if debugmemore
  cout << " Processing ExeedMaxNumJetsOpposedToMu() " << endl;
#endif
  return (NjetAboveThresh(et_jet_cut, delta_phi, 
			  the_mu, ev) > max_jets_aboveThresh);

}//ExceedMaxNumEnergeticJets




// check if muon satisfies quality requirements for all tracking algorithms, fill goodQual
//-------------------------------------------------------------------
void CheckQuality(const wprime::Muon* mu, 
		  bool goodQual[], float pttrk_cut,
		  float chi2_cut, float muon_etacut)
{
//-------------------------------------------------------------------
#if debugmemore
  cout << " Processing HasQuality() " << endl;
#endif

  //See twiki: https://twiki.cern.ch/twiki/bin/view/CMS/ExoticaWprime
  //for the latest quality cuts
  
  //This is consider now an old cut. Put on-hold for the moment
  bool isHard = mu->tracker.p.Pt() >= pttrk_cut;
  isHard = true;

  bool checkqual = false;
  
  //This cut is on-hold. In the barrel, it
  // is effectively the same as using
  //# hits in the tracker:
  //
  //((mu->tracker.Npixel_layer + mu->tracker.Nstrip_layer) >= 10)
  //    && ((mu->tracker.Npixel_layerNoMeas + 
  //         mu->tracker.Nstrip_layerNoMeas) < 5)
 
  bool muonID = 
      (mu->tracker.Ntrk_hits > 10)
      && (TMath::Abs(mu->tracker.d0) < 0.2) 
      && mu->AllTrackerMuons
      && mu->AllGlobalMuons  
      && mu->global.Nmuon_hits > 0
      && mu->Nmatches > 1
      && mu->tracker.Npixel_hits > 0
      ;
  
   
  checkqual = ((mu->global.chi2 / mu->global.ndof)<chi2_cut)
      && TMath::Abs(mu->global.p.Eta()) < muon_etacut;
  
  // old value: muon's pt within range
  // new value: old value .AND. quality cuts
  goodQual[0] = goodQual[0] && checkqual && isHard && muonID;


  checkqual = ((mu->tracker.chi2 / mu->tracker.ndof) 
	       < chi2_cut)
    && TMath::Abs(mu->tracker.p.Eta()) < muon_etacut;
  goodQual[1] = goodQual[1] && checkqual && isHard && muonID;



  checkqual = ((mu->tpfms.chi2 / mu->tpfms.ndof) 
	       < chi2_cut) 
    && TMath::Abs(mu->tpfms.p.Eta()) < muon_etacut;
  goodQual[2] = goodQual[2] && checkqual && isHard && muonID;

  checkqual = ((mu->cocktail.chi2 / mu->cocktail.ndof) 
	       < chi2_cut) 
    && TMath::Abs(mu->cocktail.p.Eta()) < muon_etacut;

  goodQual[3] = goodQual[3] && checkqual && isHard && muonID;


  checkqual = ((mu->picky.chi2 / mu->picky.ndof) 
	       < chi2_cut) 
    && TMath::Abs(mu->picky.p.Eta()) < muon_etacut;

  goodQual[4] = goodQual[4] && checkqual && isHard && muonID;


  checkqual = ((mu->tmr.chi2 / mu->tmr.ndof) 
	       < chi2_cut) 
    && TMath::Abs(mu->tmr.p.Eta()) < muon_etacut;

  goodQual[5] = goodQual[5] && checkqual && isHard && muonID;

}//------HasQuality




