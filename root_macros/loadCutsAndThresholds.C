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
bool PassedHLT(const wprime::Event* ev, const wprime::Muon*, bool [])
{
  //-------------------------------------------------------------------
#if debugmemore
  cout << " Executing PassedHLT()" << endl;
#endif
  
  // here the triggers that are to be used
  bool HLT = (ev->HLT_Mu9 == 1) || (ev->HLT_Mu11 == 1);
  return HLT;
}//-------PassedHLT()



// check if muon is in pt-range for the different algorithms, fill isThere
// always returns true if muon != NULL
//-------------------------------------------------------------------
bool MuonPtWithinRange(const wprime::Event*, const wprime::Muon* mu, bool isThere[])
{
//-------------------------------------------------------------------
#if debugmemore
  cout << " Processing MuonPtWithinRange() " << endl;
#endif
  
  // Make sure there is a muon
  if(mu == 0) 
    {
      for(int i = 0; i != Num_trkAlgos; ++i)
	isThere[i] = false;

      return false;
    }

  const TLorentzVector * P[Num_trkAlgos] = {
    &(mu->global.p), &(mu->tracker.p), &(mu->tpfms.p), &(mu->cocktail.p),
    &(mu->picky.p),  &(mu->tmr.p)};

  for(int algo = 0; algo != Num_trkAlgos; ++algo)
    {
      float pt = (P[algo])->Pt();
      isThere[algo] = (pt >= minPtMu);
    }

  return true;

}//-------MuonPtWithinRange



// true if only one muon with track pT > the threshold
//-------------------------------------------------------------------
bool OnlyOneHighTrackPtMuon(const wprime::Event* ev, const wprime::Muon*, bool [])
{
//-------------------------------------------------------------------
#if debugmemore
 cout << " Processing OnlyOneHighTrackPtMuon() " << endl;
#endif

 return (NmuAboveThresh(OneMuPtTrackCut, ev) == 1);
}//----------ThereIsOneMuonOnly



// true if isolation requirements satisfied for muon
//-------------------------------------------------------------------
bool IsolatedMuon(const wprime::Event*, const wprime::Muon* the_mu, bool [])
{
//-------------------------------------------------------------------
#if debugmemore
  cout << " Processing IsolatedMuon() " << endl;
#endif

  return (CombRelIsolation(the_mu,deltaRIsoIndex) <= CombRelCut);
  //return (SumPtIsolation(theMu,deltaRIsoIndex) <= iso_cut);

}//--------IsolatedMuon

// true if there is no significant jet activity in event
// (wrapper for ExceedMaxNumJetsOpposedToMu)
bool NoJetActivity(const wprime::Event* ev, const wprime::Muon* the_mu, bool [])
//-------------------------------------------------------------------
{
  //  if(TMass(the_mu->global.p, ev->pfmet) < 200)return false;
#if 0
  float ratio = the_mu->tracker.p.Pt()/(ev->pfmet.Mod());
  if(ratio < 0.5 || ratio > 1.5) return false;
  TVector2 muon_T(the_mu->tracker.p.Px(), the_mu->tracker.p.Py());
  float delta_phi = ev->pfmet.DeltaPhi(muon_T);
  if(TMath::Abs(delta_phi) < 2.5)return false;
  else return true;
#endif
  return !ExceedMaxNumJetsOpposedToMu(MaxNjetsAboveThresh, EtJetCut, 
				      Delta_Phi, the_mu,ev);

} // --------------NoJetActivity

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




// check if muon satisfies quality requirements for all tracking algorithms.
// fill goodQual; always return true
//-------------------------------------------------------------------
bool GoodQualityMuon(const wprime::Event*, const wprime::Muon* mu, bool goodQual[])
{
//-------------------------------------------------------------------
#if debugmemore
  cout << " Processing GoodQualityMuon() " << endl;
#endif

  //See twiki: https://twiki.cern.ch/twiki/bin/view/CMS/ExoticaWprime
  //for the latest quality cuts
  
  //This is consider now an old cut. Put on-hold for the moment
  bool isHard = mu->tracker.p.Pt() >= PtTrackCut;
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
  
  const wprime::Track * tk[Num_trkAlgos] = {
    &(mu->global), &(mu->tracker), &(mu->tpfms), &(mu->cocktail),
    &(mu->picky),  &(mu->tmr)};

  for(int algo = 0; algo != Num_trkAlgos; ++algo)
    {
      checkqual = (( (tk[algo])->chi2 / (tk[algo])->ndof) < Chi2Cut)
	&& TMath::Abs( (tk[algo])->p.Eta()) < Muon_Eta_Cut;

      // old value: muon's pt within range
      // new value: old value .AND. quality cuts
      goodQual[algo] = goodQual[algo] && checkqual && isHard && muonID;
    }
  return true;

}//------GoodQualityMuon


// determine before event-loop the order in which cuts are to be executed
void setupCutOrder(selection_map & cuts)
{
#if debugme
  cout << "\n Cuts will be applied in this order: " << endl;
#endif

  for(int cut_i = 0; cut_i != Num_selection_cuts; ++cut_i)
    { // loop over selection cuts
      string arg = cuts_desc_short[cut_i];
#if debugme
      cout << " Cut #" << (cut_i+1) << ": " << cuts_desc_long[cut_i]
	   << " (" << arg << ") " << endl;
#endif      
      if(arg == "hlt")cuts[arg] = &PassedHLT;
      else if(arg == "ptrange")cuts[arg] = &MuonPtWithinRange;
      else if(arg == "qual")cuts[arg] = &GoodQualityMuon;
      else if(arg == "1mu")cuts[arg] = &OnlyOneHighTrackPtMuon;
      else if(arg == "iso")cuts[arg] = &IsolatedMuon;
      else if(arg == "jet")cuts[arg] = &NoJetActivity;
      else
	{
	  cout << " Oops! Don't understand how to prepare for cut nicknamed as "
	       << arg << endl;
	  abort();
	}
    } // loop over selection cuts

  cout << endl;

}

