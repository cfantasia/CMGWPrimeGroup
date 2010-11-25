// wprime stuff

#include "UserCode/CMGWPrimeGroup/root_macros/loadCutsAndThresholds.h"

vector<string> exactTriggers;
vector<string> versionedTriggers;

int PtMtCutIndex = -1;

bool applyCorrection = false;
bool hadronicMETcalculated = false;

TVector2 hadronicMETcached;

TH1D * hRecoilPerp = 0;
TH2D * hRecoilParalvsVBPt = 0;
TH1D ** histRecoilParal = NULL;

void setRecoilPerp(TH1D * h){hRecoilPerp = h;}
void setRecoilParalvsVBPt(TH2D * hh){hRecoilParalvsVBPt = hh;}
void setApplyCorrection(bool flag){applyCorrection = flag;}
void setHadronicMETCalculated(bool flag){hadronicMETcalculated = flag;}

void setRecoilProjections()
{
  assert(hRecoilParalvsVBPt != 0);
  assert(histRecoilParal == 0);
  int N = hRecoilParalvsVBPt->GetXaxis()->GetNbins();
  histRecoilParal = new TH1D *[N];

  for(int bin_no = 1; bin_no <= N; ++bin_no)
    {
      // Get projection in the W pt bin
      histRecoilParal[bin_no-1] = new TH1D(*(hRecoilParalvsVBPt->ProjectionY("_pbinWpt_paral", bin_no, bin_no)));

      if(histRecoilParal[bin_no-1]->Integral()==0) {
	cout << " *** Couldn't correct for recoil for bin_no = " << bin_no
	     << endl;

      }
    }
}

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


// get hadronic MET component (that needs to be corrected 
// if applyCorrection=true)from Z data; this will be done according to hadronic 
// activity from Z->mumu reconstructed events
TVector2 getHadronicMET(const wprime::Event * ev)
{
  if (!applyCorrection)
    return ev->pfmetaddmu;

  if(hadronicMETcalculated)
    return hadronicMETcached;

  int nW = ev->w_mc->GetLast() + 1;
  TLorentzVector * W_p4 = 0;
  for(int j = 0; j != nW; ++j)
    {
      wprime::MCParticle * w = (wprime::MCParticle *) ev->w_mc->At(j);
      if(w->status != 3)continue;
      W_p4 = &(w->p);
    }

  // correct hadronic MET by taking into account recoil for given W pt

  // Find bin corresponding to W pt
  int binWpt = (hRecoilParalvsVBPt->GetXaxis())->FindBin(W_p4->Pt());
  // protect against outliers: EITHER use the last bin OR return MET w/o correction
  if(binWpt > hRecoilParalvsVBPt->GetXaxis()->GetNbins())
    binWpt = hRecoilParalvsVBPt->GetXaxis()->GetNbins();


#if 0
  double pt_mean_paral = histRecoilParal[binWpt-1]->GetMean();
  double pt_sigma_paral = histRecoilParal[binWpt-1]->GetRMS();
  double pt_mean_perp = hRecoilPerp->GetMean();
  double pt_sigma_perp = hRecoilPerp->GetRMS();
#endif

  // Shoot random MET (parallel) from the Z histograms
  double dataSampledMETParal= histRecoilParal[binWpt-1]->GetRandom();
  // Shoot perpendicular component
  double dataSampledMETPerp = hRecoilPerp->GetRandom();

  // Rotate back from system in which the boson is in the x axis
  double cosW = W_p4->Px()/W_p4->Pt();
  double sinW = W_p4->Py()/W_p4->Pt();
  double dataSampledMEx = 
    cosW*dataSampledMETParal - sinW*dataSampledMETPerp;
  double dataSampledMEy = 
    sinW*dataSampledMETParal + cosW*dataSampledMETPerp;

  hadronicMETcached.Set(dataSampledMEx, dataSampledMEy);
  setHadronicMETCalculated(true);
  return hadronicMETcached;    
}


// Get new MET: there are two corrections to be made:
// (a) the hadronic MET component (that needs to be corrected 
// if applyCorrection=true)from Z data; this will be done according to hadronic 
// activity from Z->mumu reconstructed events
// (b) the muon-pt component that needs to be updated if we switch to one
// of the dedicated high-pt muon reconstructors
TVector2 getNewMET(const wprime::Event * ev, const TLorentzVector & mu_p) 
{    

  return getHadronicMET(ev) - TVector2(mu_p.Px(), mu_p.Py());
}


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
bool PassedHLT(const wprime::Event* ev, const wprime::Muon*, bool [], bool [])
{
  //-------------------------------------------------------------------
#if debugmemore
  cout << " Executing PassedHLT()" << endl;
#endif

  if(exactTriggers.empty())
    {
      char exact[1024]; char versioned[1024];
      // consider all HLT_Mu triggers up to 25 GeV
      for(int thr = 9; thr <= 25; thr += 2) 
	{
	  sprintf(exact, "HLT_Mu%d", thr);
	  sprintf(versioned, "HLT_Mu%d_v", thr);
	  exactTriggers.push_back(exact);
	  versionedTriggers.push_back(versioned);
	}
    }

  int nhlt = ev->hlt->GetLast() + 1;

  vector<string>::const_iterator it = exactTriggers.begin();
  vector<string>::const_iterator it2 = versionedTriggers.begin();
  
  while(it != exactTriggers.end() && it2 != versionedTriggers.end())
    { // loop over desired triggers

      for(int j = 0; j!= nhlt;++j){ // loop over triggers stored in event
	wprime::TrigInfo * trig = (wprime::TrigInfo *) ev->hlt->At(j);
	string name = trig->name.c_str();
	
	if(trig->fired != 1 || trig->hltpre != 1)continue;
	
	if((name == *it) || (name.find(*it2) != string::npos))
	  return true;
    
      } //loop over triggers stored in event

      ++it; ++it2;

    } // loop over desired triggers

  return false;


}//-------PassedHLT()



// check if muon is in pt/Mt-range for the different algorithms, fill isThere
// always returns true if muon != NULL
//-------------------------------------------------------------------
bool MuonPtMtWithinRange(const wprime::Event* ev, const wprime::Muon* mu, 
		       bool isTherePt[], bool isThereMt[])
{
//-------------------------------------------------------------------
#if debugmemore
  cout << " Processing MuonPtMtWithinRange() " << endl;
#endif

  ++PtMtCutIndex; 
  if(PtMtCutIndex == NumPtMtThresholds)PtMtCutIndex = 0;

  // Make sure there is a muon
  if(mu == 0) 
    {
      for (int mual = MuAlgo_MIN; mual <= MuAlgo_MAX; ++mual)
	{
	  isTherePt[mual] = isThereMt[mual] = false;
	}
      return false;
    }
  
  const TLorentzVector * P[Num_trkAlgos] = {
    &(mu->global.p), &(mu->tracker.p), &(mu->tpfms.p), &(mu->cocktail.p),
    &(mu->picky.p),  &(mu->tmr.p), &(mu->dyt.p)};
  
  for (int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
    {
      TVector2 MET = getNewMET(ev, *(P[algo]));
      
      float pt = (P[algo])->Pt();
      float Mt = TMass(*P[algo], MET);
      
      isTherePt[algo] = isTherePt[algo] &&(pt>=PtThreshold[PtMtCutIndex]);
      isThereMt[algo] = isThereMt[algo] &&(Mt>=MtThreshold[PtMtCutIndex]);
    }

  return true;
  
}//-------MuonPtWithinRange


// true if only one muon with track pT > the threshold
//-------------------------------------------------------------------
bool OnlyOneHighTrackPtMuon(const wprime::Event* ev, const wprime::Muon*, 
			    bool [], bool [])
{
//-------------------------------------------------------------------
#if debugmemore
 cout << " Processing OnlyOneHighTrackPtMuon() " << endl;
#endif

 return (NmuAboveThresh(OneMuPtTrackCut, ev) == 1);
}//----------ThereIsOneMuonOnly



// true if isolation requirements satisfied for muon
//-------------------------------------------------------------------
bool IsolatedMuon(const wprime::Event*, const wprime::Muon* the_mu, 
		  bool [], bool [])
{
//-------------------------------------------------------------------
#if debugmemore
  cout << " Processing IsolatedMuon() " << endl;
#endif

  return (CombRelIsolation(the_mu,deltaRIsoIndex) <= CombRelCut);
  //return (SumPtIsolation(theMu,deltaRIsoIndex) <= iso_cut);

}//--------IsolatedMuon

// check if muon, MET pass kinematic cuts, update goodQualMt (always return true)
bool KinematicCuts(const wprime::Event* ev, const wprime::Muon* mu, 
		   bool goodQualMt[])
{
#if debugmemore
  cout << " Processing KinematicCuts() " << endl;
#endif

  const TLorentzVector * P[Num_trkAlgos] = {
    &(mu->global.p), &(mu->tracker.p), &(mu->tpfms.p), &(mu->cocktail.p),
    &(mu->picky.p),  &(mu->tmr.p), &(mu->dyt.p)};
  
  for (int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
    {
      if(!goodQualMt[algo])continue; // this is only for CPU-savings

      TVector2 MET = getNewMET(ev, *(P[algo]));
      float ratio = P[algo]->Pt()/MET.Mod();

      TVector2 muon_T(P[algo]->Px(), P[algo]->Py());
      float delta_phi = MET.DeltaPhi(muon_T);

      if(ratio < 0.4 || ratio > 1.5 || TMath::Abs(delta_phi) < 2.5)
	goodQualMt[algo] = false;

    }

  return true;      
}

// call NoJetActivity for muon-pt analysis and KinematicCuts for Mt analysis
bool NoJetActivityKinematicCuts(const wprime::Event* ev, const wprime::Muon* mu,
				bool goodQualPt[], bool goodQualMt[])
{
  if(!NoJetActivity(ev, mu))
    {
      for (int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
	goodQualPt[algo] = false;
    }
  KinematicCuts(ev, mu, goodQualMt);

  return true;
}


// true if there is no significant jet activity in event
// (wrapper for ExceedMaxNumJetsOpposedToMu)
bool NoJetActivity(const wprime::Event* ev, const wprime::Muon* the_mu)
//-------------------------------------------------------------------
{
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
bool GoodQualityMuon(const wprime::Event*, const wprime::Muon* mu, 
		     bool goodQualPt[], bool goodQualMt[])
{
//-------------------------------------------------------------------
#if debugmemore
  cout << " Processing GoodQualityMuon() " << endl;
#endif

  //See twiki: https://twiki.cern.ch/twiki/bin/view/CMS/ExoticaWprime
  //for the latest quality cuts
  
  //This is consider now an old cut. Put on-hold for the moment
  //  bool isHard = mu->tracker.p.Pt() >= PtTrackCut;

  bool checkqual = false;
  
  bool muonID = 
    (mu->global.Ntrk_hits > 10)
    && (TMath::Abs(mu->global.d0) < 0.02) 
    && mu->AllTrackerMuons
    && mu->AllGlobalMuons  
    && mu->global.Nmuon_hits > 0
    && mu->Nmatches > 1
    && mu->global.Npixel_hits > 0
    && ((mu->global.chi2)/(mu->global.ndof) < Chi2Cut)
      ;
  
  const wprime::Track * tk[Num_trkAlgos] = {
    &(mu->global), &(mu->tracker), &(mu->tpfms), &(mu->cocktail),
    &(mu->picky),  &(mu->tmr), &(mu->dyt)};

  for (int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
    {
      checkqual = //(( (tk[algo])->chi2 / (tk[algo])->ndof) < Chi2Cut) &&
	TMath::Abs( (tk[algo])->p.Eta()) < Muon_Eta_Cut;
      
      // old value: muon's pt within range
      // new value: old value .AND. quality cuts
      goodQualPt[algo] = goodQualPt[algo]&& checkqual && muonID;
      goodQualMt[algo] = goodQualMt[algo]&& checkqual && muonID;
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
      else if(arg.find("thr") != string::npos)
	cuts[arg] = &MuonPtMtWithinRange;
      else if(arg == "qual")cuts[arg] = &GoodQualityMuon;
      else if(arg == "1mu")cuts[arg] = &OnlyOneHighTrackPtMuon;
      else if(arg == "iso")cuts[arg] = &IsolatedMuon;
      else if(arg == "jetmet")cuts[arg] = &NoJetActivityKinematicCuts;
      else
	{
	  cout << " Oops! Don't understand how to prepare for cut nicknamed as "
	       << arg << endl;
	  abort();
	}
    } // loop over selection cuts

  cout << endl;

}

