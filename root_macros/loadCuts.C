// wprime stuff
#include "UserCode/CMGWPrimeGroup/root_macros/loadCuts.h"

// returns # of (global) muons with tracker-pt above <tracker_muon_pt>
unsigned NmuAboveThresh(float tracker_muon_pt, const wprime::Event * ev)
{
  unsigned N = 0;

  int nmuons = ev->mu->GetLast() + 1;
  for(int j = 0; j != nmuons; ++j)
    { // loop over muons
      wprime::Muon * mu = (wprime::Muon *) ev->mu->At(j);
      if(mu->tracker.p.Pt() > tracker_muon_pt)
	++N;
    } // loop over muons

  return N;
}

// returns # of jets with Et above threshold
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
}

// returns # of jets with Et above threshold with angle greater than delta_phi from muon
unsigned NjetAboveThresh(float threshold, float delta_phi, 
			 const wprime::Muon * mu, const wprime::Event * ev)
{
  unsigned N = 0;
  
  int njets = ev->jet->GetLast() + 1;
  for(int j = 0; j != njets; ++j)
    { // loop over jets
      TLorentzVector * jet = (TLorentzVector *) ev->jet->At(j);
      if(jet->Et() > threshold && jet->DeltaPhi(mu->tracker.p) > delta_phi)
	++N;
    } // loop over jets

  return N;

}


