#include "UserCode/CMGWPrimeGroup/root_macros/util.h"

// all cuts are in in GeV
const float EtJetCut = 100; 
const float SumPtCut = 10; // 3; 
const float PtTrackCut = 60 ;
const float OneMuPtTrackCut = 20; 

const float Chi2Cut = 10;

// minimum angle between muon and jet (for jet-activity veto)
const float Delta_Phi = TMath::Pi() - 0.3;

// eta-cut for muon track
const float Muon_Eta_Cut = 1.8;


// muon-pt histogram parameters
const unsigned  nBinPtMu = 45; // 400; // 45; // 18; 200; 380; 
const float minPtMu = 100; // 100;
const float  maxPtMu = 1500; // 800; 2000;


void getEff(float & eff, float & deff, float Num, float Denom)
{
  eff = Num/Denom;
  deff = TMath::Sqrt(eff * (1-eff)/Denom);
}
