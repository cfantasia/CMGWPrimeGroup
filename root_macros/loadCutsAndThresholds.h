#ifndef _load_cuts_h__
#define _load_cuts_h__

#include <TROOT.h>
#include <TClonesArray.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1F.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "UserCode/CMGWPrimeGroup/interface/wprimeEvent.h"
#include "UserCode/CMGWPrimeGroup/root_macros/wprime_histo_constants.h"


// $$$$$$$$$$$$$$$$$$$$$$$ Main Study
// +++++++++++++++++++++++++++++++Analysis Thresholds:
const float EtJetCut = 100; 
const unsigned MaxNjetsAboveThresh = 0;
const float SumPtCut = 10; // Cone DeltaR =0.3; 
const float CombRelCut = 0.15; // Cone DeltaR =0.3; 
const unsigned deltaRIsoIndex = 2; //for the isolation container
const float PtTrackCut = 60 ;
const float OneMuPtTrackCut = 20; 
const float Chi2Cut = 10;
const float Delta_Phi = TMath::Pi() - 0.3;//min angle muon/jet for jet veto
const float Muon_Eta_Cut = 1.8;

#define debugme  0
#define debugmemore 0
#define dumpHighPtMuons 0

// +++++++++++++++++++++++++++++++muon-pt histogram parameters
const unsigned  nBinPtMu = 140;//45; // 400; // 45; // 18; 200; 380; 
const float minPtMu = 100; // 100;
const float  maxPtMu = 1500; // 800; 2000;
// +++++++++++++++++++++++++++++++muon-eta histogram parameters
const unsigned nBinEtaMu = 28;
const float minEtaMu = -2.4;
const float maxEtaMu = 2.4;
// +++++++++++++++++++++++++++++++muon-phi histogram parameters
const unsigned nBinPhiMu = 18;
const float minPhiMu = -3.6;
const float maxPhiMu = 3.6;
// +++++++++++++++++++++++++++++++muon-jet delphi histogram parameters
const unsigned nBinDPhiMu = 35;
const float minDPhiMu = 0;
const float maxDPhiMu = 3.5;
// +++++++++++++++++++++++++++++++muon  iso histogram parameters
const unsigned nBinIsoMu = 25;
const float minIsoMu = 0;
const float maxIsoMu = 0.5;
// +++++++++++++++++++++++++++++++tmass histogram parameters
const unsigned nBinTmMu = 195;
const float minTmMu = 50;
const float maxTmMu = 2000;
// +++++++++++++++++++++++++Declare histograms 
TH1F * hPT[Num_histo_sets][Num_trkAlgos] = {{0}};
TH1F * hETA[Num_histo_sets][Num_trkAlgos] = {{0}};
TH1F * hPHI[Num_histo_sets][Num_trkAlgos] = {{0}};
TH1F * hMJDPHI[Num_histo_sets][Num_trkAlgos] = {{0}};
TH1F * hISO[Num_histo_sets][Num_trkAlgos] = {{0}};
TH1F * hTM[Num_histo_sets][Num_trkAlgos] = {{0}};

// $$$$$$$$$$$$$$$$$$$$$ Charge Asymmetry (option = 2)
// +++++++++++++++++++++++++Declare histograms
TH1F * hPTplus[Num_trkAlgos] = {0};
TH1F * hPTminus[Num_trkAlgos] = {0};
// ++++++++++++++++++++++++++++++++Useful constants
// if 1, create charge-asymmetry histograms after all cuts
// if 2, also create histograms before last cut, and so on
const int Num_histo_sets_chargePt = 1; 






// +++++++++++++++++++++++++Declare auxiliary methods

// Calculate efficiencies
void getEff(float & eff, float & deff, float Num, float Denom);

// Get the hardest muon (based on tracker-pt) [theMu is the index in ev->mu array]
void GetTheHardestMuon(const wprime::Event * ev, int & theMu);

///returns DeltaPhi between an object and the leading jet in the event
float XJetDPhi(const TLorentzVector& lv, const wprime::Event * ev);

//transverse mass of an object with a met object
float TMass(const TLorentzVector& lv, const TVector2& themet);

// returns muon tracker isolation
float SumPtIsolation(const wprime::Muon* the_mu, unsigned detR_iso_index);

//computes muon combined relative isolation value
float CombRelIsolation(const wprime::Muon* the_mu,unsigned detR_iso_index);

// returns # of (global) muons with tracker-pt above <tracker_muon_pt>
unsigned NmuAboveThresh(float tracker_muon_pt, const wprime::Event * ev);

// returns # of jets with Et above threshold
unsigned NjetAboveThresh(float threshold, const wprime::Event * ev);

// returns # of jets with Et above 
// threshold with angle greater than
// delta_phi from muon
unsigned NjetAboveThresh(float threshold, float delta_phi, 
                         const wprime::Muon * mu, const wprime::Event * ev);


// +++++++++++++++++++++++++ Declare cut methods
// true if HLT conditions are met
bool PassedHLT(const wprime::Event* ev);

// check if muon is in pt-range for the different algorithms, fill isThere
void CheckMuonPtInRange(const wprime::Muon* mu, bool isThere[],
			float min_MuPt, float max_MuPt);

// true if only one muon with track pT > the threshold
bool OnlyOneHighTrackPtMuon(const wprime::Event* ev, float one_mu_pt_trkcut);

// true if isolation requirements satisfied for muon
bool Isolation(const wprime::Muon* the_mu, unsigned detR_iso_index, 
               float iso_cut);

// true if energetic Jet(s) found back to back with muon 
bool ExceedMaxNumJetsOpposedToMu(unsigned max_jets_aboveThresh, float et_jet_cut,
				 float delta_phi, const wprime::Muon* the_mu,
				 const wprime::Event* ev);

// check if muon satisfies quality requirements for all tracking algorithms, 
// fill goodQual
void CheckQuality(const wprime::Muon* mu, bool goodQual[],
		  float pttrk_cut, float chi2_cut, float muon_etacut);


#endif // #define _load_cuts_h__
