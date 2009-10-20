#ifndef _common_util_h_
#define _common_util_h_

#include "UserCode/CMGWPrimeGroup/interface/wprimeEvent.h"

#include <string>

// cone size and SumPt for isolation
/*
const unsigned nBinCone = wprime::N_CONESIZE;
const float minCone = wprime::MIN_CONE;
const float maxCone = wprime::MAX_CONE;
*/

// all cuts are in in GeV
extern const float EtJetCut;
extern const float SumPtCut;
extern const float PtTrackCut;
extern const float OneMuPtTrackCut;

extern const float Chi2Cut;

// minimum angle between muon and jet (for jet-activity veto)
extern const float Delta_Phi;

// eta-cut for muon track
extern const float Muon_Eta_Cut;

static const int Num_trkAlgos = 3; // global, tracker, tev_1st

static const string algo_desc[Num_trkAlgos] = 
  {" global", " tracker", " TeV-1st"};

// muon-pt histogram parameters
extern const unsigned  nBinPtMu;
extern const float minPtMu;
extern const float  maxPtMu;

void getEff(float & eff, float & deff, float Num, float Denom);

#endif
