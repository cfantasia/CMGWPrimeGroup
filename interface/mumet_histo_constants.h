#ifndef _mumet_histo_constants_h__
#define _mumet_histo_constants_h__

#include <string>

// ++++++++++++++++++++++++++++++++Useful constants

const int Num_trkAlgos = 7; // global, tracker, tpfms, cocktail, picky, tmr, dyt
const std::string algo_desc_short[Num_trkAlgos] = {"gbl","trk","tpfms","ckt","pic","def","dyt"};
// use this for histogram descriptions
const std::string algo_desc_long[Num_trkAlgos] = {"global", "tracker", "TPFMS","cocktail","picky","default", "DYT"};

const int Num_mumet_cuts = 5; // one new set of histograms after each cut
// use this for histogram names
const std::string mumet_cuts_desc_short[Num_mumet_cuts] = {"hlt","qual","1mu","iso", "met"};

// use this for histogram descriptions
const std::string mumet_cuts_desc_long[Num_mumet_cuts]= {"Single-muon HLT", "Quality", "1 muon only","Isolation", "MET kinematic cuts"};


// +++++++++++++++++++++++++++++++muon-pt histogram parameters
const unsigned  nBinPtMu = 140;//45; // 400; // 45; // 18; 200; 380; 
const float minPtMu = 0;
const float  maxPtMu = 800; // 800; 2000;
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
const unsigned nBinTmMu = 340;
const float minTmMu = 0;
const float maxTmMu = 1700;


#endif // #ifndef _mumet_histo_constants_h__
