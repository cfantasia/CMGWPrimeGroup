#ifndef _wprime_histo_constants_h__
#define _wprime_histo_constants_h__

#include <string>

// ++++++++++++++++++++++++++++++++Useful constants
const int Num_selection_cuts = 8; // one new set of histograms after each cut
const int Num_trkAlgos = 7; // global, tracker, tpfms, cocktail, picky, tmr, dyt
// making distributions for two flavors: muon-pt and Mt
const int Num_flavors = 2;
const int PT_INDEX = 0;
const int MT_INDEX = 1;

const string FLAVOR_NAME[Num_flavors] = {" Muon-pt analysis", " Mt analysis"};

// lower, upper limits for algorithms to consider when processing data
// default: processing only cocktail tracking (3)
const int MuAlgo_MIN = 3;
const int MuAlgo_MAX = 3;

// lower, upper limits for cuts that specify Pt/Mt thresholds
const int THRESH_MIN = 5;
const int THRESH_MAX = 7;

// use this for histogram names
const std::string algo_desc_short[Num_trkAlgos] = {"gbl","trk","tpfms","ckt","pic","tmr","dyt"};
// use this for histogram descriptions
const std::string algo_desc_long[Num_trkAlgos] = {"global", "tracker", "TPFMS","cocktail","picky","TMR", "DYT"};

// use this for histogram names
const std::string cuts_desc_short[Num_selection_cuts] = {"hlt","qual","1mu","iso", "jetmet", "thr1", "thr2", "thr3"};
// use this for histogram descriptions
const std::string cuts_desc_long[Num_selection_cuts]= {"Single-muon HLT", "Quality", "1 muon only","Isolation", "Jet Veto OR kinematic cuts", "Threshold 1", "Threshold 2", "Threshold 3"}; 

// last histogram after all cuts to be used for fit
const std::string final_histo_desc = cuts_desc_short[Num_selection_cuts - 1];

#endif // #ifndef _wprime_histo_constants_h__
