#ifndef _wprime_histo_constants_h__
#define _wprime_histo_constants_h__

#include <string>

// ++++++++++++++++++++++++++++++++Useful constants
const int Num_selection_cuts = 6; // one new set of histograms after each cut
const int Num_trkAlgos = 6; // global, tracker, tpfms, cocktail, picky, tmr

// use this for histogram names
const std::string algo_desc_short[Num_trkAlgos] = {"gbl","trk","tpfms","ckt","pic","tmr"};
// use this for histogram descriptions
const std::string algo_desc_long[Num_trkAlgos] = {"global", "tracker", "TPFMS","cocktail","picky","TMR"};


//const string cuts_desc_short[Num_selection_cuts] =
//{"hlt","ptrange","1mu","iso", "jet", "qual"};
const std::string cuts_desc_short[Num_selection_cuts] = {"hlt","ptrange", "qual","1mu","iso", "jet"};
//const string cuts_desc_long[Num_selection_cuts]= {"HLT_Mu9", "Pt within
//range", "1 muon only","Isolation", "Jet Veto", "Quality"};
const std::string cuts_desc_long[Num_selection_cuts]= {"HLT_Mu9", "Pt within range", "Quality", "1 muon only","Isolation", "Jet Veto"};

// last histogram after all cuts to be used for fit
const std::string final_histo_desc = cuts_desc_short[Num_selection_cuts - 1];

#endif // #ifndef _wprime_histo_constants_h__
