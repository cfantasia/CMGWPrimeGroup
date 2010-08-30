#ifndef _wprime_histo_constants_h__
#define _wprime_histo_constants_h__

// ++++++++++++++++++++++++++++++++Useful constants
const int Num_histo_sets = 6; // one new set of histograms after each cut
const int Num_trkAlgos = 6; // global, tracker, tpfms, cocktail, picky, tmr

// use this for histogram names
const string algo_desc_short[Num_trkAlgos] = {"gbl","trk","tpfms","ckt","pic","tmr"};
const string cuts_desc_short[Num_histo_sets] = {"hlt","ptrange","1mu","iso", "jet", "qual"};
// use this for histogram descriptions
const string algo_desc_long[Num_trkAlgos] = {"global", "tracker", "TPFMS","cocktail","picky","TMR"};
const string cuts_desc_long[Num_histo_sets]= {"HLT_Mu9", "Pt within range", 
					     "1 muon only","Isolation", 
					     "Jet Veto", "Quality"};

#endif // #ifndef _wprime_histo_constants_h__
