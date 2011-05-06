#ifndef _TeVMuon_tracking_h_
#define _TeVMuon_trackingh_

const unsigned Num_MuTeVtrkAlgos = 7; // global, tracker, tpfms, cocktail, picky, tmr, dyt
const std::string algo_desc_short[Num_MuTeVtrkAlgos] = {"gbl","trk","tpfms","ckt","pic","def","dyt"};
// use this for histogram descriptions
const std::string algo_desc_long[Num_MuTeVtrkAlgos] = {"global", "tracker", "TPFMS","cocktail","picky","default", "DYT"};

#endif // #define _TeVMuon_tracking_h_
