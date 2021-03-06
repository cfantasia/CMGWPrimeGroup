#ifndef _TeVMuon_tracking_h_
#define _TeVMuon_tracking_h_

const unsigned Num_MuTeVtrkAlgos = 9; // global, tracker, tpfms, cocktail, picky, tmr, dyt, pat, standalone
const std::string algo_desc_short[Num_MuTeVtrkAlgos] = {"gbl","trk","tpfms","ckt","pic","def","dyt","pat", "std"};
// use this for histogram descriptions
const std::string algo_desc_long[Num_MuTeVtrkAlgos] = {"global", "tracker", "TPFMS","cocktail","picky","default", "DYT", "PAT", "standalone"};

const unsigned kGLOBAL = 0;
const unsigned kINNER = 1;
const unsigned kTPFMS = 2;
const unsigned kCOCKTAIL = 3;
const unsigned kPICKY = 4;
const unsigned kTEV = 5;
const unsigned kDYT = 6;
const unsigned kPAT = 7;
const unsigned kSTANDALONE = 8;
#endif // #define _TeVMuon_tracking_h_
