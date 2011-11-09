#ifndef _elmet_histo_constants_h__
#define _elmet_histo_constants_h__

#include <string>

// ++++++++++++++++++++++++++++++++Useful constants

const int Num_elmet_cuts = 5; // one new set of histograms after each cut
// use this for histogram names
const std::string elmet_cuts_desc_short[Num_elmet_cuts] = {"hlt","qual","1el","iso", "met"};

// use this for histogram descriptions
const std::string elmet_cuts_desc_long[Num_elmet_cuts]= {"Single-electron HLT", "Quality", "1 electron only","Isolation", "MET kinematic cuts"};


// +++++++++++++++++++++++++++++++electron-pt histogram parameters
const unsigned  nBinEtEle = 750;//45; // 400; // 45; // 18; 200; 380; 
const float minEtEle = 0;
const float  maxEtEle = 1500; // 800; 2000;
// +++++++++++++++++++++++++++++++electron-eta histogram parameters
const unsigned nBinEtaEle = 28;
const float minEtaEle = -2.5;
const float maxEtaEle = 2.5;
// +++++++++++++++++++++++++++++++electron-phi histogram parameters
const unsigned nBinPhiEle = 18;
const float minPhiEle = -3.6;
const float maxPhiEle = 3.6;
// +++++++++++++++++++++++++++++++electron-jet delphi histogram parameters
const unsigned nBinDPhiEle = 35;
const float minDPhiEle = 0;
const float maxDPhiEle = 3.5;
// +++++++++++++++++++++++++++++++electron  iso histogram parameters
const unsigned nBinIsoEle = 25;
const float minIsoEle = 0;
const float maxIsoEle = 0.5;
// +++++++++++++++++++++++++++++++tmass histogram parameters
const unsigned nBinTmEle = 1250;
const float minTmEle = 0;
const float maxTmEle = 2500;


#endif // #ifndef _elmet_histo_constants_h__
