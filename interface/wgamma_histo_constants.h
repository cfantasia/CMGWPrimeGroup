#ifndef _wgamma_histo_constants_h__
#define _wgamma_histo_constants_h__

#include <string>

  // ++++++++++++++++++++++++++++++++Useful constants

 
  const int Num_photon_cuts = 5;

  const std::string photon_cuts_desc_short[Num_photon_cuts] = {"pt","isolation","hovere", "etawidth","trackveto"};
 
  // use this for histogram descriptions
  const std::string photon_cuts_desc_long[Num_photon_cuts]= {"Photon above Threshold", "Photon Isolation", "Hadronic/EM", "Eta Width", "Track Veto"};

  // +++++++++++++++++++++++++++++++tmass histogram parameters
  const unsigned nBinWTmMu = 680;
  const float minWTmMu = 0;
  const float maxWTmMu = 1700;

  // +++++++++++++++++++++++++++++++++photon pt histogram parameters
  const unsigned nBinPhoPt = 140;
  const float minPhotonPt = 0;
  const float maxPhotonPt = 800;

  // +++++++++++++++++++++++++++++++++photon eta histogram parameters
  const unsigned nBinPhoEta = 28;
  const float minPhotonEta = -2.4;
  const float maxPhotonEta = 2.4;

  // +++++++++++++++++++++++++++++++++photon phi histogram parameters
  const unsigned nBinPhoPhi = 18;
  const float minPhotonPhi = -3.6;
  const float maxPhotonPhi = 3.6;

  // +++++++++++++++++++++++++++++++++photon phi histogram parameters
  const unsigned nBinMWgamma = 680;
  const float minMWgamma = 0;
  const float maxMWgamma = 1700;

 

#endif // #ifndef _wgamma_histo_constants_h__
