#ifndef _fit_wprime_h_
#define _fit_wprime_h_

#include <TH1F.h>

// fitting range is defined by fXMIN, fXMAX;
// background is not described well by simple functions 
// (exponential, landau, RBW) for pT < 150; limit fit to region above fXMIN
// fXMAX should be something reasonable (ie. up to when we run out of statistics)
const float fXMIN = 120;
const float fXMAX = 700; // 850; // 1200;//850;
// in GeV
const Double_t width_W = 2.196;
const Double_t mass_W = 80.398;

// set bin-size, should be equal to histogram bin-size we fit
void setBinSize(float bin_size);
// set resolution function
void setResolution(TH1F * g);

// set landau (true) or RBW (false) background
void setLandauBgd(bool flag);

// this is the "auxiliary", non-normalized relativistic Breight Wigner function; 
Double_t RBW_aux(Double_t * x, Double_t * par);

// if dN/dcos(theta) = (1 + costheta^2) or (1 +- costheta)^2
// use this to normalize part-(a) of function mySig below
Double_t mySig_aux(Double_t * x, Double_t * par);

// relativistic Breit Wigner (normalized) function;
// to be used for description of background (between fXMIN, fXMAX)
Double_t myRBW(Double_t * x, Double_t * par);

// convolution of (a) dN/dcos(theta) = (1 + costheta^2) [or (1 +- costheta)^2]
// and (b) relativistic Breit Wigner
Double_t mySig(Double_t * x, Double_t * par);

// normalized Landau function
Double_t myLandau(Double_t * x, Double_t * par);

Double_t myExp(Double_t * x, Double_t * par);

// convolution of (a) dN/dcos(theta) = (1 + costheta^2) [or (1 +- costheta)^2]
// and (b) MC-true mass distribution of wprime
// Double_t mySig2(Double_t * x, Double_t * par)

Double_t my_gauss(Double_t * x, Double_t * par);

/* smear function "func2smear" (global pointer); 
   need to set pointer before calling smear_func 
   (should be done by spefic smearing function, e.g. smeared_exp) */
Double_t smear_func(Double_t * x, Double_t * par);

Double_t smeared_sig(Double_t * x, Double_t * par);

Double_t myBgd(Double_t * x, Double_t * par);

Double_t mySigBgd(Double_t * x, Double_t * par);

#endif //_fit_wprime_h_
