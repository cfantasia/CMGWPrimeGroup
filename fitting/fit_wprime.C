#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "WprimeFitter.hpp"

using namespace RooFit;

void fit_wprime()
{
  
  WprimeFitter a(wprime_MuMET);
  a.setNpseudoExperiments(1000);
  a.doOneMassPointOnly(true);

  // method WprimeFitter::run() loops over all mass points listed in 
  // wprimeFitter_signalDescription.h; user is supposed
  // to call with different values of scale-factors of interest; 
  // for scale factor <f>, the signal cross-section is the one predicted by
  // SSM, scaled down by <f>; 
  // resolution distributions/models are created once (ie. do not depend on 
  // signal cross-section); signal distributions are (re)created every time
  // method ::run() is called, presumably with a different scale factor


  // background option = 1 -> 1/(x+b)^c (DEFAULT)
  // background option = 2 -> 1/(x^2 + b*x + c)^d  
  a.setBackgroundOption(2);

  a.setScaleFactor(20);
  a.run();

}

