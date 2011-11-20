#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "WprimeFitter.hpp"

using namespace RooFit;

void fit_wprime()
{
  
  WprimeFitter a(wprime_ElMET);
  a.setNpseudoExperiments(1000);
  a.doOneMassPointOnly(false);

  // method WprimeFitter::run() loops over all mass points listed in 
  // wprimeFitter_signalDescription.h; user is supposed
  // to call with different values of scale-factors of interest; 
  // for scale factor <f>, the signal cross-section is the one predicted by
  // SSM, scaled down by <f>; 
  // resolution distributions/models are created once (ie. do not depend on 
  // signal cross-section); signal distributions are (re)created every time
  // method ::run() is called, presumably with a different scale factor

  a.setScaleFactor(10);
  a.run();

}

