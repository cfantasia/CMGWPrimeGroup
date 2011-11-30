#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "WprimeFitter.hpp"

using namespace RooFit;

void fit_wprime()
{
  
  WprimeFitter a(wprime_MuMET);
  a.setNpseudoExperiments(1000);
  a.doOneMassPointOnly(false);
  a.doRunFits(true);
  a.findOnlyMedian(false);
  
  /* method WprimeFitter::run() loops over all mass points listed in 
     wprimeFitter_signalDescription.h; for each mass point it 
     varies the nominal/SSM W' cross-section by scanning over the parameter 
     space for the scale factor f; for each mass-point and value f the code 
     determines the LLR(H0, H1)  distribution and the 1-p_tail integral 
     for Z_expected = median of bgd-only LLR distribution. 
     When the 95 % CL = 1-p_tail point has been found, 
     this is the upper limit for the given mass point; 
     
     Some code details: 
     o resolution distributions/models are created once (ie. do not depend on 
     signal cross-section); 
     o signal distributions are (re)created every time
     method ::runPseudoExperiments() is called, presumably with a 
     different scale factor
  */

  // background option = 1 -> 1/(x+b)^c (DEFAULT)
  // background option = 2 -> 1/(x^2 + b*x + c)^d  
  a.setBackgroundOption(1);

  a.run();

}

