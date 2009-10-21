#ifndef _results_h_
#define _results_h_

#include <TObject.h>

class Results : public TObject {
 public:
  /// chi^2
   float chi2;
   /// experiment number
   Int_t exp_no;
   /// Ndof
   Int_t Ndof;
   /// probability
   float Prob;
   /// fit status (TMinuit::GetStatus())
   Int_t status;
   /// Bgd parameters
   float p0;
   float dp0;
   float p1;
   float dp1;
   float p2;
   float dp2;
   /// integral of bgd function in fit range		
     /// float bgd_tot; 
   /// Sig parameters
   float mass;
   float dmass;
   float width;
   float dwidth;
   float Nsig;
   float dNsig;
   float f; // fudge factor: increases predicted width by f
   float df;

     Results();
     ~Results();
   
   ClassDef(Results, 1)
};

#endif // #define _results_h_

