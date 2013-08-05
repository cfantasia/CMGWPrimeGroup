//Author: Cory Fantasia 2012
//Purpose: Plot resolution of fixed pt
//Usage: root -b -l -q plotRecoMass.C

#include "common.h"
//#include "UserCode/CMGWPrimeGroup/root_macros/common.h"

void
plotRecoMass(){

  TFile *f = TFile::Open("../../../WprimeWZ.root", "read"); assert(f);

  for(int mass=200; mass<=2000; mass+=100){
    string sample = Form("WprimeToWZTo3LNu_M-%i", mass);
    
    TH1F hRes("hRes", "Gen Res", 4000, mass-2000., mass+2000.);
    //hRes.StatOverflows(kTRUE);

    get_sum_of_hists(f, sample,
                     "tEvts_MET", "WZMass",
                     "", hRes);

    int bin1 = hRes.FindFirstBinAbove(hRes.GetMaximum()/2);
    int bin2 = hRes.FindLastBinAbove(hRes.GetMaximum()/2) + 1;
    double fwhm = hRes.GetBinCenter(bin2) - hRes.GetBinCenter(bin1);
      
    printf("Mass %i Mean %.2f Res %.2f FWHM %.2f\n", mass, hRes.GetMean(), hRes.GetRMS(), fwhm);
  }
}
