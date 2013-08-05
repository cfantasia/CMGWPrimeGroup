//Author: Cory Fantasia 2013
//Purpose: Plot limit vs number of toys for test purposes
//Usage: root -b -l -q plotLimitVsNToys.C+

#include "../root_macros/common.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include <fstream>

double getMedian(TTree* tree, const string & cuts, double & lo95,double & lo68,double & hi68,double & hi95, int ntoys);

void
plotLimitVsNToys(){
  TGraph* gxsec = new TGraph("../Limits/xSec_WZ.dat", "%lg %lg");

  TGraphAsymmErrors* gNToys = new TGraphAsymmErrors();
  string inputFile = Form("higgsCombine-TestNToys.mH1900.root");
  
  TFile *fIn = TFile::Open(inputFile.c_str(), "read"); if(!fIn) return;
  TTree* tree = (TTree*) fIn->Get("limit"); if(!tree) return;
  
  for(int ntoys=1, ipoint=0; ntoys<=10000; ntoys*=2, ipoint++){

    double medianLimit,lo95,lo68,hi68,hi95;
    medianLimit = lo95 = lo68 = hi68 = hi95 = 0;
    
    medianLimit = getMedian(tree, "iToy>0", lo95,lo68,hi68,hi95, ntoys);

    printf("point %i ntoys %i med %f high %f low %f\n", ipoint, ntoys, medianLimit, hi68,lo68); 

    medianLimit *= gxsec->Eval(1900);
    lo95 *= gxsec->Eval(1900);
    lo68 *= gxsec->Eval(1900);
    hi68 *= gxsec->Eval(1900);
    hi95 *= gxsec->Eval(1900);

    gNToys->SetPoint      (ipoint, ntoys, medianLimit);
    gNToys->SetPointEYhigh(ipoint, hi95);
    gNToys->SetPointEYlow (ipoint, lo95);

  }
  gNToys->Draw("alp");
}

double
getMedian(TTree* tree, const string & cuts, double & lo95,double & lo68,double & hi68,double & hi95, int ntoys){
  vector<float> limitHistory;
  
  int nLimits = tree->Draw("limit", cuts.c_str());
  if (nLimits <= 0){
    cout<<" Failed to find any lines in tree\n";
    return 0.;
  }
  nLimits = min(ntoys, nLimits);
  cout<<"Using this many toys "<<nLimits<<endl;
  Double_t*  limit = tree->GetVal(0);
  for(int ievt=0; ievt<nLimits; ++ievt){
    limitHistory.push_back(limit[ievt]);
  }
  sort(limitHistory.begin(), limitHistory.end());
  
  double medianLimit = (nLimits % 2 == 0 ? 0.5*(limitHistory[nLimits/2-1]+limitHistory[nLimits/2]) : limitHistory[nLimits/2]);
  hi95 = limitHistory[min<int>(nLimits-1,  ceil(0.975 * nLimits))];
  hi68 = limitHistory[min<int>(nLimits-1,  ceil(0.84  * nLimits))];
  lo68 = limitHistory[min<int>(nLimits-1, floor(0.16  * nLimits))];
  lo95 = limitHistory[min<int>(nLimits-1, floor(0.025 * nLimits))];
    //cout << "   68% expected band : " << lo68 << " < r < " << hi68 << endl;
    //cout << "   95% expected band : " << lo95 << " < r < " << hi95 << endl;

  cout<<"limit is "<<medianLimit<<endl;

  return medianLimit;
}
