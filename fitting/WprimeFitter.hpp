#ifndef _WprimeFitter_hpp_
#define _WprimeFitter_hpp_


#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "RooDataHist.h"
#include "RooClassFactory.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
#include "TAxis.h"
#include "TF1.h"
#include "TLegend.h"
#include "RooPlot.h"
#include "TripleGauss.h"
#include "RooProduct.h"
#include "JacobianRBWPdf.h"
#include "RooBgdPdf.h"
#include "RooFFTConvPdf.h"
#include "RooDLLSignificanceMCSModule.h"
#include "RooRandomizeParamMCSModule.h"
#include "RooMCStudy.h"

#include "TripleGauss.h"
#include "RooBgdPdf.h"
#include "JacobianRBWPdf.h"

#include "TFile.h"
#include "TH1F.h"

#include "wprimeFitter_signalDescriptions.h"


class WprimeFitter{
 public: 
  WprimeFitter(channel ch);
  ~WprimeFitter();

  void setNpseudoExperiments(unsigned N){NpseudoExp_ = N;}
  //  void doRunFits(bool flag){runFits_ = flag;}
  void setScaleFactor(float factor){scale_factor_ = factor;}
  void doOneMassPointOnly(bool flag){oneMassPointOnly_ = flag;}

  void run();

 private:
  unsigned NpseudoExp_;
  bool runFits_;
  float scale_factor_;
  bool oneMassPointOnly_;
  float lumi_ipb_;

  TripleGauss * resolution[Nsignal];
  
  channel channel_;
  
  TH1F * LLR[Nsignal];
  TH1F * Nsig_h[Nsignal];
  
  void getBins(TH1F * & h, float xmin, float xmax, 
	       int & bin_first, int  & bin_last)
  {
    assert(h != 0);
    bin_first = h->FindBin(xmin);
    bin_last = h->FindBin(xmax);
  }

  // put_here contains 2 parameters for the parameterization of RooBgdPdf;
  // need to be generalized to deal with different background models
  void getBackground(RooRealVar & mt, RooDataHist & mt_BGD, float * put_here);

  // fitting range for signal distribution;
  float fXMIN; float fXMAX;
  // plotting range
  float pXMIN; float pXMAX;
  // full range
  float XMIN; float XMAX;
  // fitting range for resolution function
  float rXMIN; float rXMAX;

  RooRealVar * mt;

  string file_SIG; string file_BGD; string file_data;
  TFile * fileSIG; TFile * fileBGD; TFile * fileData;

  string bgd_name; string sig_name; string res_name;
  TH1F * bgd_hist; TH1F * sig_hist[Nsignal]; TH1F * res_hist[Nsignal];
  RooDataHist * mt_BGD;

  void init();
  void modelBackground();
  void modelResolutions();
  void getInputHistograms();

  // # of bins for signal, background histograms in ROOT file
  int Nbins; 
  // # of background events in full histogram range
  float Nbgd; 
  // # of sig events in full histogram range (scaled down by scale factor)
  float Nsig; 

  RooBgdPdf * BgdPdf;
  RooRealVar * nbgd; 
  Nexp Nevt[Nsignal];

  RooRealVar * bb;
  RooRealVar * cc;

  RooRealVar * mean1[Nsignal];
  RooRealVar * sigma1[Nsignal];
  RooRealVar * ff1[Nsignal];
  RooRealVar * mean2[Nsignal];
  RooRealVar * sigma2[Nsignal];
  RooRealVar * ff2[Nsignal];
  RooRealVar * mean3[Nsignal];
  RooRealVar * sigma3[Nsignal];
    
  float Zexpected;
  
  void getLLR();
  void initFit();
  void runPseudoExperiments(int sig_i, RooAbsPdf * model, 
			    RooAbsPdf & SigBgdPdf, RooRealVar &nsig);

};

#endif // #define _WprimeFitter_hpp_
