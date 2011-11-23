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
#include "RooBgdPdf2.h"
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

  // background option = 1 -> 1/(x+b)^c (DEFAULT)
  // background option = 2 -> 1/(x^2 + b*x + c)^d
  void setBackgroundOption(int option)
  {
    assert(option == 1 || option == 2);
    bgd_option_ = option; 
    backgroundModeled_ = false;
  }

  void run();

 private:
  bool backgroundModeled_;
  unsigned NpseudoExp_;
  bool runFits_;
  float scale_factor_;
  bool oneMassPointOnly_;
  float lumi_ipb_;
  int bgd_option_;
  TripleGauss * resolution[Nsignal_points];
  
  channel channel_;
  
  TH1F * LLR[Nsignal_points];
  TH1F * Nsig_h[Nsignal_points];
  
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

  // fitting range for signal (+bgd) distribution;
  float fXMIN; float fXMAX;
  // fitting range for bgd distribution;
  float bXMIN; float bXMAX;
  // plotting range
  float pXMIN; float pXMAX;
  // full range
  float XMIN; float XMAX;
  // fitting range for resolution function
  float rXMIN; float rXMAX;

  RooRealVar * mt;

  string file_SIG; string file_BGD; string file_data;
  TFile * fileSIG; TFile * fileBGD; TFile * fileData;

  string bgd_name; string data_name; string sig_name; string res_name;
  TH1F * bgd_hist; TH1F * data_hist; TH1F * sig_hist[Nsignal_points]; TH1F * res_hist[Nsignal_points];
  RooDataHist * mt_BGD;
  RooDataHist * mt_DATA;

  void init();
  void modelBackground();
  void modelBackgroundOption1();
  void modelBackgroundOption2();
  void modelResolutions();
  void getInputHistograms();

  // # of bins for signal, background histograms in ROOT file
  int Nbins; 
  // # of background events in full histogram range
  float Nbgd; 
  // # of sig events in full histogram range (scaled down by scale factor)
  float Nsig; 
  // chi2 for each mass point, used to find observed limit
  float chi2[Nsignal_points];

  RooAbsPdf * BgdPdf;
  RooRealVar * nbgd; 
  Nexp Nevt[Nsignal_points];

  RooRealVar * bb;
  RooRealVar * cc;
  RooRealVar * dd;

  RooRealVar * mean1[Nsignal_points];
  RooRealVar * sigma1[Nsignal_points];
  RooRealVar * ff1[Nsignal_points];
  RooRealVar * mean2[Nsignal_points];
  RooRealVar * sigma2[Nsignal_points];
  RooRealVar * ff2[Nsignal_points];
  RooRealVar * mean3[Nsignal_points];
  RooRealVar * sigma3[Nsignal_points];
    
  float Zexpected;
  
  void getLLR();
  void initFit();
  void runPseudoExperiments(int sig_i, RooAbsPdf * model, 
			    RooAbsPdf & SigBgdPdf, RooRealVar &nsig);

};

#endif // #define _WprimeFitter_hpp_
