#ifndef _WprimeFitter_hpp_
#define _WprimeFitter_hpp_

#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooFormulaVar.h"
#include "RooUniform.h"

#include "RooHistPdf.h"

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
#include "TH2F.h"
#include "wprimeFitter_signalDescriptions.h"


class WprimeFitter{
 public: 
  WprimeFitter(channel ch);
  ~WprimeFitter();

  void setNpseudoExperiments(unsigned N){NpseudoExp_ = N;}
  void doRunFits(bool flag){runFits_ = flag;}
  void doOneMassPointOnly(int MassPoint)
  {
    assert(MassPoint >= 0 && MassPoint < Nsignal_points);
    MassPoint_ = MassPoint;
  }

  void findOnlyMedian(bool flag){findOnlyMedian_ = flag;}
  void debugMe(bool flag){debugMe_ = flag;}

  // background option = 1 -> 1/(x+b)^c (DEFAULT)
  // background option = 2 -> 1/(x^2 + b*x + c)^d
  void setBackgroundOption(int option)
  {
    assert(option == 1 || option == 2);
    bgd_option_ = option; 
    backgroundModeled_ = false;
  }
  void skipLimitCalculation(bool flag){skipLimitCalculation_ = flag;}

  void run();

 private:
  // if true, will only calculate the 1-p_tail probability for SSM x-section,
  // but will not attempt to extract a cross-section limit
  bool skipLimitCalculation_; 
  int MassPoint_;
  vector<int> massPoints;
  bool debugMe_;
  bool findOnlyMedian_;
  bool backgroundModeled_;
  unsigned NpseudoExp_;
  bool runFits_;
  bool oneMassPointOnly_;
  float lumi_ipb_;
  int bgd_option_;
  TripleGauss * resolution[Nsignal_points];

  float eff_corr_TP_;
  
  channel channel_;
  
  TH1F * LLR[Nsignal_points]; // these correspond to ensemble of PEs generated for S+B
  TH1F * LLR_bgdOnly[Nsignal_points]; // these correspond to ensemble of PEs generated for bgd-only
  TH1F * Nsig_h[Nsignal_points];
  TH1F * Nbgd_h[Nsignal_points];
  
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
  float rXMIN_const; float rXMAX_const;
  float rXMIN[Nsignal_points]; float rXMAX[Nsignal_points];

  RooRealVar * mt;

  string file_SIG; string file_BGD; string file_data;
  TFile * fileSIG; TFile * fileBGD; TFile * fileData;

  string bgd_name; string data_name; string sig_name; string res_name;
  TH1F * bgd_hist; TH1F * data_hist; TH1F * sig_hist[Nsignal_points]; TH1F * res_hist[Nsignal_points];
  RooDataHist * mt_BGD;
  RooDataHist * mt_DATA;

  RooDataHist * mt_SIG;
  void init();
  void modelBackground();
  void modelBackgroundOption1();
  void modelBackgroundOption2();
  void modelResolutions();
  void getInputHistograms();

  void calculateExpectedZvalues(int sig_i, ofstream & tracking);
  void calculateExpectedLimits(int sig_i, ofstream & tracking);
  void calculateObservedLimit(int sig_i, ofstream & tracking);
  void calculateZvalues(int sig_i);
  void getNsig(int sig_i, float scale_factor);
  // # of bins for signal, background histograms in ROOT file
  int Nbins; 
  // # of background events in full histogram range
  float Nbgd; 
  // # of sig events in full histogram range (scaled down by scale factor)
  float Nsig; 

  RooAbsPdf * BgdPdf;
  RooRealVar * nbgd; 
  Nexp Nevt[Nsignal_points];

  RooFitResult * bgd_fr;

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
    
  // corresponds to median of LLR distribution 
  // for bgd-only ensemble of pseudo-experiments
  float Zexpect_0; 
  // correspond to +-1/2 sigma coverage points
  float Zexpect_Minus1, Zexpect_Minus2, Zexpect_Plus1, Zexpect_Plus2;
  // corresponds to expected limit and +-1/2 sigma coverage points
  float observedLimit, medianExpLimit,
    upper1sigLimit, lower1sigLimit, upper2sigLimit, lower2sigLimit;

  void initZvalues()
  {Zexpect_0 = Zexpect_Minus1 = Zexpect_Minus2 = Zexpect_Plus1 = Zexpect_Plus2
      = -1;}

  void initLimits()
  {
    observedLimit = medianExpLimit =  upper1sigLimit = lower1sigLimit = 
      upper2sigLimit = lower2sigLimit = -9999;
  }

  // try adjusting scale-factor depending on cl_test value and step_i;
  // for step_i = 0 (1, 2), scale_factor is increased by 10 (1, 0.1);
  // this method will be called till cl95 point is found
  void adjustScaleFactor(float & scale_factor, float cl_test, int & step_i);
    
  void getLLR();
  void initFit();
  void runPseudoExperiments(int sig_i, RooAbsPdf * model, 
			    RooAbsPdf & SigBgdPdf, RooRealVar &nsig,
			    bool bgdOnly);

  // if bdOnly=true, method will calculate LLR for bgd-only ensemble
  // otherwise, will calculate cl95 that corresponds to Zexpect and return
  // scale-factor (ie. x-sec(SSM)/scale-factor) for which this is achieved
  float runPseudoExperiments(int sig_i, ofstream & tracking, bool bgdOnly, float Zexpect=-1);

};

#endif // #define _WprimeFitter_hpp_
