#include "WprimeFitter.hpp"
#include "RooDLLSignificanceMCSModule2.h"
#include <fstream>

#include "TText.h"

using namespace RooFit;

const int Nsteps=4;
const float step_size[Nsteps]={100, 10.,1., 0.1};

const float eff_corr_TP_muon = 0.974;
const float eff_corr_TP_electron = 0.960;

WprimeFitter::WprimeFitter(channel ch)
{
  channel_ = ch;
  if(channel_ != wprime_MuMET && channel_ != wprime_ElMET)
    {
      cerr << " Unknown channel = " << ch << endl;
      abort();
    }
  
  init();
  runFits_ = true; 
  getInputHistograms();
  
  mt_BGD = new RooDataHist("mt_BGD","total BGD", *mt, Import(*bgd_hist));
  mt_DATA = new RooDataHist("mt_DATA","total DATA", *mt, Import(*data_hist));
  modelResolutions();
}

void WprimeFitter::init()
{
  NpseudoExp_ = 0; 
  bgd_option_ = 1; 
  backgroundModeled_ = findOnlyMedian_ = debugMe_ = 
    skipLimitCalculation_ = false;
  MassPoint_ = -1;

  for(int i = 0; i != Nsignal_points; ++i)
    {
      LLR[i] = LLR_bgdOnly[i] = Nsig_h[i] = Nbgd_h[i] = 
	sig_hist[i] = res_hist[i] = 0;
      resolution[i] = 0;
      mean1[i] = mean2[i] = mean3[i] = 0;
      sigma1[i] = sigma2[i] = sigma3[i] = 0;
      ff1[i] = ff2[i] = 0;
    }
  
  // will need setter methods for these parameters...
  fXMIN = 380; fXMAX = 2500;

  if(channel_ == wprime_MuMET)
    {bXMIN = 380; bXMAX = 1500;} // works best with bgd-option=1, for bgd-option=2 use range (380,900)
  else if (channel_ == wprime_ElMET)
    {bXMIN = 220; bXMAX = 1500;}// works fine with bgd-option=1 

  if(fXMIN < bXMIN)fXMIN = bXMIN;

  pXMIN = min(fXMIN, bXMIN); pXMAX = fXMAX;
  rXMIN_const = -800; rXMAX_const = 800;
  if(channel_ == wprime_MuMET){
    for(int i=1; i<Nsignal_points; ++i){
      rXMIN[i] = rXMIN_const; rXMAX[i] = rXMAX_const;
    }
    rXMIN[0]=-500; rXMAX[0]=500;
  }
  else if (channel_ == wprime_ElMET){
    rXMIN[0]=rXMIN[2]=rXMIN[3]=rXMIN[4]=rXMIN[7]=rXMIN[9]=rXMIN[13]=rXMIN_const;
    rXMAX[0]=rXMAX[2]=rXMAX[3]=rXMAX[4]=rXMAX[7]=rXMAX[9]=rXMAX[13]=rXMAX_const;
    rXMIN[1]=rXMIN[8]=rXMIN[14]=rXMIN[15]=rXMIN[16]=-600;
    rXMAX[1]=rXMAX[8]=rXMAX[14]=rXMAX[15]=rXMAX[16]=600;
    rXMIN[6]=-500;
    rXMAX[6]=500;
    rXMIN[5]=rXMIN[10]=rXMIN[11]=rXMIN[12]=-500; //still needs improvement
    rXMAX[5]=rXMAX[10]=rXMAX[11]=rXMAX[12]=500;  //still needs improvement
  }

  mt = new RooRealVar("Mt", "M_{T} GeV/c^{2}", pXMIN, pXMAX);
  mt->setRange("mt_fit", fXMIN, fXMAX);
  mt->setRange("mt_bgdfit", bXMIN, bXMAX);
  mt->setRange("mt_plot", pXMIN, pXMAX);
  mt->setRange("mt_full", XMIN, XMAX);
  mt->setRange("resol_fit", rXMIN_const, rXMAX_const);
  mt->setBins(10000, "fft");
  
  switch(channel_)
    {
    case wprime_MuMET:
      file_SIG = file_SIG_mu; file_BGD = file_BGD_mu; file_data = file_data_mu;
      bgd_name = bgd_name_mu; data_name = data_name_mu; res_name = resHist_name_mu; 
      sig_name = mtHist_name_mu;
      eff_corr_TP_ = eff_corr_TP_muon;
      break;
    case wprime_ElMET:
      file_SIG = file_SIG_el; file_BGD = file_BGD_el; file_data = file_data_el;
      bgd_name = bgd_name_el; data_name = data_name_el; res_name = resHist_name_el;
      sig_name = mtHist_name_el;
      eff_corr_TP_ = eff_corr_TP_electron;
      break;
    default:
      ; // do nothing     
    }
  
}

void WprimeFitter::getInputHistograms()
{
  fileSIG   = TFile::Open(file_SIG.c_str() );
  fileBGD = TFile::Open(file_BGD.c_str() );
  fileData = TFile::Open(file_data.c_str() );
  
  if(!fileSIG || fileSIG->IsZombie()){
    cerr << " Ooops! Couldn't find file " << file_SIG << endl;
    abort();
  }
  if(!fileBGD || fileBGD->IsZombie()){
    cerr << " Ooops! Couldn't find file " << file_BGD << endl;
    abort();
  }
  if(!fileData || fileData->IsZombie()){
    cerr << " Ooops! Couldn't find file " << file_data << endl;
    abort();
  }
  
  bgd_hist = (TH1F*)fileBGD->Get(bgd_name.c_str());
  Nbins = bgd_hist->GetXaxis()->GetNbins();
  data_hist = (TH1F*)fileData->Get(data_name.c_str());
  
  for(int sig_i = 0; sig_i != Nsignal_points; ++sig_i)
    {
      string name = dirname[sig_i] + "/" + res_name;
      res_hist[sig_i] = (TH1F*) fileSIG->Get(name.c_str());
      if(!res_hist[sig_i])
	{
	  cerr << " Failed to get resolution histogram " << name << endl;
	  abort();
	}
      
      name = dirname[sig_i] + "/" + sig_name;
      sig_hist[sig_i] = (TH1F*) fileSIG->Get(name.c_str());
      if(!sig_hist[sig_i])
	{
	  cerr << " Failed to get signal histogram " << name << endl;
	  abort();
	}
      
      
      int Nbins2 = sig_hist[sig_i]->GetXaxis()->GetNbins();
      if(Nbins != Nbins2)
	{
	  cout << " Sample # " << sig_i << " (" << desc[sig_i] << ") "<< endl; 
	  cout << " Cannot handle different # of bins for sig & bgd " << endl;
	  cout << " # of signal hist bins = " << Nbins2 << endl;
	  cout << " # of background hist bins = " << Nbins << endl;
	  abort();
	}
      
    }
  
  string lumi_histo = "lumi_ipb";
  TH1F * lumi = (TH1F*) fileSIG->Get(lumi_histo.c_str());
  if(!lumi)
    {
      cerr << " Can't find histo " << lumi_histo << endl;
      abort();
    }
  lumi_ipb_ = lumi->GetBinContent(1);
  cout << " Distributions correspond to integrated luminosity = " 
       << lumi_ipb_ << " ipb " << endl;
  
}

WprimeFitter::~WprimeFitter()
{
  // if(BgdPdf) delete BgdPdf;
  
  //  delete [] resolution;
  // delete mt_BGD;
  //  delete mt;
}

void WprimeFitter::initFit()
{
  for(int sig_i = 0; sig_i != Nsignal_points; ++sig_i)
    {
      Nevt[sig_i].Ntot =  Nevt[sig_i].Nsig = Nevt[sig_i].Nbgd = 0;
      //      if(LLR[sig_i]) delete LLR[sig_i];
      // if(Nsig_h[sig_i]) delete Nsig_h[sig_i];
    }
}

// try adjusting scale-factor depending on cl_test value and step_i;
// for step_i = 0 (1, 2), scale_factor is increased by 10 (1, 0.1);
// this method will be called till cl95 point is found
void WprimeFitter::adjustScaleFactor(float & scale_factor, float cl_test, 
				     int & step_i)
{
  if(cl_test > 0.95)
    scale_factor += step_size[step_i];
  else{
    scale_factor += -1.*step_size[step_i] + step_size[step_i+1];
    step_i++;
  }
}
 

void WprimeFitter::run()
{
  if(!backgroundModeled_)
    {
      modelBackground();
      backgroundModeled_ = true;
    }
  
  massPoints.clear();
  if(MassPoint_ < 0)
    {
      for(int i = 0; i != Nsignal_points; ++i)
	massPoints.push_back(i);
    }
  else
    massPoints.push_back(MassPoint_);
    
  if(!runFits_)return;
  initFit();
  
  ofstream limits, tracking;
  string filename = "nLimit";
  string filename2 = "tracking";
  if(MassPoint_ > 0){
    filename += Form("_%.1f", WprimeMass[MassPoint_]/1000);
    filename2 += Form("_%.1f", WprimeMass[MassPoint_]/1000);
  }
  if(channel_ == wprime_MuMET){
    filename += "_MuMET.txt";
    filename2 += "_MuMET.txt";
  }
  else if(channel_ == wprime_ElMET){
    filename += "_ElMET.txt";
    filename2 += "_ElMET.txt";
  }

  limits.open(filename.c_str()); tracking.open(filename2.c_str());
  limits << "Mass/F:Lumi/F:ObsLimit/F:ExpLimit/F:ExpLimitP1/F:ExpLimitM1/F:ExpLimitP2/F:ExpLimitM2/F\n";
  tracking << "Mass\t  Zexpect\tCL\tScaleFactor\n";
  
  typedef vector<int>::const_iterator It;
  for (It sig_i = massPoints.begin(); sig_i != massPoints.end(); ++sig_i)
    {// loop over mass points

      calculateExpectedZvalues(*sig_i, tracking);
      calculateExpectedLimits(*sig_i, tracking);
      calculateObservedLimit(*sig_i, tracking);

      limits << WprimeMass[*sig_i] << '\t' << lumi_ipb_ << '\t' 
	     << observedLimit << '\t' << medianExpLimit 
	     << '\t' << upper1sigLimit << '\t' << lower1sigLimit << '\t' 
	     << upper2sigLimit << '\t' << lower2sigLimit << endl;  
      
    } // loop over mass points
  
  limits.close();
  tracking.close();

  int bin_max = bgd_hist->FindBin(pXMAX);
  cout << " # of expected bgd events: " << endl;
  cout << " Between " << pXMIN << " and " << pXMAX << " GeV = "
       << bgd_hist->Integral(bgd_hist->FindBin(pXMIN), bin_max) << endl;
  cout << " Between " << fXMIN << " and " << pXMAX << " GeV = "
       << bgd_hist->Integral(bgd_hist->FindBin(fXMIN), bin_max) << endl;
  cout << " Between " << bXMIN << " and " << pXMAX << " GeV = "
       << bgd_hist->Integral(bgd_hist->FindBin(bXMIN), bin_max) << endl;
  
  if(skipLimitCalculation_ && massPoints.size() == 1)    getLLR();
}

void WprimeFitter::getNsig(int sig_i, float scale_factor)
{
  // expected # of events given by sig_hist integral
  // first correction is from T&P (eff_corr_TP_)
  Nsig = sig_hist[sig_i]->Integral(0, Nbins+1);
  
  if(debugMe_){
    cout << "\n Running with MC sample: " << desc[sig_i] 
	 << " and scale factor = " << scale_factor << endl;
    cout << " Nbgd = " << Nbgd << endl;
    cout << " By assuming SSM cross-section, Nsig = " << Nsig << endl;
  }

  Nsig = Nsig/scale_factor;
  if(debugMe_){
    cout << " After scaling down by factor " << scale_factor 
	 << ", Nsig = " << Nsig << endl;
  }

  // correct expected # of signal events according to 
  // efficiency calculated from T&P
  Nsig = Nsig * eff_corr_TP_;
  if(debugMe_){
    cout << " After correcting for T&P-calculated efficiency, Nsig = "
	 << Nsig << endl;
  }

}

void WprimeFitter::calculateExpectedZvalues(int sig_i, ofstream & tracking)
{
  cout << "\n Calculating expected Z values for bgd-only PEs " << endl;
  const bool bgdOnly = true;
  runPseudoExperiments(sig_i, tracking, bgdOnly);
  calculateZvalues(sig_i); 
}

void WprimeFitter::calculateExpectedLimits(int sig_i, ofstream & tracking)
{
  const int NZMAX = 5;
  float Zexp[NZMAX] = {Zexpect_0, Zexpect_Minus1, Zexpect_Minus2, 
		       Zexpect_Plus1, Zexpect_Plus2};
  float xsec_limit[NZMAX] = {-1};
  int NZmax = NZMAX;
  if(findOnlyMedian_)NZmax = 1;

  const bool bgdOnly = true;
  initLimits();
  for(int z = 0; z != NZmax; ++z)
    { // loop over Z-values
      cout << "\n Calculating CL95 EXPECTED limit for Z value # " 
	   << z << " out of "
	   << NZmax << " possible values " << endl;
      float scale_factor = runPseudoExperiments(sig_i, tracking, !bgdOnly, 
						Zexp[z]);
      
      xsec_limit[z] = xsec[sig_i]/scale_factor;
    } // loop over Z-values
  
  medianExpLimit = xsec_limit[0]; 
  lower1sigLimit = xsec_limit[1]; lower2sigLimit = xsec_limit[2];
  upper1sigLimit = xsec_limit[3]; upper2sigLimit = xsec_limit[4];
}

void WprimeFitter::calculateObservedLimit(int sig_i, ofstream & tracking)
{
  // COPY + PASTE FROM runPseudoExperiments - NEED A BETTER WAY
  const float mass_ = WprimeMass[sig_i];
  const float width_ = (4./3.)*(mass_/M_W)*G_W;
  // signal mass and width
  RooRealVar mass("Mass", "W' mass", mass_);//, 0, 10000);
  RooRealVar width("Width", "W' width", width_);
  
  // signal modeling: JacobianRBW
  JacobianRBWPdf sig_model("sig", "Signal", *mt, mass, width);
  
  RooFFTConvPdf SigPdf("SigPdf","JacobianRBW X resolution", *mt, 
		       sig_model, *(resolution[sig_i]));
  SigPdf.setBufferFraction(0.6);

  // COPY + PASTE FROM runPseudoExperiments - NEED A BETTER WAY
  
  RooRealVar nsig("nsig", "# of signal events", Nsig, 0, 10000000);
  RooAddPdf SigBgdPdf("SigBgdPdf", "SigBgdPdf", RooArgList(SigPdf,*BgdPdf),
		      RooArgList(nsig, *nbgd));
  RooFitResult * rf_h1 = SigBgdPdf.fitTo(*mt_DATA, Range("mt_fit"), Save());
  RooFitResult * rf_h0 = BgdPdf->fitTo(*mt_DATA, Range("mt_fit"), Save());
  double dLL = rf_h0->minNll() - rf_h1->minNll();
  float Z_observed = dLL >= 0 ? sqrt(2*dLL) : -sqrt(-2*dLL);
  cout << " ************************************************* " << endl;
  cout << " Fit results for data distribution: " << endl;
  cout << " 2*logLike(H1) = " << 2*rf_h1->minNll();
  cout << "\n 2*logLike(H0)= " << 2*rf_h0->minNll() << endl << endl;
  cout << " 2*Delta-LL = " << dLL << endl;
  cout << " Z_observed = " << Z_observed << endl;
  cout << " ************************************************* " << endl;
  
  tracking << WprimeMass[sig_i] << '\t' << 2*rf_h0->minNll() << '\t' 
	   << 2*rf_h1->minNll() <<endl;
  


  const bool bgdOnly = true;
  float sf = runPseudoExperiments(sig_i, tracking, !bgdOnly, Z_observed);
  
  observedLimit = xsec[sig_i]/sf;
  
}

// if bgdOnly=true, method will calculate LLR for bgd-only ensemble
// otherwise, will calculate cl95 that corresponds to Zexpect and return
// scale-factor (ie. x-sec(SSM)/scale-factor) for which this is achieved
float WprimeFitter::runPseudoExperiments(int sig_i, ofstream & tracking,
					 bool bgdOnly, float Zexpect)
{      
  const float mass_ = WprimeMass[sig_i];
  const float width_ = (4./3.)*(mass_/M_W)*G_W;
  // signal mass and width
  RooRealVar mass("Mass", "W' mass", mass_);//, 0, 10000);
  RooRealVar width("Width", "W' width", width_);
  
  // signal modeling: JacobianRBW
  JacobianRBWPdf sig_model("sig", "Signal", *mt, mass, width);
  
  RooFFTConvPdf SigPdf("SigPdf","JacobianRBW X resolution", *mt, 
		       sig_model, *(resolution[sig_i]));
  SigPdf.setBufferFraction(0.6);

  //make signal model
  TH1F * hSig = (TH1F*)sig_hist[sig_i]->Clone("hSig");
  
  mt_SIG = new RooDataHist("mt_SIG","Hist Siganl", *mt, Import(*hSig));
  
  RooHistPdf Model_s("Model_sb","Signal pdf",*mt, *mt_SIG,0) ;
  RooHistPdf Model_b("Model_b","bgd only pdf",*mt, *mt_BGD,0) ;
  // for interference sample
  // intHistogram ->BG
  // virtual SIG histogram is needed, because of MCStudy()
  //    RooHistPdf Model_inf("Model_inf","interference pdf",*mt, *mt_SigInt,0) ;
  
  int step_i=0; float cl_test = 1., scale_factor=1.; 
  //int steps = 0;   
  do{
    //steps ++;
    if(bgdOnly || skipLimitCalculation_)
      scale_factor = 1;
    else
      adjustScaleFactor(scale_factor, cl_test, step_i);
	
    if(bgdOnly) Nsig=0;
    else
      getNsig(sig_i, scale_factor);
    RooRealVar nsig("nsig", "# of signal events", Nsig, 0, 10000000);

    
    //for Modeling
    RooAddPdf SigBgdhistPdf("SigBgdPdf", "SigBgdPdf", RooArgList(Model_s,Model_b),
			    RooArgList(nsig, *nbgd));
    
    //for Fitting Function
    RooAddPdf SigBgdPdf("SigBgdPdf", "SigBgdPdf", RooArgList(SigPdf,*BgdPdf),
    			RooArgList(nsig, *nbgd));

    cout << "\n Will run PE ensemble for sample " << desc[sig_i] 
	 << " bgdOnly flag = " << bgdOnly 
	 << " and scale factor = " << scale_factor << endl;
    runPseudoExperiments(sig_i, &SigBgdhistPdf, SigBgdPdf, nsig, bgdOnly);
    
    if(bgdOnly) break;

    assert(Zexpect >= 0);

    cl_test = 1 - LLR[sig_i]->Integral(1, LLR[sig_i]->FindBin(Zexpect)+1)/LLR[sig_i]->Integral();
    if(debugMe_)
      cout << "*** 1 - P_tail(" << Zexpect <<") = " << 100.0*cl_test 
	   << "% CL for scale_factor = " << scale_factor << " ***" 
	   << endl << endl;

    tracking << WprimeMass[sig_i] << '\t' << Zexpect << '\t' 
	     << cl_test << '\t' 
	     << scale_factor << endl;

    if(skipLimitCalculation_)break;

  }while(step_i != Nsteps-1 || cl_test > 0.95 ) ;
  
  return scale_factor;
}

void WprimeFitter::getLLR()
{
  TLegend * lg = new TLegend(0.32, 0.595, 0.573, 0.896);
  lg->SetTextSize(0.04);
  lg->SetBorderSize(0);
  lg->SetFillColor(0);
  
  typedef vector<int>::const_iterator It;
  for (It sig_i = massPoints.begin(); sig_i != massPoints.end(); ++sig_i)
    { // loop over mass points
      LLR[*sig_i]->SetLineColor(color[*sig_i]);
      float cl = -999;

      new TCanvas(); gPad->SetLogy();
      LLR[*sig_i]->GetXaxis()->SetTitle("LLR");
      LLR[*sig_i]->SetMaximum(3000); LLR[*sig_i]->SetMinimum(0.1); 
      LLR[*sig_i]->Draw("hist");
    
      LLR_bgdOnly[*sig_i]->SetLineColor(kBlack);
      LLR_bgdOnly[*sig_i]->Draw("hist same");
	  
      lg->AddEntry(LLR[*sig_i], desc[*sig_i].c_str());
      lg->AddEntry(LLR_bgdOnly[*sig_i], "SM");

      cout<< " Sample # " << *sig_i << ": " << desc[*sig_i]; 
      
      int bin = LLR[*sig_i]->FindBin(Zexpect_0);
      cl = 1 - LLR[*sig_i]->Integral(1, bin+1)/LLR[*sig_i]->Integral();
      cout << ", 1 - P_tail("<< Zexpect_0 << ") = " << 100.0*cl << "% CL"
	   << endl;       
    } // loop over mass points
  
  lg->Draw();
  
  cout << endl;
  
  for (It sig_i = massPoints.begin(); sig_i != massPoints.end(); ++sig_i)
    { // loop over mass points
      if(debugMe_)
	{
	  cout<< " Sample # " << *sig_i << ": " << desc[*sig_i]; 
	  cout << ", Nbgd = " << Nevt[*sig_i].Nbgd << ", Nsig = " 
	       << Nevt[*sig_i].Nsig << ", Ntot = " << Nevt[*sig_i].Ntot 
	       << endl;
	}
      if(debugMe_ && skipLimitCalculation_){

	string title = string("nsig") + desc[*sig_i];
	new TCanvas();  gPad->SetLogy();
	Nsig_h[*sig_i]->SetMaximum(10000.); Nsig_h[*sig_i]->SetMinimum(0.1); 
	Nsig_h[*sig_i]->GetXaxis()->SetTitle(title.c_str());
	Nsig_h[*sig_i]->Draw();	
	
	title = string("nbgd") + desc[*sig_i];
	new TCanvas();  gPad->SetLogy();
	Nbgd_h[*sig_i]->SetMaximum(10000.); Nbgd_h[*sig_i]->SetMinimum(0.1); 
	Nbgd_h[*sig_i]->GetXaxis()->SetTitle(title.c_str());
	Nbgd_h[*sig_i]->Draw();	
	
      }


    } // loop over mass points
  
}

void WprimeFitter::runPseudoExperiments(int sig_i, RooAbsPdf * model, 
					RooAbsPdf & SigBgdPdf, 
					RooRealVar &nsig, bool bgdOnly)
{
  Nevt[sig_i].Nsig = Nsig;
  Nevt[sig_i].Nbgd = Nbgd;
  
  float Ntot = 0;
  Ntot = Nsig+Nbgd;
  Nevt[sig_i].Ntot = Nsig+Nbgd;      

  RooMCStudy * mcs = new RooMCStudy(*model, *mt, FitModel(SigBgdPdf),
				    Silence(), Extended(kTRUE), //Binned(), 
				    FitOptions(Range("mt_fit"),Extended(kTRUE),
					       PrintEvalErrors(0)));
  
  RooDLLSignificanceMCSModule2 sigModule(nsig,0);
  mcs->addModule(sigModule);
  
  
  // correction due to lumi uncertainty (4.5%)
  float dNtot = sqrt(Ntot + 0.045 * Nsig);
  
  RooRandomizeParamMCSModule randModule ;
  randModule.sampleSumGauss(RooArgSet(nsig,*nbgd),Ntot, dNtot) ;
  mcs->addModule(randModule) ;
 
  if(debugMe_)
    mcs->generateAndFit(NpseudoExp_, 0, kTRUE);
  else
    mcs->generateAndFit(NpseudoExp_); // model from histograms (Numerical)

  if(debugMe_){
    const float mass_ = WprimeMass[sig_i];
    const float width_ = (4./3.)*(mass_/M_W)*G_W;
    RooRealVar mass("Mass", "W' mass", mass_);//, 0, 10000);
    RooRealVar width("Width", "W' width", width_);
    JacobianRBWPdf sig_model("sig", "Signal", *mt, mass, width);
    RooFFTConvPdf SigPdf("SigPdf","JacobianRBW X resolution", *mt, 
			 sig_model, *(resolution[sig_i]));
    SigPdf.setBufferFraction(0.6);

    for(int pe_num=0; pe_num<10; pe_num++){
      RooRealVar nsigH0("nsigH0", "# of signal events from H0 fit", 0);
      RooRealVar nsigH1("nsigH1", "# of signal events from H1 fit", 
			mcs->fitParams(0)->getRealValue("nsig"));
      RooRealVar nbgdH0("nbgdH0", "# of background events from H0 fit", 
			mcs->fitParams(0)->getRealValue("nbgd_H0"));
      RooRealVar nbgdH1("nbgdH1", "# of background events from H1 fit", 
			mcs->fitParams(0)->getRealValue("nbgd"));
      RooAddPdf SigBgdPdfH0("SigBgdPdfH0", "SigBgdPdfH0", RooArgList(SigPdf,*BgdPdf),
			    RooArgList(nsigH0, nbgdH0));
      RooAddPdf SigBgdPdfH1("SigBgdPdfH1", "SigBgdPdfH1", RooArgList(SigPdf,*BgdPdf),
			    RooArgList(nsigH1, nbgdH1));
      RooPlot* xframe3 = mt->frame(Range("mt_fit"), Title("Transverse mass with H0 and H1 fits for PE #0"));
      model->plotOn(xframe3, Name("model"));
      SigBgdPdfH0.plotOn(xframe3, Name("fitH0"));
      SigBgdPdfH1.plotOn(xframe3, Name("fitH1"));
      //xframe3->SetMaximum(10000); xframe3->SetMinimum(0.1);
      new TCanvas(); gPad->SetLogy();
      xframe3->Draw();
    }
  }


  if(0){
    //  if(debugMe_){
    //Find average number of entries above some threshold
    float totalEvents1 = 0;  float totalEvents2 = 0; float totalEvents3 = 0;
    for(unsigned sample_i=0; sample_i<NpseudoExp_; ++sample_i){
      TH1F * signal_hist = (TH1F*) mcs->genData(sample_i)->createHistogram("Mt");
      totalEvents1 += signal_hist->Integral(signal_hist->FindBin(fXMIN),signal_hist->GetXaxis()->GetNbins()+1);
      totalEvents2 += signal_hist->Integral(signal_hist->FindBin(bXMIN),signal_hist->GetXaxis()->GetNbins()+1);
      totalEvents3 += signal_hist->Integral(signal_hist->FindBin(220),signal_hist->GetXaxis()->GetNbins()+1);
    }
    cout << "Average number of generated events above " << fXMIN << " GeV is " << ((float)totalEvents1)/NpseudoExp_ << endl;
    cout << "Average number of generated events above " << bXMIN << " GeV is " << ((float)totalEvents2)/NpseudoExp_ << endl;
    cout << "Average number of generated events above " << 220 << " GeV is " << ((float)totalEvents3)/NpseudoExp_ << endl;
  }
  
  // Make some plots
  TH1F* z_sig = (TH1F*) mcs->fitParDataSet().createHistogram("significance_nullhypo_nsig");
  TH1F* nll0 = (TH1F*) mcs->fitParDataSet().createHistogram("nll_nullhypo_nsig");
  TH1F* nll1 = (TH1F*) mcs->fitParDataSet().createHistogram("nll_H1_nsig");
  TH1F* dll = (TH1F*) mcs->fitParDataSet().createHistogram("dll_nullhypo_nsig");
  TH1F * nsig_h = (TH1F*) mcs->fitParDataSet().createHistogram("nsig");
  TH1F * nbgd_h = (TH1F*) mcs->fitParDataSet().createHistogram("nbgd");

  if(!z_sig)
    {
      cout << " oops, z_sig = " << z_sig << " for sig_i = " << sig_i << endl;
      abort();
    }
  if(!nsig_h)
    {
      cout << " oops, nsig_h = " << nsig_h << " for sig_i = " << sig_i << endl;
      abort();
    }
  if(!nbgd_h)
    {
      cout << " oops, nbgd_h = " << nbgd_h << " for sig_i = " << sig_i << endl;
      abort();
    }
  if(!nll0)
    {
      cout << " oops, nll0 = " << nll0 << " for sig_i = " << sig_i << endl;
      abort();
    }
  if(!nll1)
    {
      cout << " oops, nll1 = " << nll1 << " for sig_i = " << sig_i << endl;
      abort();
    }
  if(!dll)
    {
      cout << " oops, dll = " << dll << " for sig_i = " << sig_i << endl;
      abort();
    }

  if(debugMe_)
    {
      cout << " Ntot = " << Ntot << endl;
      cout << " Without lumi uncertainty, dNtot = " << sqrt(Ntot) 
	   << " with lumi uncertainty, dNtot = " << dNtot << endl << endl;
      cout << " PEs created by";
      if(bgdOnly)
	cout << " bgd-oly";
      else
	cout << " signal and bgd";
      cout << " distributions " << endl;
      
      cout << " Input nsig = " << Nevt[sig_i].Nsig << ", Fit average nsig = "
	   << nsig_h->GetMean() << " +- " << nsig_h->GetRMS() << endl;
      cout << " Input nbgd = " << Nevt[sig_i].Nbgd << ", Fit average nbgd = "
	   << nbgd_h->GetMean() << " +- " << nbgd_h->GetRMS() << endl;
    } 

  if(debugMe_ && skipLimitCalculation_)
    {
      RooPlot * frame = mcs->plotPull(nsig, Bins(40), FitGauss());
      string title = "Distribution of N_{sig} pull value " + desc[sig_i];
      frame->SetTitle(title.c_str());
      new TCanvas(); frame->Draw();
      RooPlot * frame2 = mcs->plotPull(*nbgd, Bins(40), FitGauss());
      title = "Distribution of N_{bgd} pull value " + desc[sig_i];
      frame2->SetTitle(title.c_str());
      new TCanvas(); frame2->Draw();
      title = "Distribution of z_{sig} value " + desc[sig_i];
      new TCanvas(title.c_str()); z_sig->Draw();
      title = "Distribution of LL(H0) value " + desc[sig_i];
      new TCanvas(title.c_str()); nll0->Draw();
      title = "Distribution of LL(H1) value " + desc[sig_i];
      new TCanvas(title.c_str()); nll1->Draw();
      title = "Distribution of #DeltaLL value " + desc[sig_i];
      new TCanvas(title.c_str()); dll->Draw();
   }

  if(bgdOnly)
    LLR_bgdOnly[sig_i] = new TH1F(*z_sig);
  else
    {
      LLR[sig_i] = new TH1F(*z_sig);
      Nsig_h[sig_i] = new TH1F(*nsig_h);
      Nbgd_h[sig_i] = new TH1F(*nbgd_h);
    }
}

void WprimeFitter::calculateZvalues(int sig_i)
{
  initZvalues(); 
  float counted_entries=0., total_entries=LLR_bgdOnly[sig_i]->Integral();

  for(int i=1; i<=LLR_bgdOnly[sig_i]->GetNbinsX(); ++i){
    counted_entries += LLR_bgdOnly[sig_i]->GetBinContent(i);
    // cout << " i = " << i << " Z = " << LLR_bgdOnly[sig_i]->GetBinCenter(i) 
    // << " cumulative = " <<  counted_entries/total_entries << endl;
    if(Zexpect_Minus2<0 && counted_entries/total_entries>0.02275) 
      Zexpect_Minus2 = LLR_bgdOnly[sig_i]->GetBinCenter(i);
    if(Zexpect_Minus1<0 && counted_entries/total_entries>0.15865) 
      Zexpect_Minus1 = LLR_bgdOnly[sig_i]->GetBinCenter(i);
    if(Zexpect_0<0 && counted_entries/total_entries>0.5) 
      Zexpect_0 = LLR_bgdOnly[sig_i]->GetBinCenter(i);
    if(Zexpect_Plus1<0 && counted_entries/total_entries>0.84135)
      Zexpect_Plus1 = LLR_bgdOnly[sig_i]->GetBinCenter(i);
    if(Zexpect_Plus2<0 && counted_entries/total_entries>0.97725)
      Zexpect_Plus2 = LLR_bgdOnly[sig_i]->GetBinCenter(i);
  }
  
  cout << " Expected Z-values from LLR of MC background ensembles for ";
  cout << " mass point " << desc[sig_i] << endl;
  cout << " -2sigma = " << Zexpect_Minus2 << endl;
  cout << " -1sigma = " << Zexpect_Minus1 << endl;
  cout << "  Median = " << Zexpect_0 << endl;
  cout << " +1sigma = " << Zexpect_Plus1 << endl;
  cout << " +2sigma = " << Zexpect_Plus2 << endl;

}

void WprimeFitter::modelBackground()
{
  if(bgd_option_ == 1)
    modelBackgroundOption1();
  else if (bgd_option_ == 2)
    modelBackgroundOption2();
  
  //////////////////////////////////
  //Estimation scale factor for MC sideband region

  // maximum should not go above fXMIN
  float xmin = 220; float xmax = fXMIN;
  // better to use region down to bXMIN if possible (220 may be too low)
  if(bXMIN < xmax)xmin = bXMIN; 
  assert(xmin < xmax);

  double Nbgdsideband  
    = bgd_hist ->Integral( bgd_hist->FindBin(xmin),bgd_hist->FindBin(xmax) );
  double Ndatasideband 
    = data_hist->Integral( data_hist->FindBin(xmin),data_hist->FindBin(xmax) );
  
  float sf_mcdata = Nbgdsideband/Ndatasideband;
  
  cout<<" Scale factor for background from sideband region = "
      << sf_mcdata <<endl;
  
  bgd_hist->Scale(1./(sf_mcdata));
  
  float NmcBGscaled = bgd_hist->Integral(bgd_hist->FindBin(pXMIN),bgd_hist->FindBin(pXMAX));
  Nbgd = NmcBGscaled;
  
  nbgd = new RooRealVar("nbgd","number of background events,",Nbgd,0,100000);
}

void WprimeFitter::modelBackgroundOption1()
{
  RooRealVar b("b", "b", 1000, -100000, 100000);
  RooRealVar c("c", "c", 15, -1000000, 1000000);
  RooBgdPdf bgd_tmp("bgd_tmp", "bgd_tmp", *mt, b, c);
  
  RooPlot* xframe2 = mt->frame(Range("mt_bgdfit"), Title("Bgd transverse mass"));
  
  RooFitResult * rf = bgd_tmp.fitTo(*mt_BGD, Range("mt_bgdfit"), RooFit::SumW2Error (kFALSE), Save());
  mt_BGD->plotOn(xframe2, Name("data"));
  bgd_tmp.plotOn(xframe2, Name("model"));
  cout << " Bgd mt fit: chi2/ndof = " 
       << xframe2->chiSquare("model", "data", 2) << endl;
  
  bgd_tmp.paramOn(xframe2,Layout(0.55));
  xframe2->SetMaximum(10000); xframe2->SetMinimum(0.1);
  new TCanvas(); gPad->SetLogy();
  
  xframe2->Draw();
  
  RooArgList pars(* bgd_tmp.getParameters(RooArgSet(*mt) ) );
  
  float bb_ = ((RooRealVar*) pars.find("b"))->getVal();
  float cc_ = ((RooRealVar*) pars.find("c"))->getVal();
  
  bb = new RooRealVar("bb", "bb", bb_);
  cc = new RooRealVar("cc", "cc", cc_);
  
  BgdPdf = new RooBgdPdf("BgdPdf", "BgdPdf", *mt, *bb, *cc);

  float rangeMin=fXMIN, rangeMax=fXMAX;
  cout << " Actual # of entries between [" << rangeMin << ", " << rangeMax 
       << "] = " << bgd_hist->Integral(bgd_hist->FindBin(rangeMin), bgd_hist->FindBin(rangeMax)) << endl;
  cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
  
}

void WprimeFitter::modelBackgroundOption2()
{
  RooRealVar b("b", "b", -500, -100000, 100000);
  RooRealVar c("c", "c", 100000, -1000000, 1000000);
  RooRealVar d("d", "d", 3, -100000, 100000);
  RooBgdPdf2 bgd_tmp("bgd_tmp", "bgd_tmp", *mt, b, c, d);
  
  RooPlot* xframe2 = mt->frame(Range("mt_bgdfit"), Title("Bgd transverse mass"));
  
  bgd_tmp.fitTo(*mt_BGD, Range("mt_bgdfit"), RooFit::SumW2Error (kFALSE), Save());
  mt_BGD->plotOn(xframe2, Name("data"));
  bgd_tmp.plotOn(xframe2, Name("model"));
  cout << " Bgd mt fit: chi2/ndof = " 
       << xframe2->chiSquare("model", "data", 3) << endl;
  
  bgd_tmp.paramOn(xframe2,Layout(0.55));
  xframe2->SetMaximum(10000); xframe2->SetMinimum(0.1);
  new TCanvas();gPad->SetLogy();
  
    xframe2->Draw();
  
  RooArgList pars(* bgd_tmp.getParameters(RooArgSet(*mt) ) );
  
  float bb_ = ((RooRealVar*) pars.find("b"))->getVal();
  float cc_ = ((RooRealVar*) pars.find("c"))->getVal();
  float dd_ = ((RooRealVar*) pars.find("d"))->getVal();
  
  bb = new RooRealVar("bb", "bb", bb_);
  cc = new RooRealVar("cc", "cc", cc_);
  dd = new RooRealVar("dd", "dd", dd_);
  
  BgdPdf = new RooBgdPdf2("BgdPdf", "BgdPdf", *mt, *bb, *cc, *dd);
}

void WprimeFitter::modelResolutions()
{
  for(int i = 0; i != Nsignal_points; ++i){ // loop over mass points
    float xmin = res_hist[i]->GetXaxis()->GetXmin();
    float xmax = res_hist[i]->GetXaxis()->GetXmax();
    
    // Observable: transverse mass resolution
    RooRealVar dmt("dMt", "M_{T}^{RECO} - M_{T}^{MC}, GeV/c^{2}", xmin, xmax);
    string title = res_hist[i]->GetTitle();
    RooDataHist res_hist2("res_hist2",title.c_str(),dmt, Import(*res_hist[i]));
    
    RooRealVar m1("m1", "m1", 0, -10000, 10000);
    RooRealVar s1("s1", "s1", 30, -100000, 100000);
    RooRealVar m2("m2", "m2", 0, -10000, 10000);
    RooRealVar s2("s2", "s2", 100, -100000, 100000);
    RooRealVar m3("m3", "m3", 0, -10000, 10000);
    RooRealVar s3("s3", "s3", 200, -100000, 100000);
    RooRealVar f1("f1", "fraction of 1st gaussian", 0.3, 0.0, 1.0);
    RooRealVar f2("f2", "fraction of 2nd gaussian", 0.3, 0.0, 1.0);
    TripleGauss res_temp("res_temp", "triple gauss", dmt, m1, s1,
			 f1, m2, s2, f2, m3, s3);
    
    res_temp.fitTo(res_hist2, Range(rXMIN[i], rXMAX[i]), Verbose(-1), Save());
    //    res_temp.fitTo(res_hist2, Range("resol_fit"), Verbose(-1), Save()); <-- THIS DOES NOT WORK FOR SOME REASON
    //    res_temp.fitTo(res_hist2, Save());
    
    string title2 = "delta W' transverse mass (" + string(desc[i]) +") ";
    RooPlot * xframe = dmt.frame(Title(title2.c_str()));
    
    res_hist2.plotOn(xframe, Name("data"));
    res_temp.plotOn(xframe, Name("model"));
    double nchi2 = xframe->chiSquare("model", "data", 8);

    char tmp[1024];
    sprintf(tmp, " Resolution fit: chi2/ndof = %.1f", nchi2);
    string chi2_text = string(tmp);
    cout << " Sample: " << desc[i] << chi2_text << endl;

    TText* txt = new TText(-600,4000,tmp) ;
    txt->SetTextSize(0.04) ;
    txt->SetTextColor(kRed) ;
    xframe->addObject(txt) ;
    
    xframe->SetMaximum(10000);xframe->SetMinimum(0.1);
    new TCanvas();gPad->SetLogy();
    
    xframe->Draw();
    RooArgList pars(* res_temp.getParameters(RooArgSet(dmt) ) );
    
    float f_mean1 = ((RooRealVar*) pars.find("m1"))->getVal();
    float f_mean2 = ((RooRealVar*) pars.find("m2"))->getVal();
    float f_mean3 = ((RooRealVar*) pars.find("m3"))->getVal();
    float f_sigma1 = ((RooRealVar*) pars.find("s1"))->getVal();
    float f_sigma2 = ((RooRealVar*) pars.find("s2"))->getVal();
    float f_sigma3 = ((RooRealVar*) pars.find("s3"))->getVal();
    float ff1_ = ((RooRealVar*) pars.find("f1"))->getVal();
    float ff2_ = ((RooRealVar*) pars.find("f2"))->getVal();
    
    mean1[i] = new RooRealVar("mean1", "mean1", f_mean1);
    sigma1[i] = new RooRealVar("sigma1", "sigma1", f_sigma1);
    ff1[i] = new RooRealVar("ff1", "fraction of 1st gaussian",ff1_);
    mean2[i] = new RooRealVar("mean2", "mean2", f_mean2);
    sigma2[i] = new RooRealVar("sigma2", "sigma2",f_sigma2);
    ff2[i] = new RooRealVar("ff2", "fraction of 2nd gaussian", ff2_);
    mean3[i] = new RooRealVar("mean3", "mean3",f_mean3);
    sigma3[i] = new RooRealVar("sigma3", "sigma3",f_sigma3);
    
    string name = string("resolution") + Form("%.1f", WprimeMass[i]/1000);
    resolution[i] = new TripleGauss(name.c_str(), "triple gauss", *mt, 
				    *mean1[i], *sigma1[i], *ff1[i], 
				    *mean2[i], *sigma2[i], *ff2[i], *mean3[i], 
				    *sigma3[i]); 
  } // loop over mass points
}
