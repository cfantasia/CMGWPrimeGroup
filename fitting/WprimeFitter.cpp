#include "WprimeFitter.hpp"
#include <fstream>
using namespace RooFit;

const int Nsteps=3;
const float step_size[Nsteps]={10.,1., 0.1};

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
  bgd_option_ = 1; backgroundModeled_ = findOnlyMedian_ = debugMe_ = false;

  for(unsigned i = 0; i != Nsignal_points; ++i)
    {
      LLR[i] = Nsig_h[i] = sig_hist[i] = res_hist[i] = 0;
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
  XMIN = 0; XMAX = 2500;
  rXMIN = -400; rXMAX = 400;

  mt = new RooRealVar("Mt", "M_{T} GeV/c^{2}", pXMIN, pXMAX);
  mt->setRange("mt_fit", fXMIN, fXMAX);
  mt->setRange("mt_bgdfit", bXMIN, bXMAX);
  mt->setRange("mt_plot", pXMIN, pXMAX);
  mt->setRange("mt_full", XMIN, XMAX);
  mt->setRange("resol_fit", rXMIN, rXMAX);
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
  
  for(unsigned sig_i = 0; sig_i != Nsignal_points; ++sig_i)
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
  for(unsigned sig_i = 0; sig_i != Nsignal_points; ++sig_i)
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
  
  int Nmax = Nsignal_points;
  if(oneMassPointOnly_)Nmax = 1;
  
  if(!runFits_)return;

  initFit();
  
  ofstream limits, tracking;
  limits.open("nLimit.txt");
  tracking.open("tracking.txt");
  limits << "Mass/F:Lumi/F:ObsLimit/F:ExpLimit/F:ExpLimitP1/F:ExpLimitM1/F:ExpLimitP2/F:ExpLimitM2/F\n";
  tracking << "Mass\t  Zexpect\tCL\tScaleFactor\n";
  
  for (int sig_i = 0; sig_i != Nmax; ++sig_i)
    {// loop over mass points


      // First, find the Z-values that corresponds to median & 1,2 sigma points
      if(sig_i == 0)
	{
	  runPseudoExperiments(sig_i, tracking);
	  calculateZvalues(); 
	  continue;
	}


      // Then, find the expected limits
      
      const int NZMAX = 5;
      float Zexp[NZMAX] = {Zexpect_0, Zexpect_Minus1, Zexpect_Minus2, 
		       Zexpect_Plus1, Zexpect_Plus2};
      float xsec_limit[NZMAX] = {-1};
      
      int NZmax = NZMAX;
      if(findOnlyMedian_)NZmax = 1;
      
      float Z_observed = -9999, observed = -9999, median = -999, 
	upper1sig = -999, lower1sig = -999, 
	upper2sig = -9999, lower2sig = -9999;

      for(int z = 0; z != NZmax; ++z)
	{ // loop over Z-values
	  cout << "\n Calculating CL95 limit for Z value # " << z << " out of "
	       << NZmax << " possible values " << endl;
	  float scale_factor = runPseudoExperiments(sig_i, tracking, Zexp[z]);
	  
	  xsec_limit[z] = xsec[sig_i]/scale_factor;
	} // loop over Z-values

      median = xsec_limit[0]; 
      lower1sig = xsec_limit[1]; lower2sig = xsec_limit[2];
      upper1sig = xsec_limit[3]; upper2sig = xsec_limit[4];

      // Then, find the observed limits


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
      // COPY + PASTE FROM runPseudoExperiments - NEED A BETTER WAY

      Double_t nChi2H0 = 0, nChi2H1 = 0, dChi2 = 0;
      // method RooPlot::chiSquare returns chi2/(Nbins-NFITPARAM)
      // (where Nbins: # of bins that corresponds to fitting range). Then:
      // if nChi2H1 = chi2_H1/Nbins and nChi2H0 = chi2_H0/Nbins,
      // delta-chi2 := chi2_H0 - chi2_H1 = (nChi2H0-nChi2H1)*Nbins
      // where Nbins can be determined by calculating the reduced chi2 twice
      // by varying the # of dof by one.
      RooPlot* xframe_H1 = mt->frame(Range("mt_fit"), Title("Data transverse mass"));
      RooRealVar nsig("nsig", "# of signal events", Nsig, 0, 10000000);
      RooAddPdf SigBgdPdf("SigBgdPdf", "SigBgdPdf", RooArgList(SigPdf,*BgdPdf),
			  RooArgList(nsig, *nbgd));
      SigBgdPdf.fitTo(*mt_DATA, Range("mt_fit"), Save());
      mt_DATA->plotOn(xframe_H1, Name("data"));
      SigBgdPdf.plotOn(xframe_H1, Name("sigbgd_fit"));
      nChi2H1 = xframe_H1->chiSquare("sigbgd_fit", "data");
      Double_t nChi2H1_b = xframe_H1->chiSquare("sigbgd_fit", "data", 1);

      // function floor(x + 0.5) returns the closest integer to (float) x
      int Nfit_bins = int (floor(nChi2H1_b/(nChi2H1_b - nChi2H1) + 0.5));

      RooPlot* xframe_H0 = mt->frame(Range("mt_fit"), Title("Data transverse mass"));
      BgdPdf->fitTo(*mt_DATA, Range("mt_fit"), Save());
      mt_DATA->plotOn(xframe_H0, Name("data"));
      BgdPdf->plotOn(xframe_H0, Name("bgd_fit"));
      nChi2H0 = xframe_H0->chiSquare("bgd_fit", "data");
      
      cout << " Chi2H1 = " << nChi2H1*Nfit_bins 
	   << " Chi2H0 = " << nChi2H0*Nfit_bins << endl;

      dChi2 = nChi2H0*Nfit_bins - nChi2H1*Nfit_bins;
      tracking << WprimeMass[sig_i] << '\t' << nChi2H0*Nfit_bins << '\t' 
	       << nChi2H1*Nfit_bins <<endl;
      
      Z_observed = dChi2 >= 0 ? sqrt(dChi2) : 0.;

      float sf = runPseudoExperiments(sig_i, tracking, Z_observed);
      
      observed = xsec[sig_i]/sf;

      limits << WprimeMass[sig_i] << '\t' << lumi_ipb_ << '\t' 
	     << observed << '\t' << median 
	     << '\t' << upper1sig << '\t' << lower1sig << '\t' << upper2sig 
	     << '\t' << lower2sig << '\n';
      
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
  
  getLLR();
}

// if sig_i==0, method will calculate LLR for bgd-only ensemble
// otherwise, will calculate cl95 that corresponds to Zexpect and return
// scale-factor (ie. x-sec(SSM)/scale-factor) for which this is achieved
float WprimeFitter::runPseudoExperiments(int sig_i, ofstream & tracking,
					 float Zexpect)
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
  
  int step_i=0; float cl_test = 1., scale_factor=1.; 
  
  do{
    if(sig_i == 0)
      scale_factor = 1;
    else
      adjustScaleFactor(scale_factor, cl_test, step_i);
	
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
	
    RooRealVar nsig("nsig", "# of signal events", Nsig, 0, 10000000);
    
    RooAddPdf SigBgdPdf("SigBgdPdf", "SigBgdPdf", RooArgList(SigPdf,*BgdPdf),
			RooArgList(nsig, *nbgd));
    
    RooAbsPdf * model = 0;
    //sig_i = 0 corresponds to bgd-only ensemble
    // need a better way to make this clearer
    if(sig_i == 0)
      model = (RooAbsPdf*) BgdPdf;
    else
      model = (RooAbsPdf*) &SigBgdPdf;
    
    cout << "\n Will run PE ensemble for sample " << desc[sig_i] << 
      " and scale factor = " << scale_factor << endl;
    runPseudoExperiments(sig_i, model, SigBgdPdf, nsig);
    
    if(sig_i == 0)break;

    assert(Zexpect >= 0);

    cl_test = 1 - LLR[sig_i]->Integral(1, LLR[sig_i]->FindBin(Zexpect)+1)/LLR[sig_i]->Integral();
    if(debugMe_)
      cout << "*** 1 - P_tail(" << Zexpect <<") = " << 100.0*cl_test 
	   << "% CL for scale_factor = " << scale_factor << " ***" 
	   << endl << endl;

    tracking << WprimeMass[sig_i] << '\t' << Zexpect << '\t' 
	     << cl_test << '\t' 
	     << scale_factor << endl;
    
  }while(step_i != Nsteps-1 || cl_test > 0.95);
   
  return scale_factor;
}

void WprimeFitter::getLLR()
{
  // starting from largest mass point, as it corresponds to largest LLR values
  // this makes sure all LLR histograms will be visible
  
  bool showPlot = true;
  
  TLegend * lg = new TLegend(0.32, 0.595, 0.573, 0.896);
  lg->SetTextSize(0.04);
  lg->SetBorderSize(0);
  lg->SetFillColor(0);
  
  unsigned Nmax = Nsignal_points;
  if(oneMassPointOnly_)Nmax = 1;
  
  for(int sig_i = Nmax-1; sig_i != -1; --sig_i)
    { // loop over mass points
      LLR[sig_i]->SetLineColor(color[sig_i]);
      float cl = -999;
      
      if(showPlot)
	{
	  if(sig_i == Nsignal_points-1)
	    {	    
	      new TCanvas(); gPad->SetLogy();
	      LLR[sig_i]->GetXaxis()->SetTitle("LLR");
	      LLR[sig_i]->SetMaximum(3000); LLR[sig_i]->SetMinimum(0.1); 
	      LLR[sig_i]->Draw("hist");
	    }
	  else
	    {
	      LLR[sig_i]->Draw("hist same");
	    }
	  
	  lg->AddEntry(LLR[sig_i], desc[sig_i].c_str());
	}
      // show every second plot; canvas gets crowded otherwise
      showPlot = !showPlot; 
      
      cout<< " Sample # " << sig_i << ": " << desc[sig_i]; 
      
      int bin = LLR[sig_i]->FindBin(Zexpect_0);
      cl = 1 - LLR[sig_i]->Integral(1, bin+1)/LLR[sig_i]->Integral();
      cout << ", 1 - P_tail("<< Zexpect_0 << ") = " << 100.0*cl << "% CL"
	   << endl;       
    } // loop over mass points
  
  lg->Draw();
  
  cout << endl;
  
  for(unsigned sig_i = 0; sig_i != Nmax; ++sig_i)
    { // loop over mass points
      if(debugMe_)
	{
	  cout<< " Sample # " << sig_i << ": " << desc[sig_i]; 
	  cout << ", Nbgd = " << Nevt[sig_i].Nbgd << ", Nsig = " 
	       << Nevt[sig_i].Nsig << ", Ntot = " << Nevt[sig_i].Ntot << endl;
	}
#if 0
      string title = string("nsig") + desc[sig_i];
      new TCanvas();  gPad->SetLogy();
      Nsig_h[sig_i]->SetMaximum(1000.); Nsig_h[sig_i]->SetMinimum(0.1); 
      Nsig_h[sig_i]->GetXaxis()->SetTitle(title.c_str());
      Nsig_h[sig_i]->Draw();	
      
      if(sig_i == 0)
	{	    
	  new TCanvas(); gPad->SetLogy();
	  Nsig_h[0]->SetMaximum(1000.);  Nsig_h[0]->SetMinimum(0.1); 
	  Nsig_h[0]->GetXaxis()->SetTitle("nsig");
	  Nsig_h[0]->Draw();
	}
      else
	{
	  Nsig_h[sig_i]->Draw("same");
	}
#endif
    } // loop over mass points
  
}

void WprimeFitter::runPseudoExperiments(int sig_i, RooAbsPdf * model, 
					RooAbsPdf & SigBgdPdf, RooRealVar &nsig)
{
  RooMCStudy * mcs= new RooMCStudy(*model, *mt, FitModel(SigBgdPdf), 
				   Binned(), Silence(), Extended(kTRUE), 
				   FitOptions(Range("mt_fit"),Extended(kTRUE),
					      PrintEvalErrors(0)));
  RooDLLSignificanceMCSModule sigModule(nsig,0);
  mcs->addModule(sigModule);
  
  Nevt[sig_i].Nsig = Nsig;
  Nevt[sig_i].Nbgd = Nbgd;
  
  float Ntot = 0;
  if(sig_i == 0)
    {
      Ntot = Nbgd;
      Nevt[sig_i].Ntot = Nbgd;
    }
  else
    {
      Ntot = Nsig+Nbgd;
      Nevt[sig_i].Ntot = Nsig+Nbgd;      
    }

  // correction due to lumi uncertainty (4.5%)
  float dNtot = sqrt(Ntot + 0.045 * Nsig);

  if(debugMe_)
    {
      cout << " Ntot = " << Ntot << endl;
      cout << " Before lumi correction, dNtot = " << sqrt(Ntot) 
	   << " after lumi correction, dNtot = " << dNtot << endl;
    }

  // sig_i = 0 corresponds to bgd-only ensemble of pseudo-experiments
  if(sig_i == 0)
    {
      Ntot = Nbgd;
      Nevt[sig_i].Ntot = Nbgd;
      
      mcs->generateAndFit(NpseudoExp_, int(Ntot));
      
    }
  else{
    //
    RooRandomizeParamMCSModule randModule ;
    // randModule.sampleGaussian(nbgd,Ntot, sqrt(Ntot)) ;
    randModule.sampleSumGauss(RooArgSet(nsig,*nbgd),Ntot, dNtot) ;
    mcs->addModule(randModule) ;  

    if(debugMe_)
      mcs->generateAndFit(NpseudoExp_, 0, kTRUE);
    else
      mcs->generateAndFit(NpseudoExp_, 0);
    
    
    /*
      new TCanvas();
      TH1* h_sig_gen = (TH1F*) mcs->fitParDataSet().createHistogram("nsig");
      TH1* h_bgd_gen = (TH1F*) mcs->fitParDataSet().createHistogram("nbgd");
      h_sig_gen->Draw();
      new TCanvas();
      h_bgd_gen->Draw();
      new TCanvas();
      TH1* hh_nsig_nbgd  = mcs->fitParDataSet().createHistogram("hh",&nsig,YVar(*nbgd)) ;
      //TH2* hh_nsig_nbgd = (TH2*) mcs->fitParDataSet().createHistogram("nsig,nbgd") ;
      hh_nsig_nbgd->Draw();
    */

    if(debugMe_){
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

  }
  
  // Make some plots
  TH1F* dll = (TH1F*) mcs->fitParDataSet().createHistogram("dll_nullhypo_nsig");
  TH1F * nsig_h = (TH1F*) mcs->fitParDataSet().createHistogram("nsig");
  if(!dll)
    {
      cout << " oops, dll = " << dll << " for sig_i = " << sig_i << endl;
      abort();
    }
  if(!nsig_h)
    {
      cout << " oops, nsig_h = " << nsig_h << " for sig_i = " << sig_i << endl;
      abort();
    }
  
  LLR[sig_i] = new TH1F(*dll);
  Nsig_h[sig_i] = new TH1F(*nsig_h);
}

void WprimeFitter::calculateZvalues()
{
  initZvalues();
  const int sig_i = 0;
  float counted_entries=0., total_entries=LLR[sig_i]->Integral();

  for(int i=1; i<=LLR[sig_i]->GetNbinsX(); ++i){
    counted_entries += LLR[sig_i]->GetBinContent(i);

    if(Zexpect_Minus2<0 && counted_entries/total_entries>0.022) 
      Zexpect_Minus2 = LLR[sig_i]->GetBinCenter(i);
    if(Zexpect_Minus1<0 && counted_entries/total_entries>0.158) 
      Zexpect_Minus1 = LLR[sig_i]->GetBinCenter(i);
    if(Zexpect_0<0 && counted_entries/total_entries>0.5) 
      Zexpect_0 = LLR[sig_i]->GetBinCenter(i);
    if(Zexpect_Plus1<0 && counted_entries/total_entries>0.841)
      Zexpect_Plus1 = LLR[sig_i]->GetBinCenter(i);
    if(Zexpect_Plus2<0 && counted_entries/total_entries>0.977)
      Zexpect_Plus2 = LLR[sig_i]->GetBinCenter(i);
  }
  
  cout << " Expected Z-values from LLR of MC background ensembles\n";
  cout << "  Median = " << Zexpect_0 << endl;
  cout << " -1sigma = " << Zexpect_Minus1 << endl;
  cout << " +1sigma = " << Zexpect_Plus1 << endl;
  cout << " -2sigma = " << Zexpect_Minus2 << endl;
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
  
  bgd_tmp.fitTo(*mt_BGD, Range("mt_bgdfit"), RooFit::SumW2Error (kFALSE), Save());
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

  //  RooArgSet prodSet(BgdPdf);
  //RooBgdPdf unNormPdf("fitted Function", "fitted Function", prodSet);
  //TF1 * bgd_func = unNormPdf.asTF(RooArgList(x), pars);
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
  for(unsigned i = 0; i != Nsignal_points; ++i){ // loop over mass points
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
    
    res_temp.fitTo(res_hist2, Range("resol_fit"), Verbose(-1), Save());
    //    res_temp.fitTo(res_hist2, Save());
    
    RooPlot * xframe = dmt.frame(Title("delta W' transverse mass"));
    
    res_hist2.plotOn(xframe, Name("data"));
    res_temp.plotOn(xframe, Name("model"));
    cout << " Sample: " << desc[i] << " Resolution fit: chi2/ndof = " 
	 << xframe->chiSquare("model", "data", 8) << endl;
    
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
