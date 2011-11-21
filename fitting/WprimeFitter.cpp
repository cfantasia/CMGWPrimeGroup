#include "WprimeFitter.hpp"

using namespace RooFit;

WprimeFitter::WprimeFitter(channel ch)
{
  channel_ = ch;
  if(channel_ != wprime_MuMET && channel_ != wprime_ElMET)
    {
      cerr << " Unknown channel = " << ch << endl;
      abort();
    }

  init();

  getInputHistograms();
  
  mt = new RooRealVar("Mt", "M_{T} GeV/c^{2}", pXMIN, pXMAX);
  mt->setRange("mt_fit", fXMIN, fXMAX);
  mt->setRange("mt_plot", pXMIN, pXMAX);
  mt->setRange("mt_full", XMIN, XMAX);
  mt->setBins(10000, "fft");

  mt_BGD = new RooDataHist("mt_BGD","total BGD", *mt, Import(*bgd_hist));
  modelBackground();
  modelResolutions();
}

void WprimeFitter::init()
{
  NpseudoExp_ = 0; scale_factor_ = -1; //runFits_ = false; 
  for(unsigned i = 0; i != Nsignal; ++i)
    {
      LLR[i] = Nsig_h[i] = sig_hist[i] = res_hist[i] = 0;
      resolution[i] = 0;
      mean1[i] = mean2[i] = mean3[i] = 0;
      sigma1[i] = sigma2[i] = sigma3[i] = 0;
      ff1[i] = ff2[i] = 0;
    }

  // will need setter methods for these parameters...
  fXMIN = 220; fXMAX = 2500;
  pXMIN = fXMIN; pXMAX = 2500;
  XMIN = 0; XMAX = 2500;
  rXMIN = -2000; rXMAX = 2000;

  switch(channel_)
    {
    case wprime_MuMET:
      file_SIG = file_SIG_mu; file_BGD = file_BGD_mu; file_data = file_data_mu;
      bgd_name = bgd_name_mu; res_name = resHist_name_mu; 
      sig_name = mtHist_name_mu;
      break;
    case wprime_ElMET:
      file_SIG = file_SIG_el; file_BGD = file_BGD_el; file_data = file_data_el;
      bgd_name = bgd_name_el; res_name = resHist_name_el;
      sig_name = mtHist_name_el;
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
  
  for(unsigned sig_i = 0; sig_i != Nsignal; ++sig_i)
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
  abort();

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
  for(unsigned sig_i = 0; sig_i != Nsignal; ++sig_i)
    {
      Nevt[sig_i].Ntot =  Nevt[sig_i].Nsig = Nevt[sig_i].Nbgd = 0;
      if(LLR[sig_i]) delete LLR[sig_i];
      if(Nsig_h[sig_i]) delete Nsig_h[sig_i];
    }
  Zexpected = 99999;
}

void WprimeFitter::run()
{
  int Nmax = Nsignal;
  if(oneMassPointOnly_)Nmax = 1;

  assert(scale_factor_ > 0);

  initFit();

  for (int sig_i = 0; sig_i != Nmax; ++sig_i)
    {// loop over mass points
      
      const float mass_ = WprimeMass[sig_i];
      const float width_ = (4./3.)*(mass_/M_W)*G_W;
      // signal mass and width
      RooRealVar mass("Mass", "W' mass", mass_);//, 0, 10000);
      RooRealVar width("Width", "W' width", width_);
      
      // signal modeling: JacobianRBW
      JacobianRBWPdf sig_model("sig", "Signal", *mt, mass, width);


      Nsig = sig_hist[sig_i]->Integral(0, Nbins+1)/scale_factor_;

      RooFFTConvPdf SigPdf("SigPdf","JacobianRBW X resolution", *mt, 
			   sig_model, *(resolution[sig_i]));
      RooRealVar nsig("nsig", "# of signal events", Nsig, 0, 10000000);

      
      RooAddPdf SigBgdPdf("SigBgdPdf", "SigBgdPdf", RooArgList(SigPdf,*BgdPdf),
			  RooArgList(nsig, *nbgd));
      
      RooAbsPdf * model = 0;
      //sig_i = 0 corresponds to bgd-only ensemble
      // need a better way to make this clear-er
      if(sig_i == 0)
	model = (RooAbsPdf*) BgdPdf;
      else
	model = (RooAbsPdf*) &SigBgdPdf;

      runPseudoExperiments(sig_i, model, SigBgdPdf, nsig);

    } // loop over mass points

  getLLR();
}

void WprimeFitter::getLLR()
{
  cout << " Zexpected from average MC background only = " 
       << Zexpected << endl;

  // starting from largest mass point, as it corresponds to largest LLR values
  // this makes sure all LLR histograms will be visible

  bool showPlot = true;

  TLegend * lg = new TLegend(0.32, 0.595, 0.573, 0.896);
  lg->SetTextSize(0.04);
  lg->SetBorderSize(0);
  lg->SetFillColor(0);
  
  for(int sig_i = Nsignal-1; sig_i != -1; --sig_i)
    { // loop over mass points
	LLR[sig_i]->SetLineColor(color[sig_i]);
	float cl = -999;

	if(showPlot)
	  {
	    if(sig_i == Nsignal-1)
	      {	    
		new TCanvas();
		LLR[sig_i]->GetXaxis()->SetTitle("LLR");
		LLR[sig_i]->SetMaximum(3000);
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
	
	int bin = LLR[sig_i]->FindBin(Zexpected);
	cl = 1 - LLR[sig_i]->Integral(1, bin+1)/LLR[sig_i]->Integral();
	cout << ", Mean = " << LLR[sig_i]->GetMean();
	cout << ", 1 - P_tail(Zdata) = " << 100.0*cl << "% CL" << endl;       
    } // loop over mass points

  lg->Draw();

  cout << endl;

  for(unsigned sig_i = 0; sig_i != Nsignal; ++sig_i)
    { // loop over mass points
      cout<< " Sample # " << sig_i << ": " << desc[sig_i]; 
      cout << ", Nbgd = " << Nevt[sig_i].Nbgd << ", Nsig = " << Nevt[sig_i].Nsig
	   << ", Ntot = " << Nevt[sig_i].Ntot << endl;
      
#if 0
      string title = string("nsig") + desc[sig_i];
      new TCanvas();
      Nsig_h[sig_i]->GetXaxis()->SetTitle(title.c_str());
      Nsig_h[sig_i]->Draw();	
      
      if(sig_i == 0)
	{	    
	  new TCanvas();
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
  RooMCStudy * mcs= new RooMCStudy(*model, *mt, FitModel(SigBgdPdf), Binned(), Silence(),
		 Extended(kTRUE), FitOptions(Extended(kTRUE), 
					     PrintEvalErrors(0)));
  RooDLLSignificanceMCSModule sigModule(nsig,0);
  mcs->addModule(sigModule);
  
  //      RooRandomizeParamMCSModule randModule ;
  // randModule.sampleGaussian(nbgd,Ntot, sqrt(Ntot)) ;
  //mcs->addModule(randModule) ;  
  
  int Ntot = int(Nsig+Nbgd);
  Nevt[sig_i].Ntot = Nsig+Nbgd;
  Nevt[sig_i].Nsig = Nsig;
  Nevt[sig_i].Nbgd = Nbgd;

  // sig_i = 0 corresponds to bgd-only ensemble of pseudo-experiments
  if(sig_i == 0)
    {
      Ntot = int(Nbgd);
      Nevt[sig_i].Ntot = Nbgd;
    }

  mcs->generateAndFit(NpseudoExp_, Ntot);

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
  if(sig_i == 0)
    {
      // must replace mean by median!!!
      Zexpected = LLR[sig_i]->GetMean();
      cout << " Zexpected from average MC background only = " 
	   << Zexpected << endl;
    }
  
}

void WprimeFitter::modelBackground()
{
  RooRealVar b("b", "b", 1000, -10000, 10000);
  RooRealVar c("c", "c", 15, -100000, 100000);
  RooBgdPdf bgd_tmp("bgd_tmp", "bgd_tmp", *mt, b, c);
  
  RooPlot* xframe2 = mt->frame(Range("mt_fit"), Title("Bgd transverse mass"));

  bgd_tmp.fitTo(*mt_BGD, Range("mt_fit"), Save());
  mt_BGD->plotOn(xframe2, Name("data"));
  bgd_tmp.plotOn(xframe2, Name("model"));
  cout << " Bgd mt fit: chi2/ndof = " 
       << xframe2->chiSquare("model", "data", 2) << endl;
  cout << " Bgd mt fit: chi2/ndof = " << xframe2->chiSquare() << endl;

  bgd_tmp.paramOn(xframe2,Layout(0.55));
  new TCanvas();
  xframe2->Draw();

  RooArgList pars(* bgd_tmp.getParameters(RooArgSet(*mt) ) );

  float bb_ = ((RooRealVar*) pars.find("b"))->getVal();
  float cc_ = ((RooRealVar*) pars.find("c"))->getVal();

  bb = new RooRealVar("bb", "bb", bb_);
  cc = new RooRealVar("cc", "cc", cc_);
      
  BgdPdf = new RooBgdPdf("BgdPdf", "BgdPdf", *mt, *bb, *cc);
  Nbgd = bgd_hist->Integral(0, Nbins+1);
  nbgd = new RooRealVar("nbgd","number of background events,",Nbgd,0,100000);

}

void WprimeFitter::modelResolutions()
{
  for(unsigned i = 0; i != Nsignal; ++i){ // loop over mass points
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
  
    // RooFitResult * fr = resolution.fitTo(res_hist2, Range("resol_fit"), Save());
    res_temp.fitTo(res_hist2, Save());
    
    RooPlot * xframe = dmt.frame(Title("delta W' transverse mass"));
  
    res_hist2.plotOn(xframe, Name("data"));
    res_temp.plotOn(xframe, Name("model"));
    cout << " Sample: " << desc[i] << " Resolution fit: chi2/ndof = " 
	 << xframe->chiSquare("model", "data", 8) << endl;
    //  new TCanvas();
    // xframe->Draw();
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
