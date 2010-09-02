#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TMinuit.h>
#include <TLatex.h>
#include <TVirtualFitter.h>

#include <string>
#include <iostream>
#include <sstream>

#include "fit_wprime.h"
#include "common_fit.h"
#include "Results.h"
#include "../root_macros/wprime_histo_constants.h"

using std::string; using std::cout; using std::endl;

// these are the "running" histograms that change on event-by-event basis for fit;
// to be drawn by original (reference) histograms
TH1F * wprime = 0;
TH1F * bgd = 0;
TH1F * tot = 0;
TF1 * ftot = 0;

float Ntot = 0;

Int_t Nsig_input = 0;

static int npoints = -999;
static int status = -9999;

// # of function parameters to be employed for fit
unsigned NPARAM_FIT = 0;

void setNPARAM_FIT(unsigned num){NPARAM_FIT = num;}

// get fit results, store in Results structure
void getResults(Results * result, TF1 * f, int exp_no)
{
  result->exp_no = exp_no;
  result->chi2 = f->GetChisquare();
  result->Ndof = f->GetNDF();
  result->Prob = f->GetProb();
  result->status = status;
  result->p0 = f->GetParameter(0);  result->dp0 = f->GetParError(0);
  result->p1 = f->GetParameter(1);  result->dp1 = f->GetParError(1);
  result->p2 = f->GetParameter(2);  result->dp2 = f->GetParError(2);
  result->Nsig_input = Nsig_input;
  result->Nsig = f->GetParameter(3); result->dNsig = f->GetParError(3); 
  result->mass = f->GetParameter(4); result->dmass = f->GetParError(4); 
  //  result->width = f->GetParameter(5); result->dwidth = f->GetParError(5); 
  result->f = f->GetParameter(5); result->df = f->GetParError(5); 

  cout << " =================== Results: ==========================" << endl;
  cout << " Nsig = " << result->Nsig << " +- " << result->dNsig;
  cout << " (Input: " << result->Nsig_input << ") " << endl;
  cout << " Mass = " << result->mass << " +- " << result->dmass << endl;
  cout << " f = " << result->f << " +- " << result->df << endl;
  cout << " Chi2/Ndof = " << result->chi2 << "/" << result->Ndof 
       << ", Prob = " << result->Prob <<", Status = " << result->status 
       << endl;
}

void makeHistogram(TH1F * & makeThis, const char * hname, const char * htitle, 
		   TH1F * & likeThis)
{
  makeThis = new TH1F(hname, htitle,
		       likeThis->GetNbinsX(),
		       likeThis->GetXaxis()->GetXmin(), 
		       likeThis->GetXaxis()->GetXmax());
  makeThis->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  // do not attach to any directory, so we can delete ourselves at end of loop
  makeThis->AddDirectory(0);
}

void populateHistogram(TH1F * makeThis, unsigned N, TH1F * likeThis, 
		       bool store_Nsig_input)
{
  if (N==0)return;
  // add poisson flucturation in the # of events for i-th pseudo-experiment
  TF1 f("mypoiss", "TMath::Poisson(x, [0])", 0, 5.*N);
  f.SetParameter(0, N);
  unsigned N2 = unsigned(f.GetRandom());
  if(store_Nsig_input)
    Nsig_input = N2;
  for(unsigned j = 0; j != N2; ++j)
    makeThis->Fill(likeThis->GetRandom());
}

// custom chi2
// could be either Almeida-Barbi-do Vale or Baker-Cousis implementation
void myFCN(Int_t& /*nPar*/, Double_t * /*grad*/, Double_t &fval, Double_t *par, 
	   Int_t /*iflag */  ) 
{
  Double_t chi2 = 0;
  npoints = 0;
  TAxis *xaxis = tot->GetXaxis();
  int bin_first = xaxis->FindBin(fXMIN); int bin_last = xaxis->FindBin(fXMAX);
  for(int bin_no = bin_first; bin_no <= bin_last; ++bin_no)
    {
      Double_t x = xaxis->GetBinCenter(bin_no);
      Double_t data = tot->GetBinContent(bin_no);
      Double_t theory = ftot->EvalPar(&x, par);
      if(theory < 0)
	theory = -theory;
      
      if(ABV_chi2)
	// Almeida, Barbi, do Vale, NIM A449 (2000), 383-395
	chi2 += 2*(theory - data) + (2*data + 1) * 
	  TMath::Log( (2*data + 1)/(2*theory + 1));      
      else
	{
	  // traditional log-likelihood function (aka Baker-Cousins 1984)
	  if (theory < 1.E-9)theory = 1.E-9;
	  if (data < 1.E-9) data = 1.E-9;
	  chi2 += 2*(theory - data + data*TMath::Log(data/theory));
	}

      ++npoints;
    }
  fval = chi2;
}

void fitIt(TH1F * data, TF1 * theory)
{
  tot = data; ftot = theory;
  if(custom_chi2)
    {
      TVirtualFitter * fitter = TVirtualFitter::Fitter(0,NPARAM_FIT);
      for (unsigned i = 0; i != NPARAM_FIT; ++i) 
	fitter->SetParameter(i, ftot->GetParName(i), 
			     ftot->GetParameter(i), 20, 0,0);

      if(fixNtot)
	{
	  fitter->SetParameter(0, ftot->GetParName(0), 
			       ftot->GetParameter(0), 0, Ntot, Ntot);
	}

      //      fitter->SetParameter(1, ftot->GetParName(1), ftot->GetParameter(1),
      //			   0, fXMIN, fXMIN);

      if(NPARAM_FIT > 3)
	{
	  if(fixFudge)
	    {
	      fitter->SetParameter(5, ftot->GetParName(5), 
				   ftot->GetParameter(5),0, Fudge, Fudge);
	    }
	  if(use_wprime_mass_limits)
	    {
	      fitter->SetParameter(4, ftot->GetParName(4), 
				   ftot->GetParameter(4), 20, 
				   0, upper_wprime_mass_limit);
	    }
	}
      
      fitter->SetFCN(myFCN);
      //      TVirtualFitter::Fitter(tot)->SetFCN(myFCN);

      double arglist[100];
      arglist[0] = 0;
      // set print level
      fitter->ExecuteCommand("SET PRINT",arglist,1);

      arglist[0] = 1;
      fitter->ExecuteCommand("SET ERR", arglist, 1);

      if(use_wprime_mass_limits)
	{
	  // limits on the wprime mass
	  arglist[0] = 5;
	  arglist[1] = 0;
	  arglist[2] = upper_wprime_mass_limit;
	  fitter->ExecuteCommand("SET LIM", arglist, 3);
	}

      if(fixNtot)
	{
	  // fix total # of events
	  arglist[0] = 1;
	  arglist[1] = Ntot;
	  arglist[2] = Ntot;
	  fitter->ExecuteCommand("SET LIM", arglist, 3);
	}

      if(fixFudge && NPARAM_FIT > 3)
	{
	  // fix "fudge factor" for W' width
	  arglist[0] = 6;
	  arglist[1] = Fudge;
	  arglist[2] = Fudge;
	  fitter->ExecuteCommand("SET LIM", arglist, 3);
	}

      // minimize
      arglist[0] = 5000; // number of function calls
      arglist[1] = 0.01; // tolerance
      fitter->ExecuteCommand("MIGRAD",arglist,2);

      //get result
      double minParams[NPARAM_FIT];
      double parErrors[NPARAM_FIT];
      for (unsigned i = 0; i != NPARAM_FIT; ++i) 
	{  
	  minParams[i] = fitter->GetParameter(i);
	  parErrors[i] = fitter->GetParError(i);
	}
      double chi2, edm, errdef; 
      int nvpar, nparx;
      status = fitter->GetStats(chi2,edm,errdef,nvpar,nparx);
      ftot->SetParameters(minParams);
      ftot->SetParErrors(parErrors);
      ftot->SetChisquare(chi2);
      ftot->SetNumberFitPoints(npoints);
      int ndf = npoints-nvpar;
      ftot->SetNDF(ndf);
      tot->GetListOfFunctions()->Add(ftot);      
    }
  else
    { // use MINUIT built-in log-likelihood
      tot->Fit(ftot, "LRV");
      double chi2, edm, errdef; 
      int nvpar, nparx;
      cout << " status = " << TVirtualFitter::Fitter(tot)->GetStats(chi2,edm,errdef,nvpar,nparx) << endl;
      cout << " Chi2/Ndof = " << chi2 << "/" <<npoints-nvpar << endl;
    }

}

void doPlots(TH1F * data, TF1 * theory, bool bgdOnlyFit)
{

  TF1 fbgd("mybgd", myBgd, fXMIN, fXMAX, 3);
  fbgd.SetParameters(theory->GetParameter(0) -theory->GetParameter(3),
		     theory->GetParameter(1), theory->GetParameter(2));
  TF1 fsig("mysig", smeared_sig, fXMIN, fXMAX, 3);
  fsig.SetParameters(theory->GetParameter(3), theory->GetParameter(4),
		     theory->GetParameter(5));

  if(!bgdOnlyFit)
    {
      Double_t par[3] ={theory->GetParameter(3), theory->GetParameter(4),
			theory->GetParameter(5)};
      cout << " Nsig in [" << fXMIN << ", " << fXMAX << "] GeV = " 
       << fsig.Integral(fXMIN, fXMAX, par, 0.1)/data->GetBinWidth(0)
	   << endl;
    }

  data->SetMarkerStyle(4);
  data->GetXaxis()->SetTitle("p_{T} (GeV/c)");  
  data->SetTitle("W ' #rightarrow #mu#nu p_{T}");
  fsig.SetLineColor(kRed);
  fbgd.SetLineColor(kBlue);
  string desc = "CMS Preliminary"; string desc2 = " 100 pb^{-1}";

  TCanvas * c1 = new TCanvas();
  data->GetXaxis()->SetRangeUser(fXMIN, fXMAX);
  data->Draw("e");
  theory->Draw("same");
  fbgd.Draw("same");
  if(!bgdOnlyFit)
    fsig.Draw("same");
  
  TLatex * l1 = new TLatex(500, 15, desc.c_str());
  l1->SetTextSize(0.06); 
  TLatex * l4 = new TLatex(500, 22, desc.c_str());
  l4->SetTextSize(0.06); 
  TLatex * l3 = new TLatex(600, 15, desc2.c_str());
  l3->SetTextSize(0.05);
  TLatex * l2 = new TLatex(600, 10, desc2.c_str());
  l2->SetTextSize(0.05); 

  l1->Draw(); l2->Draw();
  c1->SetLogy();
  c1->SaveAs("muonpt_fit_log.gif");     
  c1->SetLogy(false);

  delete c1;
  c1 = new TCanvas();
  data->GetXaxis()->SetRangeUser(fXMIN, fXMAX);
  data->Draw("e");
  theory->Draw("same");
  fbgd.Draw("same");
  fsig.Draw("same");
  l4->Draw(); l3->Draw();
  c1->SaveAs("muonpt_fit_lin.gif");     
  delete c1;


  data->GetListOfFunctions()->Add(&fbgd);
  data->GetListOfFunctions()->Add(&fsig);

  TFile * fout = new TFile("fit_sigbgd.root", "recreate");
  data->Write();
  fout->Close();

}

void fitData(TH1F * data, TF1 * & theory, Results * result, int exp_no, 
	     float mass, float evt_sig, bool bgdOnlyFit)
{
  if(bgdOnlyFit)
    setNPARAM_FIT(3);
  else
    setNPARAM_FIT(6);

  TAxis *xaxis = data->GetXaxis();
  int bin_first = xaxis->FindBin(fXMIN); 
  int bin_last = xaxis->FindBin(fXMAX);
  Ntot = data->Integral(bin_first, bin_last);
  // need to check that fXMIN, fXMAX are near bin boundaries, otherwise Ntot
  // may not give correct # of signal+background events under fitting region
  // if fixNtot = true. Should add an explicit check maybe? Print on screen for now
  cout << " Fitting range: " << endl;
  cout << " first-bin # = " << bin_first << " (" 
       << xaxis->GetBinLowEdge(bin_first)
       << ", " <<  xaxis->GetBinUpEdge(bin_first) << ")" << endl;
  cout << " last-bin # = " << bin_last << " (" 
       << xaxis->GetBinLowEdge(bin_last)
       << ", " <<  xaxis->GetBinUpEdge(bin_last) << ")" << endl;

  cout << " Total # of events within bin range = " << Ntot << endl;

  setBinSize(data->GetBinWidth(0));
  setLandauBgd(isBgdLandau);
  
  if(bgdOnlyFit)
    {
      theory = new TF1("ftot", myBgd, fXMIN, fXMAX, NPARAM_FIT);
      theory->SetParameters(Ntot, fXMIN, 93.0);
    }
  else
    {
      theory = new TF1("ftot", mySigBgd, fXMIN, fXMAX, NPARAM_FIT);
      theory->SetParameters(Ntot, fXMIN, 93.0, evt_sig, mass, Fudge);
      theory->SetParName(3, "W ' evt count");
      theory->SetParName(4, "W ' Mass");
      theory->SetParName(5, "W ' Width scale factor");
    }

  theory->SetParName(0, "Total evt count");
  if(isBgdLandau)
    {
      theory->SetParName(1, "Landau mpv");
      theory->SetParName(2, "Landau sigma");
    }
  else
    {
      theory->SetParName(1, "Bgd RBW M");
      theory->SetParName(2, "Bgd RBW  #Gamma");
    }

  for(unsigned i = 0; i != NPARAM_FIT; ++i)
    theory->SetParError(i, 20);

  if(fixNtot)
    theory->FixParameter(0, Ntot);

  //  theory->FixParameter(1, fXMIN);

  if(!bgdOnlyFit)
    {
      if(use_wprime_mass_limits)
	theory->SetParLimits(4, 0, upper_wprime_mass_limit);
      
      if(fixFudge)
	theory->FixParameter(5, Fudge);
      else
	theory->SetParError(5, 0.1);
    }

  if(exp_no > -1)
    //if(exp_no == 3)
    {

      if(doFits)
	{
	  fitIt(data, theory);
	  getResults(result, theory, exp_no);
	}
  
      cout << " Nsig = " << result->Nsig << " (" << 100.*result->Nsig/Ntot 
	   << "% of total = " << Ntot << ")" << endl;
      doPlots(data, theory, bgdOnlyFit);
    }
  
  if(!custom_chi2)
    delete theory;  
}



// signal-free histogram corresponds to <mass_point_for_data>
// this is only for consistency, as no signal is added for mass_option=0
//TH1F * ref_hist[num_ref_plots] = {wp20_orig, wp08_orig, wp10_orig, 
//			    wp11_orig, wp12_orig, wp13_orig, 
//			    wp14_orig, wp15_orig, wp20_orig, 
//			    bgd_orig};
void fitSigBgd_eventLoop(unsigned mass_option, const unsigned * N_evt2gen, 
			 Results * result, TH1F ** ref, int exp_no, 
			 bool bgdOnlyFit)
{
  // bgdOnlyFit: if true (false) do background-only (signal+background) fit


  // variable defined only to make things clearer when reading the code
  const unsigned bgd_index = num_ref_plots - 1;

  string htitle = "Muon Pt"; string htitle_wp = histo_desc[mass_option];

  htitle += htitle_wp;
  makeHistogram(wprime, "wprime", htitle.c_str(), ref[mass_option]);
  htitle = "Muon Pt, Bgd-only";
  makeHistogram(bgd, "bgd", htitle.c_str(), ref[bgd_index]);

  float very_small = 0.001;
  float diff = wprime->GetBinWidth(0) - bgd->GetBinWidth(0);
  assert(TMath::Abs(diff) < very_small);

  bool store_Nsig_input = true;
  populateHistogram(wprime, N_evt2gen[mass_option], ref[mass_option],
		    store_Nsig_input);
  store_Nsig_input = false;
  populateHistogram(bgd, N_evt2gen[bgd_index], ref[bgd_index], 
		    store_Nsig_input);

  float evt_sig = 1.0*N_evt2gen[mass_option];
  float mass = -9999;
  if(mass_option == 0)
    mass = all_masses[mass_point_for_data];
  else
    mass = all_masses[mass_option];

  htitle = algo_desc_long[algo_option] + htitle_wp;
  // this is the total (signal + background) distribution 
  // (wprime and W + QCD + top + Z/DY)
  TH1F * data = 0;
  makeHistogram(data, "tot", htitle.c_str(), ref[bgd_index]);
  data->Add(bgd);
  if(mass_option > 0)
    data->Add(wprime);

  TF1 * theory = 0;
  fitData(data, theory, result, exp_no, mass, evt_sig, bgdOnlyFit);

  delete wprime; delete bgd;

}
