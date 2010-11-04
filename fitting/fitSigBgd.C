#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TRandom.h>
#include <TMinuit.h>
#include <TTree.h>
#include <TVirtualFitter.h>
#include <TBranch.h>

#include <string>
#include <iostream>
#include <sstream>

#include "fit_wprime.h"
#include "common_fit.h"
#include "Results.h"
#include "../root_macros/wprime_histo_constants.h"

using std::string; using std::cout; using std::endl;

// see common_fit.h for the mass points to which this array corresponds
unsigned N_evt2gen[num_ref_plots]; 

void getInputHistograms();
void makeReferenceHistograms();

extern void fitSigBgd_eventLoop(unsigned mass_option, const unsigned * N_evt2gen,
				Results * result, TH1F ** ref, int exp_no,
				bool bgdOnlyFit);

extern void setLumi_ipb(float lumi);


TFile * file0 = 0;
TFile * gfile = 0;

// these are the original histograms (full statistics)
// from these, will draw random distributions according for given luminosity
TH1F * w_orig = 0;
TH1F * qcd_orig = 0;
TH1F * z_orig = 0;
TH1F * top_orig = 0;
TH1F * wp08_orig = 0;
TH1F * wp10_orig = 0;
TH1F * wp11_orig = 0;
TH1F * wp12_orig = 0;
TH1F * wp13_orig = 0;
TH1F * wp14_orig = 0;
TH1F * wp15_orig = 0;
TH1F * wp20_orig = 0;
TH1F * bgd_orig = 0;

TH1F * lumi_ipb = 0;

// resolution function
TH1F * g0[mass_points] = {0};

bool badHisto(TH1F * h, string s)
{
  if(!h)
    {
      cout << " *** Oops! Didn't find " << s << " histogram... " << endl;
      return true;
    }

  return false;
}

void cleanup(Results * result, TFile * output_file)
{
  delete result;	
  delete bgd_orig;
  delete output_file;
}

// see common_fit.h for the mass points to which mass_option corresponds
// N_EXP: # of pseudo-experiments
int fitSigBgd(unsigned mass_option, unsigned N_EXP, bool bgdOnlyFit)
{
  if(mass_option > 11)
    {
      cout << " Ooops! Only understand " << 
	" option = 0 (no signal) through 11 (2.0 TeV) "<< endl;
      abort();
    }

  gStyle->SetFillColor(1);
  
  string input_file = "Wprime_analysis.root";
  string resolution_file = "gsmear.root";

  file0 = TFile::Open(input_file.c_str());
  if(!file0 || file0->IsZombie())
    {
      cout << " *** Ooops! Cannot find input file " << input_file << endl;
      return -1;
    }

  gfile = TFile::Open(resolution_file.c_str());
  if(!gfile || gfile->IsZombie())
    {
      cout << " *** Ooops! Cannot find muon-pt resolution file " 
	   << resolution_file << endl;
      return -1;
    }


  string bgdOnly_suffix;
  if(bgdOnlyFit)
    bgdOnly_suffix = "H0";
  else
    bgdOnly_suffix = "H1";

  string file_name = "output_" + data_desc[mass_option] + "_" + 
    bgdOnly_suffix + ".root";
  string tree_name = "Fit results for " + histo_desc[mass_option] + 
    " wprime experiments";

  TFile * output_file = new TFile(file_name.c_str(), "recreate");
  TTree * root_tree = new TTree("tree", tree_name.c_str());
  root_tree->SetDirectory(output_file);
  //  TVirtualFitter::SetDefaultFitter("Fumili2");
  TVirtualFitter::SetDefaultFitter("Minuit");
  
  Results * result = new Results();
  root_tree->Branch("myTuple", "Results", &result, 8000, 2);
  
  getInputHistograms();
  makeReferenceHistograms();
  setResolution(g0[mass_option]);

  TH1F * ref_hist[num_ref_plots] = {0, wp08_orig, wp10_orig, 
				    wp11_orig, wp12_orig, wp13_orig, 
				    wp14_orig, wp15_orig, wp20_orig, 
				    bgd_orig};
  // signal-free histogram corresponds to <mass_point_for_data> (see common_fit.h)
  // this is only for consistency, as no signal is added for mass_option=0
  ref_hist[0] = ref_hist[mass_point_for_data];

  for(unsigned exp_no = 1; exp_no <= N_EXP; ++exp_no)
    { // loop over pseudo-experiments
      cout<< " ************************************************************* "
	  << endl;
      cout << " PSEUDO-EXPERIMENT # " << exp_no << endl;
      cout<< " ************************************************************* "
	  << endl;
      fitSigBgd_eventLoop(mass_option, N_evt2gen, result, ref_hist, 
			  exp_no, bgdOnlyFit);
      root_tree->Fill();
      root_tree->AutoSave();
    } // loop over pseudo-experiments
      
  output_file->Close();
  cleanup(result, output_file);

  file0->Close(); gfile->Close();
  return 0;
}

void getInputHistograms()
{
  string histo_lumi = "lumi_ipb";
  lumi_ipb = (TH1F* )file0->Get(histo_lumi.c_str());
  if(badHisto(lumi_ipb, histo_lumi))
    return;

  string histo = "hPT" + algo_desc_short[algo_option] + "_"
    + final_histo_desc;

  string histoW = "W/" + histo;
  w_orig = (TH1F* ) file0->Get(histoW.c_str());
  if(badHisto(w_orig, "W"))
    return;
  
  string histoQ = "QCD/" + histo;
  qcd_orig = (TH1F* )file0->Get(histoQ.c_str());
  if(badHisto(qcd_orig, "QCD"))
    return;
  
  string histoZ = "Z/" + histo;
  z_orig = (TH1F* )file0->Get(histoZ.c_str());
  if(badHisto(z_orig, "Z"))
    return;
  
  string histoT = "Top/" + histo;
  top_orig = (TH1F* )file0->Get(histoT.c_str());
  if(badHisto(top_orig, "Top"))
    return;

  string histo08 = "wprime0.8/" + histo;
  wp08_orig = (TH1F* )file0->Get(histo08.c_str());
  if(badHisto(wp08_orig, "wprime0.8"))
    return;

  string histo10 = "wprime1.0/" + histo;
  wp10_orig = (TH1F* )file0->Get(histo10.c_str());
  if(badHisto(wp10_orig, "wprime1.0"))
    return;

  string histo11 = "wprime1.1/" + histo;
  wp11_orig = (TH1F* )file0->Get(histo11.c_str());
  if(badHisto(wp11_orig, "wprime1.1"))
    return;

  string histo12 = "wprime1.2/" + histo;
  wp12_orig = (TH1F* )file0->Get(histo12.c_str());
  if(badHisto(wp12_orig, "wprime1.2"))
    return;

  string histo13 = "wprime1.3/" + histo;
  wp13_orig = (TH1F* )file0->Get(histo13.c_str());
  if(badHisto(wp13_orig, "wprime1.3"))
    return;

  string histo14 = "wprime1.4/" + histo;
  wp14_orig = (TH1F* )file0->Get(histo14.c_str());
  if(badHisto(wp14_orig, "wprime1.4"))
    return;

  string histo15 = "wprime1.5/" + histo;
  wp15_orig = (TH1F* )file0->Get(histo15.c_str());
  if(badHisto(wp15_orig, "wprime1.5"))
    return;
  
  string histo20 = "wprime2.0/" + histo;
  wp20_orig = (TH1F* )file0->Get(histo20.c_str());
  if(badHisto(wp20_orig, "wprime2.0"))
    return;

  string suffix = algo_desc_short[algo_option];

  string gname;

  // i=0 corresponds to signal-free distribution
  for(unsigned i = 0; i < mass_points; ++i)
    {
      string tmp = gname_pre[i];
      if (i == 0) tmp = gname_pre[mass_point_for_data];

      gname = tmp + "_" + suffix;
      g0[i] = (TH1F* )gfile->Get(gname.c_str());
      if(badHisto(g0[i], gname))
	return;
    }
  
}

void makeReferenceHistograms()
{
  string histo_name = "bgd" + algo_desc_short[algo_option];
  string histo_title = "Bgd-only" + algo_desc_long[algo_option]; 
  bgd_orig = new TH1F(histo_name.c_str(), histo_title.c_str(), 
		 top_orig->GetNbinsX(), top_orig->GetXaxis()->GetXmin(), 
		 top_orig->GetXaxis()->GetXmax());
  bgd_orig->SetDirectory(0);
  bgd_orig->Add(z_orig);
  bgd_orig->Add(qcd_orig);  
  bgd_orig->Add(top_orig);
  bgd_orig->Add(w_orig);
  bgd_orig->GetXaxis()->SetTitle("Muon p_{T} (GeV/c)");

  // # of events to generate per sample
  assert(num_ref_plots == 10);

  setLumi_ipb(lumi_ipb->GetBinContent(1));
// use to scale histograms in Wprime_analysis.root
  float k = integ_lumi/lumi_ipb->GetBinContent(1); 
  cout << " MC distributions made with " << lumi_ipb->GetBinContent(1)
       << " pb^{-1}, scale factor = " << k << endl;
  N_evt2gen[0] = 1;
  N_evt2gen[1] = int(wp08_orig->Integral()*k + 0.5);
  N_evt2gen[2] = int(wp10_orig->Integral()*k + 0.5);
  N_evt2gen[3] = int(wp11_orig->Integral()*k + 0.5);
  N_evt2gen[4] = int(wp12_orig->Integral()*k + 0.5);
  N_evt2gen[5] = int(wp13_orig->Integral()*k + 0.5);
  N_evt2gen[6] = int(wp14_orig->Integral()*k + 0.5);
  N_evt2gen[7] = int(wp15_orig->Integral()*k + 0.5);
  N_evt2gen[8] = int(wp20_orig->Integral()*k + 0.5);
  N_evt2gen[9] = int((w_orig->Integral() + qcd_orig->Integral()
		      + z_orig->Integral() + top_orig->Integral())*k + 0.5);
  
  float xmin = wp10_orig->GetXaxis()->GetXmin();
  float xmax = wp10_orig->GetXaxis()->GetXmax();
  cout << "\n\n Total # of events expected after all cuts for " << integ_lumi 
       << " pb^-1 between " << xmin << " and " << xmax << " GeV" << endl;
  // # of events expected after cuts for various wprime mass points
  for(unsigned i = 0; i != num_ref_plots; ++i)
    cout << histo_desc[i] << " : " << N_evt2gen[i] << endl;

}
