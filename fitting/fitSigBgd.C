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

using std::string; using std::cout; using std::endl;

unsigned N_evt2gen[num_ref_plots]; // W'(1.0), W'(1.5), W'(2.0), Bgd

// this corresponds to the histogram after all cuts have been applied
string set_cuts = "qual";

string algo[N_algos] = {"glb", "trk", "tev"};

void getInputHistograms();
void makeReferenceHistograms();

extern void fitSigBgd_eventLoop(unsigned mass_option, const unsigned * N_evt2gen,
				Results * result, TH1F ** ref, int exp_no);

TFile * file0 = 0;
TFile * gfile = 0;

// these are the original histograms (full statistics)
// from these, will draw random distributions according for given luminosity
TH1F * w_orig = 0;
TH1F * qcd_orig = 0;
TH1F * z_orig = 0;
TH1F * top_orig = 0;
TH1F * wp10_orig = 0;
TH1F * wp15_orig = 0;
TH1F * wp20_orig = 0;
TH1F * bgd_orig = 0;

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

// mass_option 1: 1.0 TeV, 2: 1.5 TeV
// N_EXP: # of pseudo-experiments
int fitSigBgd(unsigned mass_option, unsigned N_EXP)
{
  if(mass_option != 1 && mass_option != 2 && mass_option != 3)
    {
      cout << " Ooops! Only understand " << 
	" option=1 (1.0 TeV) or 2 (1.5 TeV) or 3 (2.0 TeV)"<< endl;
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

  string file_name; string tree_name = "Fit results for ";
  switch(mass_option)
    {
    case 1:
      file_name = "output_1.0TeV.root";
      tree_name += "1.0 TeV wprime experiments";
      break;
    case 2:
      file_name = "output_1.5TeV.root";
      tree_name += "1.5 TeV wprime experiments";
      break;
    case 3:
      file_name = "output_2.0TeV.root";
      tree_name += "2.0 TeV wprime experiments";
      break;
    default:
      ; // do nothing
    }

  TFile * output_file = new TFile(file_name.c_str(), "recreate");
  TTree * root_tree = new TTree("tree", tree_name.c_str());
  root_tree->SetDirectory(output_file);
  //  TVirtualFitter::SetDefaultFitter("Fumili2");
  TVirtualFitter::SetDefaultFitter("Minuit");
  
  Results * result = new Results();
  root_tree->Branch("myTuple", "Results", &result, 8000, 2);
  
  //  for(unsigned i = 0; i != N_algos; ++i)
  getInputHistograms();
  makeReferenceHistograms();
  setResolution(g0[mass_option - 1]);
  TH1F * ref_hist[num_ref_plots] = {wp10_orig, wp15_orig, wp20_orig, 
				    bgd_orig};

  for(unsigned j = 1; j <= N_EXP; ++j)
    { // loop over pseudo-experiments
      cout<< " ************************************************************* "
	  << endl;
      cout << " PSEUDO-EXPERIMENT # " << j << endl;
      cout<< " ************************************************************* "
	  << endl;
      fitSigBgd_eventLoop(mass_option, N_evt2gen, result, ref_hist, j);
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
  string histo = "hPT" + algo[algo_option] + "_" + set_cuts;

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

  string histo10 = "wprime10/" + histo;
  wp10_orig = (TH1F* )file0->Get(histo10.c_str());
  if(badHisto(wp10_orig, "wprime10"))
    return;
  
  string histo15 = "wprime15/" + histo;
  wp15_orig = (TH1F* )file0->Get(histo15.c_str());
  if(badHisto(wp15_orig, "wprime15"))
    return;
  
  string histo20 = "wprime20/" + histo;
  wp20_orig = (TH1F* )file0->Get(histo20.c_str());
  if(badHisto(wp20_orig, "wprime20"))
    return;

  g0[0] = (TH1F* )gfile->Get("g10");
  if(badHisto(g0[0], "g10"))
    return;
  g0[1] = (TH1F* )gfile->Get("g15");
  if(badHisto(g0[1], "g15"))
    return;
  g0[2] = (TH1F* )gfile->Get("g20");
  if(badHisto(g0[2], "g20"))
    return;
}

void makeReferenceHistograms()
{
  string histo_name = "bgd" + algo[algo_option];
  string histo_title = "Bgd-only" + desc[algo_option]; 
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
  assert(num_ref_plots == 4);
  float k = integ_lumi/100;
  N_evt2gen[0] = int(wp10_orig->Integral()*k + 0.5);
  N_evt2gen[1] = int(wp15_orig->Integral()*k + 0.5);
  N_evt2gen[2] = int(wp20_orig->Integral()*k + 0.5);
  N_evt2gen[3] = int((w_orig->Integral() + qcd_orig->Integral()
		      + z_orig->Integral() + top_orig->Integral())*k + 0.5);
  
  float xmin = wp10_orig->GetXaxis()->GetXmin();
  float xmax = wp10_orig->GetXaxis()->GetXmax();
  cout << "\n\n Total # of events expected after all cuts for " << integ_lumi 
       << " pb^-1 between " << xmin << " and " << xmax << " GeV" << endl;
  // # of events expected after cuts for wprime with mass 1, 1.5, 2 TeV and bgd
  for(unsigned i = 0; i != num_ref_plots; ++i)
    cout << histo_desc[i] << " : " << N_evt2gen[i] << endl;

}
