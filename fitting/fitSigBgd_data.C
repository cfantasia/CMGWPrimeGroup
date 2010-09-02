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

void getInputHistograms();

extern void fitData(TH1F * data, TF1 * & theory, Results * result, int exp_no, float mass, float evt_sig, bool bgdOnlyFit);

TFile * file0 = 0;
TFile * gfile = 0;

// resolution function
TH1F * g0 = 0;

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
  delete output_file;
}

TH1F * data = 0;
TF1 * theory = 0;

int fitSigBgd_data()
{
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

  string file_name = "output_Wprime_7TeV.root"; 
  string tree_name = "Fit results for W->mu pt spectrum";

  TFile * output_file = new TFile(file_name.c_str(), "recreate");
  TTree * root_tree = new TTree("tree", tree_name.c_str());
  root_tree->SetDirectory(output_file);
  //  TVirtualFitter::SetDefaultFitter("Fumili2");
  TVirtualFitter::SetDefaultFitter("Minuit");
  
  Results * result = new Results();
  root_tree->Branch("myTuple", "Results", &result, 8000, 2);
  
  getInputHistograms();
  setResolution(g0);

  float mass = 1100; float evt_sig = 1.0; bool bgdOnlyFit = false;
  fitData(data, theory, result, 0, mass, evt_sig, bgdOnlyFit);
  root_tree->Fill();
  root_tree->AutoSave();
      
  output_file->Close();
  cleanup(result, output_file);

  file0->Close(); gfile->Close();
  return 0;
}

void getInputHistograms()
{
  string histo = "data/hPT" + algo_desc_short[algo_option] + "_" + 
    final_histo_desc;

  data = (TH1F* ) file0->Get(histo.c_str());
  if(badHisto(data, "data"))
    return;

  string suffix = algo_desc_short[algo_option];
  // need to specify what wprime mass point to be used for resolution function
  string gname = "g11_" + suffix;
  g0 = (TH1F* )gfile->Get(gname.c_str());
  if(badHisto(g0, gname.c_str()))
    return;
}

