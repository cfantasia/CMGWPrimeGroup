#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLatex.h>
#include <THStack.h>
#include <TF1.h>

#include <string>
#include <iostream>

#include "UserCode/CMGWPrimeGroup/interface/TeVMuon_tracking.h"
#include "UserCode/CMGWPrimeGroup/interface/mumet_histo_constants.h"
#include "UserCode/CMGWPrimeGroup/interface/elmet_histo_constants.h"

using std::string; using std::cout; using std::endl;

// select tracking algorithm for analyses using muons
// from 0 to Num_MuTeVtrkAlgos-1 (see TeVMuon_tracking.h)
// NB: this is currently common for Muon-pt and transverse mass distributions
const unsigned tracking_option = 3;

// algorithm description (relevant for analysis using muons)
string algo; 

// histogram title;
string desc;

TFile *_file0 = 0;

// option = 1 (muon-pt), option = 2 (transverse mass)
void doPlots(int option);

float Lumi_ipb = -1;

// 1 -> Mu+MET
// 2 -> El+MET
const unsigned analysis_channel = 1;

// electron channel: if want colored histogram set "wClr" to true. 
bool wClr = true;

//const unsigned NbgdSamplesMuMET = 19;
const unsigned NbgdSamplesMuMET = 1;
const string bgdNamesMuMET[NbgdSamplesMuMET] = {"WMuNu_highPt"};
TH1F * bgdMuMET[NbgdSamplesMuMET] = {0};

const unsigned NbgdSamplesElMET = 0;
//const unsigned NbgdSamplesElMET = 27;
// color distributions for MC background
//  1      -> Wenu    : kAzure+1
//  2-4,   -> QCD     : kGray+2
//  5-7,   -> ttbar   : kMagenta
//  8-10,  -> DY      : kOrange
// 11-12,  -> Diboson : kGreen
// 13-15,  -> Wlepton : kTeal
// 16-     -> Gamma   : kRed
const string bgdNamesElMET[NbgdSamplesElMET];
#if 0
const string bgdNamesElMET[NbgdSamplesElMET] = {
  "WtoEnu_highPt",     
  "QCD_20to30_highPt", "QCD_30to80_highPt",  "QCD_80to170_highPt",
  "ttjet_highPt",      "TtoBLNus_highPt",    "TtoBLNutW_highPt", 
  "DYtoEE_highPt",     "DYtoMuMu_highPt",    "DYtoTauTau_highPt",  
  "WtoMunu_highPt",    "WtoTaunu_highPt",  
  "ZZ_highPt",         "WW_highPt",          "WZ_highPt",
  "G_0to15_highPt",    "G_15to30_highPt",    "G_30to50_highPt",    
  "G_50to80_highPt",   "G_80to120_highPt",
  "G_120to170_highPt", "G_170to300_highPt",  "G_300to470_highPt",  
  "G_470to800_highPt",
  "G_800to1400_highPt","G_1400to1800_highPt","G_1800_highPt"
};
#endif
const Color_t bgdColorElMET[7] = {kAzure+1,kGray+2,kMagenta,kOrange,kGreen,kTeal,kRed};
Color_t bgdClr;

TH1F * bgdElMET[NbgdSamplesElMET];// = {0};

int nbins_output = 125;
float xmin_output=0, xmax_output=2500;

void plotMT_interference()
{
  assert(tracking_option <= Num_MuTeVtrkAlgos-1);

  string input_file = "";
  if(analysis_channel == 1)
    input_file = "Wprime_analysis_MuMET.root";
  else  if(analysis_channel == 2)
    input_file = "Wprime_analysis_ElMET.root";
  else
    abort();

  _file0 = TFile::Open(input_file.c_str());
  if(!_file0 || _file0->IsZombie())
    {
      cout << " *** Ooops! Cannot find input file " << input_file << endl;
      return;
    }

  string histo_lumi = "lumi_ipb";
  TH1F * lumi_ipb = (TH1F* )_file0->Get(histo_lumi.c_str());
  Lumi_ipb = lumi_ipb->GetBinContent(1);

  gStyle->SetOptStat(00000);

  if(analysis_channel == 1)
    algo = algo_desc_short[tracking_option];
  else if(analysis_channel == 2)
    algo = "";

  // doPlots(1); // lepton-pt
  doPlots(2); // lepton + MET transverse mass
}

bool badHisto(TH1F * h, string s)
{
  if(!h)
    {
      cout << " *** Oops! Didn't find " << s << " histogram... " << endl;
      return true;
    }

  return false;
}

// option = 1 (muon-pt), option = 2 (transverse mass)
void doPlots(int option)
{
  //gStyle->SetFillColor(1);

  // interested in plot made after all cuts have been applied
  string final_histo_desc;
  if(analysis_channel == 1)
    final_histo_desc = mumet_cuts_desc_short[Num_mumet_cuts-1];
  else if(analysis_channel == 2)
    final_histo_desc = elmet_cuts_desc_short[Num_elmet_cuts-1];

  string prefix = "";
  if(option == 1 && analysis_channel == 1)
    prefix = "hPT";
  if(option == 1 && analysis_channel == 2)
    prefix = "hET";
  else if(option == 2)
    prefix = "hTM";
 
  TH1F ** bgdSamples = 0;
  const string * bgdNames = 0;
  unsigned NbgdSamples = 0;

  if(analysis_channel == 1)
    {
      bgdSamples = &bgdMuMET[0];
      bgdNames = &bgdNamesMuMET[0];
      NbgdSamples = NbgdSamplesMuMET;
    }
  else if(analysis_channel == 2)
    {
      bgdSamples = &bgdElMET[0];
      bgdNames = &bgdNamesElMET[0];
      NbgdSamples = NbgdSamplesElMET;
    }

  string histo = prefix + algo + "_" + final_histo_desc;  
  for(unsigned i = 0; i != NbgdSamples; ++i)
    {
      string histo_i = bgdNames[i] + "/" + histo;
      bgdSamples[i] = (TH1F* ) _file0->Get(histo_i.c_str());
      if(badHisto(bgdSamples[i], histo_i))
	return;
    }

  string mass = "1.5";//W' mass

  string histo_wp = "wprime" + mass + "/" + histo;
  TH1F * wp = (TH1F* )_file0->Get(histo_wp.c_str());
  if(badHisto(wp, "wprime" + mass))
    return;

  string histo_wp_noint = "wprime" + mass + "_noint/" + histo;
  TH1F * wp_noint = (TH1F* )_file0->Get(histo_wp_noint.c_str());
  if(badHisto(wp_noint, "wprime" + mass + "_noint"))
    return;

  string histo_wp_oppsign = "wprime" + mass + "_oppsign/" + histo;
  TH1F * wp_oppsign = (TH1F* )_file0->Get(histo_wp_oppsign.c_str());
  if(badHisto(wp_oppsign, "wprime" + mass + "_oppsign"))
    return;

  string histo_wp_samesign = "wprime" + mass + "_samesign/" + histo;
  TH1F * wp_samesign = (TH1F* )_file0->Get(histo_wp_samesign.c_str());
  if(badHisto(wp_samesign, "wprime" + mass + "_samesign"))
    return;

  string hname = "tot_bgd";
  double xmin = -1; double xmax = -1; 
  double xmax_ratio = -1; double xmax_cumu = -1;
  string title = "INVALID";
  string var_plotted = "INVALID2";
  char data_ipb[1024]; sprintf(data_ipb, " data (%4.1f pb^{-1})", Lumi_ipb);
  char lumi_value[1024]; sprintf(lumi_value, "%4.1fipb", Lumi_ipb);
  char lumi_value2[1024]; sprintf(lumi_value2, "%4.1f pb^{-1}", Lumi_ipb);
  string cms_prelim = "CMS Preliminary 2011"; 
  char lumi_sqrts_[1024]; 
  sprintf(lumi_sqrts_, "L_{int} = %4.1f pb^{-1}, #sqrt{s} = 7 TeV", Lumi_ipb);
  string lumi_sqrts = lumi_sqrts_;
  float x_offset = 0;

  int Nbins = bgdSamples[0]->GetNbinsX();

  if(option == 1)
    {
      hname += "_mupt";
      desc = " p_{T} distribution";
      xmin = 150; xmax = 1200; 
      title = "Muon p_{T} (GeV/c)";
      var_plotted = "MuPt";
      x_offset = -100;
    }
  else if(option == 2)
    {
      hname += "_TM";
      string desc0 = "";
      if(analysis_channel == 1)
	desc0 = "#mu";
      else if(analysis_channel == 2)
	desc0 = "e";
      desc = desc0 + "&ME_{T} transverse mass: 2011 data (" 
	+ string(lumi_value2) + ")";
      xmin = 200; xmax = 2000; 
      title = "M_{T} (GeV/c^{2})";
      var_plotted = "TM";
    }
  string file = var_plotted + "_" + lumi_value + ".gif";

  // ============== CREATE DISTRIBUTIONS HERE ======================
  // this is the total background distribution (W + QCD + top + Z/DY)
  TH1F * bgd = new TH1F(hname.c_str(), desc.c_str(), 
			Nbins, bgdSamples[0]->GetXaxis()->GetXmin(), 
			bgdSamples[0]->GetXaxis()->GetXmax());
  TH1F * bgd_output = new TH1F("MT_bgd","MT_bgd",nbins_output,xmin_output,xmax_output);
  THStack *hsbgd =new THStack(hname.c_str(),desc.c_str());//+++++++++
  for(unsigned i = 0; i != NbgdSamples; ++i){
    //for(int i = NbgdSamples - 1; i != -1; i--){
    if(analysis_channel == 2){
      //  1, 2-4, 5-7, 8-10,11-12,13-15,16-
      if(i<1){bgdClr=bgdColorElMET[0];}
      else if(i>1 && i<5){bgdClr=bgdColorElMET[1];}
      else if(4<i && i<8){bgdClr=bgdColorElMET[2];}
      else if(7<i && i<11){bgdClr=bgdColorElMET[3];}
      else if(10<i && i<13){bgdClr=bgdColorElMET[4];}
      else if(12<i && i<16){bgdClr=bgdColorElMET[5];}
      else if(i>15){bgdClr=bgdColorElMET[6];}

      bgdSamples[i]->SetLineColor(bgdClr);bgdSamples[i]->SetFillColor(bgdClr);
      hsbgd->Add(bgdSamples[i]);//+++++++++
    }
    bgd->Add(bgdSamples[i]);
    bgd_output->Add(bgdSamples[i]);
  }
  
  const int fill_style_sig = 3001;
  const int fill_style_bgd = 3001;

  bgd->SetLineColor(kAzure+1);
  bgd->SetFillColor(kAzure+1);
  bgd->SetFillStyle(fill_style_bgd);

  // =============== PLOT DISTRIBUTIONS HERE ========================
  
  TCanvas * c1 = new TCanvas();
  c1->SetLogy();

  if(wClr && analysis_channel == 2)hsbgd->Draw("e");
  else bgd->Draw("e");

  wp->SetLineColor(kRed);
  wp_noint->SetLineColor(kRed+1);
  wp_oppsign->SetLineColor(kRed+2);
  wp_samesign->SetLineColor(kRed+3);
  wp->SetFillColor(kRed);
  wp_noint->SetFillColor(kRed+1);
  wp_oppsign->SetFillColor(kRed+2);
  wp_samesign->SetFillColor(kRed+3);
  wp->SetFillStyle(fill_style_sig);
  wp_noint->SetFillStyle(fill_style_sig);
  wp_oppsign->SetFillStyle(fill_style_sig);   
  wp_samesign->SetFillStyle(fill_style_sig);
  wp->Draw("same");
  wp_noint->Draw("same");
  wp_oppsign->Draw("same");
  wp_samesign->Draw("same");
  
  TLegend * lg = new TLegend(0.52, 0.67, 0.82, 0.89);
  lg->SetTextSize(0.03);
  lg->SetBorderSize(0);
  lg->SetFillColor(0);
  lg->AddEntry(bgd, "Total bgd", "F");
  string histo_entry_title = "W ' (" + mass + " TeV)";
  lg->AddEntry(wp, histo_entry_title.c_str(), "F");
  histo_entry_title = "W ' (" + mass + " TeV), no interference";
  lg->AddEntry(wp_noint, histo_entry_title.c_str(), "F");
  histo_entry_title = "W ' (" + mass + " TeV), negative interference";
  lg->AddEntry(wp_oppsign, histo_entry_title.c_str(), "F");
  histo_entry_title = "W ' (" + mass + " TeV), positive interference";
  lg->AddEntry(wp_samesign, histo_entry_title.c_str(), "F");
  lg->Draw();
  
  c1->SaveAs(file.c_str());

  TFile *bgd_output_file = new TFile("Wprime_bgd_mu.root","recreate");
  bgd_output->Write();
  bgd_output_file->Close();  
}
