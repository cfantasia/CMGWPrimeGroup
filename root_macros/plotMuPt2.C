#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

#include <string>
#include <iostream>

#include "wprime_histo_constants.h"

using std::string; using std::cout; using std::endl;

// select tracking algorithm
// from 0 to Num_trkAlgos-1 (see wprime_histo_constants.h)
const unsigned tracking_option = 2;

// algorithm description
string algo; 

// histogram title;
string desc;

bool plot_without_qual_cuts = false;

string set_cuts;

void doPlots(TFile * _file0);

void plotMuPt2()
{
  assert(Num_histo_sets == 6);

  if(plot_without_qual_cuts)
    // pick histograms before quality cuts have been applied
    set_cuts = cuts_desc_short[4];  
  else
    // pick histograms after all quality cuts have been applied
    set_cuts = cuts_desc_short[5];  
    
  string input_file = "Wprime_analysis.root";
  TFile *_file0 = TFile::Open(input_file.c_str());
  if(!_file0 || _file0->IsZombie())
    {
      cout << " *** Ooops! Cannot find input file " << input_file << endl;
      return;
    }

  gStyle->SetOptStat(00000);

  algo = algo_desc_short[tracking_option];
  desc = algo_desc_long[tracking_option] + " muons";
  if(plot_without_qual_cuts)
    desc +=  ", before quality cuts";
  else
    desc +=  ", after quality cuts";

  desc += " (STARTUP, 0.84 pb^{-1})";
  doPlots(_file0);
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

void doPlots(TFile * _file0)
{

  //  gStyle->Reset();
  gStyle->SetFillColor(1);


  TCanvas * c1 = new TCanvas();
  c1->SetLogy();

  string histo = "hPT" + algo + "_" + set_cuts;
  string histoW = "W/" + histo;
  TH1F * w = (TH1F* ) _file0->Get(histoW.c_str());
  if(badHisto(w, "W"))
    return;

  string histoQ = "QCD/" + histo;
  TH1F * qcd = (TH1F* )_file0->Get(histoQ.c_str());
  if(badHisto(qcd, "QCD"))
    return;

  string histoZ = "Z/" + histo;
  TH1F * z = (TH1F* )_file0->Get(histoZ.c_str());
  if(badHisto(z, "Z"))
    return;
  
  string histoT = "Top/" + histo;
  TH1F * top = (TH1F* )_file0->Get(histoT.c_str());
  if(badHisto(top, "Top"))
    return;

  string histo10 = "wprime1.0/" + histo;
  TH1F * wp10 = (TH1F* )_file0->Get(histo10.c_str());
  if(badHisto(wp10, "wprime1.0"))
    return;

  string histo15 = "wprime1.5/" + histo;
  TH1F * wp15 = (TH1F* )_file0->Get(histo15.c_str());
  if(badHisto(wp15, "wprime1.5"))
    return;

  string histo20 = "wprime2.0/" + histo;
  TH1F * wp20 = (TH1F* )_file0->Get(histo20.c_str());
  if(badHisto(wp20, "wprime2.0"))
    return;

  string hdata = "data/" + histo;
  TH1F * data = (TH1F*) _file0->Get(hdata.c_str());
  if(badHisto(data, "data"))
    return;

  // this is the total background distribution (W + QCD + top + Z/DY)
  TH1F * bgd = new TH1F("tot_bgd", desc.c_str(), 
			top->GetNbinsX(), top->GetXaxis()->GetXmin(), 
			top->GetXaxis()->GetXmax());
  bgd->Add(z);
  bgd->Add(qcd);  
  bgd->Add(top);
  bgd->Add(w);

  string new_title = algo_desc_long[tracking_option] + " p_{T} distribution";
  data->SetTitle(new_title.c_str());
  data->SetMarkerStyle(4);
  data->SetMarkerSize(1.3);
  data->GetXaxis()->SetTitle("Muon p_{T} (GeV/c)");
  if(data->GetMinimum() < 0.00001)data->SetMinimum(0.00001);
  data->Draw("e");
  bgd->Draw("same");
 


  const int fill_style_sig = 3956;

  wp10->SetLineColor(kRed);
  wp15->SetLineColor(kGreen);
  wp20->SetLineColor(kBlue);
  wp10->SetFillColor(kRed);
  wp15->SetFillColor(kGreen);
  wp20->SetFillColor(kBlue);

  wp10->SetFillStyle(fill_style_sig);
  wp15->SetFillStyle(fill_style_sig);
  wp20->SetFillStyle(fill_style_sig);
  
   
  wp10->Draw("same");
  wp15->Draw("same");
  wp20->Draw("same");

  TLegend * lg = new TLegend(0.59, 0.67, 0.89, 0.89);
  lg->SetTextSize(0.03);
  lg->SetBorderSize(0);
  lg->SetFillColor(0);
  lg->AddEntry(bgd, "Total bgd", "F");
  lg->AddEntry(wp10, "W ' (1.0 TeV)", "F");
  lg->AddEntry(wp15, "W ' (1.5 TeV)", "F");
  lg->AddEntry(wp20, "W ' (2.0 TeV)", "F");
  lg->AddEntry(data, "data (0.84 pb^{-1})", "LP");
  lg->Draw();
  
  string file; string file2;
  if(plot_without_qual_cuts)  
    {
      file = "MuonPt_noqual.gif";
      file2 = "MuonPt_noqual.eps";
    }
  else
    {
      file = "MuonPt_qual.gif";
      file2 = "MuonPt_qual.eps";
    }

  c1->SaveAs(file.c_str());
  c1->SaveAs(file2.c_str());
  //  delete c1;

}
