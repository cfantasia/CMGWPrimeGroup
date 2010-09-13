#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TStyle.h>

#include <string>
#include <iostream>

#include "wprime_histo_constants.h"

using std::string; using std::cout; using std::endl;

// select tracking algorithm
// from 0 to Num_trkAlgos-1 (see wprime_histo_constants.h)
const unsigned tracking_option = 2;

void doPlots(unsigned i, TFile * _file0);

void plotMuPt()
{
  assert(tracking_option >=0 && tracking_option <= Num_trkAlgos-1);

  string input_file = "Wprime_analysis.root";
  TFile *_file0 = TFile::Open(input_file.c_str());
  if(!_file0 || _file0->IsZombie())
    {
      cout << " *** Ooops! Cannot find input file " << input_file << endl;
      return;
    }

  gStyle->SetOptStat(00000);
  for(int i = 0; i != Num_histo_sets; ++i)
    doPlots(i, _file0);
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

// i must be between 0 and Num_histo_sets-1
void doPlots(unsigned i, TFile * _file0)
{

  //  gStyle->Reset();
  gStyle->SetFillColor(1);


  string algo = algo_desc_short[tracking_option]; 
  string desc = algo_desc_long[tracking_option] + " muons, " + 
    cuts_desc_long[i] + " (2.88 pb^{-1})";

  TCanvas * c1 = new TCanvas();
  c1->SetLogy();
  THStack *hs = new THStack("hs",desc.c_str());  

  string histo = "hPT" + algo + "_" + cuts_desc_short[i];
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

  top->SetLineColor(kMagenta);
  z->SetLineColor(kCyan);
  qcd->SetLineColor(kBlack);
  w->SetLineColor(kOrange);
  wp10->SetLineColor(kRed);
  wp15->SetLineColor(kGreen);
  wp20->SetLineColor(kBlue);
  

  top->SetFillColor(kMagenta);
  z->SetFillColor(kCyan);
  qcd->SetFillColor(kBlack);
  w->SetFillColor(kOrange);
  wp10->SetFillColor(kRed);
  wp15->SetFillColor(kGreen);
  wp20->SetFillColor(kBlue);

  const int fill_style_bgd = 3544;
  const int fill_style_sig = 3956;

  top->SetFillStyle(fill_style_bgd);
  z->SetFillStyle(fill_style_bgd);
  qcd->SetFillStyle(fill_style_bgd);
  w->SetFillStyle(fill_style_bgd);
  
  wp10->SetFillStyle(fill_style_sig);
  wp15->SetFillStyle(fill_style_sig);
  wp20->SetFillStyle(fill_style_sig);
  
  hs->Add(z);
  hs->Add(qcd);  
  hs->Add(top);
  hs->Add(w);

  // this is needed when background is eliminated and the y-axis is linear (as opposed to log)
  if(data->GetMaximum() < wp10->GetMaximum())data->SetMaximum(wp10->GetMaximum());
  // this is needed when background is zero in the tails and wprime (say at 2.0 TeV) is not displayed
  if(data->GetMinimum() < 0.00001)data->SetMinimum(0.00001);

  if (i == Num_histo_sets-1)
    {
      string new_title = algo_desc_long[tracking_option] + " p_{T} distribution";
      data->SetTitle(new_title.c_str());
    }
  data->SetMarkerStyle(4);
  data->SetMarkerSize(1.3);
  data->GetXaxis()->SetTitle("Muon p_{T} (GeV/c)");
  data->Draw("e");
  hs->Draw("same");
  
  wp10->Draw("same");
  wp15->Draw("same");
  wp20->Draw("same");
  
  TLegend * lg = new TLegend(0.59, 0.67, 0.89, 0.89);
  lg->SetTextSize(0.03);
  lg->SetBorderSize(0);
  lg->SetFillColor(0);
  lg->AddEntry(z, "Drell Yan (Z/Z*)", "F");
  lg->AddEntry(qcd, "QCD", "F");
  lg->AddEntry(top, "Top", "F");
  lg->AddEntry(w, "W/W*", "F");
  lg->AddEntry(wp10, "W ' (1.0 TeV)", "F");
  lg->AddEntry(wp15, "W ' (1.5 TeV)", "F");
  lg->AddEntry(wp20, "W ' (2.0 TeV)", "F");
  lg->AddEntry(data, "data (2.88 pb^{-1})", "LP");
  lg->Draw();

  string file = cuts_desc_short[i] + ".gif";
  c1->SaveAs(file.c_str());
  //  delete c1;

}
