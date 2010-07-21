#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TStyle.h>

#include <string>
#include <iostream>

using std::string; using std::cout; using std::endl;


string set_cuts[5] = {"all", "1mu", "iso", "jveto", "qual"};

// select tracking algorithm
// 1: global muons, 2: tracker-only, 3: tracker + 1st muon station
const unsigned tracking_option = 3;

string cut_desc[5] = {", no cuts", ", after one-muon-only requirement",
		      ", after isolation cut", 
		      ", after jet-activity veto", 
		      ", after quality cuts"};

void doPlots(unsigned i, TFile * _file0);

void plotMuPt()
{
  string input_file = "Wprime_analysis.root";
  TFile *_file0 = TFile::Open(input_file.c_str());
  if(!_file0 || _file0->IsZombie())
    {
      cout << " *** Ooops! Cannot find input file " << input_file << endl;
      return;
    }

  gStyle->SetOptStat(00000);
  for(unsigned i = 0; i != 5; ++i)
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

void doPlots(unsigned i, TFile * _file0)
{

  //  gStyle->Reset();
  gStyle->SetFillColor(1);


  string algo; string desc;
  if (tracking_option == 1)
    {
      algo = "glb"; // global muons
      desc =  "Global muons";
    }
  else if (tracking_option == 2)
    {
      algo = "trk"; // tracker-only muons
      desc =  "Tracker-only muons";
    }
  else if (tracking_option == 3)
    {
      algo = "tev"; // tracker + 1st muon-station muons
      desc =  "TPFMS muons";
    }

  desc += cut_desc[i] + " (100 pb^{-1})";

  TCanvas * c1 = new TCanvas();
  c1->SetLogy();
  THStack *hs = new THStack("hs",desc.c_str());  

  string histo = "hPT" + algo + "_" + set_cuts[i];
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
  if(hs->GetMaximum() < wp10->GetMaximum())hs->SetMaximum(wp10->GetMaximum());
  if(data->GetMaximum() < wp10->GetMaximum())data->SetMaximum(wp10->GetMaximum());
  // this is needed when background is zero in the tails and wprime (say at 2.0 TeV) is not displayed
  if(hs->GetMinimum() < 0.01)hs->SetMinimum(0.01);
  if(data->GetMinimum() < 0.01)data->SetMinimum(0.01);

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
  lg->AddEntry(data, "data (255 nb^{-1})", "LP");
  lg->Draw();

  string file = set_cuts[i] + ".gif";
  c1->SaveAs(file.c_str());
  //  delete c1;

}
