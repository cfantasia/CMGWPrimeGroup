#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLatex.h>

#include <string>
#include <iostream>

#include "UserCode/CMGWPrimeGroup/interface/TeVMuon_tracking.h"
#include "UserCode/CMGWPrimeGroup/interface/mumet_histo_constants.h"

using std::string; using std::cout; using std::endl;

// select tracking algorithm
// from 0 to Num_MuTeVtrkAlgos-1 (see TeVMuon_tracking.h)
// NB: this is currently common for Muon-pt and transverse mass distributions
const unsigned tracking_option = 3;

// algorithm description
string algo; 

// histogram title;
string desc;

TFile *_file0 = 0;

// option = 1 (muon-pt), option = 2 (transverse mass)
void doPlots(int option);

float Lumi_ipb = -1;

void plotMT()
{
  assert(tracking_option <= Num_MuTeVtrkAlgos-1);

  string input_file = "Wprime_analysis.root";
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

  algo = algo_desc_short[tracking_option];

  //doPlots(1); // muon-pt
  doPlots(2); // transverse mass
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
  string final_histo_desc = mumet_cuts_desc_short[Num_mumet_cuts-1];

  string prefix = "";
  if(option == 1)
    prefix = "hPT";
  else if(option == 2)
    prefix = "hTM";


  const int NbgdSamples = 17;
  TH1F * bgd_samples[NbgdSamples] = {0};
  string bgdNames[NbgdSamples] = {
    "DYmumu_lowPt", "DYmumu_highPt", "DYtautau_lowPt", "DYtautau_lowPt",
    "QCD_lowPt", "WMinusMu_lowPt", "WMinusMu_highPt", "WPlusMu_lowPt",
    "WPlusMu_highPt", "WPlusTau_lowPt", "WPlusTau_highPt", "WW_lowPt",
    "WW_highPt", "WZ_lowPt", "WZ_highPt", "ZZ_lowPt", "ZZ_highPt"};

  string histo = prefix + algo + "_" + final_histo_desc;  
  for(int i = 0; i != NbgdSamples; ++i)
    {
      string histo_i = bgdNames[i] + "/" + histo;
      bgd_samples[i] = (TH1F* ) _file0->Get(histo_i.c_str());
      if(badHisto(bgd_samples[i], bgdNames[i]))
	return;
    }

#if 0
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
#endif

  string hdata = "data/" + histo;
  TH1F * data = (TH1F*) _file0->Get(hdata.c_str());
  if(badHisto(data, "data"))
    return;

  string hdata2 = "data2/" + histo;
  TH1F * data2 = (TH1F*) _file0->Get(hdata2.c_str());
  if(badHisto(data2, "data2"))
    return;

  data->Add(data2);

  string hname = "tot_bgd";
  string hname_ratio = "ratio";
  string hname_cumu="bgd_cumu";
  string hname_cumu_data = "data_cumu";
  double xmin = -1; double xmax = -1; 
  double xmax_ratio = -1; double xmax_cumu = -1;
  string title = "INVALID";
  string var_plotted = "INVALID2";
  char data_ipb[1024]; sprintf(data_ipb, " data (%4.2f pb^{-1})", Lumi_ipb);
  char lumi_value[1024]; sprintf(lumi_value, "%4.2fipb", Lumi_ipb);
  char lumi_value2[1024]; sprintf(lumi_value2, "%4.2f pb^{-1}", Lumi_ipb);
  string cms_prelim = "CMS Preliminary 2011"; 
  char lumi_sqrts_[1024]; 
  sprintf(lumi_sqrts_, "L_{int} = %4.2f pb^{-1}, #sqrt{s} = 7 TeV", Lumi_ipb);
  string lumi_sqrts = lumi_sqrts_;
  float x_offset = 0;

  int Nbins = data->GetNbinsX();

  if(option == 1)
    {
      hname += "_mupt";
      hname_ratio += "_mupt";
      desc = " p_{T} distribution";
      xmin = 25; xmax = 600; xmax_cumu = xmax_ratio = 600;
      title = "Muon p_{T} (GeV/c)";
      var_plotted = "MuPt";
      x_offset = -100;
    }
  else if(option == 2)
    {
      hname += "_TM";
      hname_ratio += "_TM";
      desc = "#mu & ME_{T} transverse mass: 2011 data (" 
	+ string(lumi_value2) + ")";
      xmin = 50; xmax = 900; xmax_cumu = 800; xmax_ratio = 800;
      title = "M_{T} (GeV/c^{2})";
      var_plotted = "TM";
    }
  string ratio_desc = "Data-to-MC ratio of " + desc;
  string cumu_desc = "Cumulative" + desc;
  string file = var_plotted + "_" + lumi_value + ".gif";
  string file_cumu = var_plotted + "_cumu_" + lumi_value + ".gif"; 
  string file_ratio = var_plotted + "_ratio_" + lumi_value + ".gif"; 

  // ============== CREATE DISTRIBUTIONS HERE ======================
  // this is the total background distribution (W + QCD + top + Z/DY)
  TH1F * bgd = new TH1F(hname.c_str(), desc.c_str(), 
			Nbins, data->GetXaxis()->GetXmin(), 
			data->GetXaxis()->GetXmax());
  for(int i = 0; i != NbgdSamples; ++i)
    bgd->Add(bgd_samples[i]);

  data->SetTitle(desc.c_str());
  data->SetMarkerStyle(8);

  // this is the data-over-MC ratio
  TH1F * ratio = new TH1F(hname_ratio.c_str(), ratio_desc.c_str(),
			  Nbins, data->GetXaxis()->GetXmin(), 
			  data->GetXaxis()->GetXmax());
  // this is the total cumulative bgd distribution (W + QCD + top + Z/DY)
  TH1F * bgd_cumu = new TH1F(hname_cumu.c_str(), cumu_desc.c_str(), 
			    Nbins, data->GetXaxis()->GetXmin(), 
			    data->GetXaxis()->GetXmax());
  TH1F* data_cumu = new TH1F(hname_cumu_data.c_str(), cumu_desc.c_str(), 
			    Nbins, data->GetXaxis()->GetXmin(), 
			    data->GetXaxis()->GetXmax());
  
  for(int i = 1; i <= Nbins; ++i)
    {
      bgd_cumu->SetBinContent(i, bgd->Integral(i, Nbins+1));
      data_cumu->SetBinContent(i, data->Integral(i, Nbins+1));
    }


  data->GetXaxis()->SetTitle(title.c_str());
  ratio->GetXaxis()->SetTitle(title.c_str());
  data_cumu->GetXaxis()->SetTitle(title.c_str());

  data->GetXaxis()->SetRangeUser(xmin, xmax);
  ratio->GetXaxis()->SetRangeUser(xmin, xmax_ratio);
  data_cumu->GetXaxis()->SetRangeUser(xmin, xmax_cumu);

  //  const int fill_style_sig = 3001;
  const int fill_style_bgd = 3001;

  bgd->SetLineColor(kAzure+1);
  bgd->SetFillColor(kAzure+1);
  bgd->SetFillStyle(fill_style_bgd);

  bgd_cumu->SetLineColor(kAzure+1);
  bgd_cumu->SetFillColor(kAzure+1);
  bgd_cumu->SetFillStyle(3004);

  // =============== PLOT DISTRIBUTIONS HERE ========================
  TCanvas * c1 = new TCanvas();
  c1->SetLogy();

  if(data->GetMinimum() < 0.00001)data->SetMinimum(0.00001);
 
  data->Draw("e");
  bgd->Draw("same");

#if 0
  wp10->SetLineColor(kRed);
  wp15->SetLineColor(kRed+1);
  wp20->SetLineColor(kRed+2);
  wp10->SetFillColor(kRed);
  wp15->SetFillColor(kRed+1);
  wp20->SetFillColor(kRed+2);
  wp10->SetFillStyle(fill_style_sig);
  wp15->SetFillStyle(fill_style_sig);
  wp20->SetFillStyle(fill_style_sig);   
  wp10->Draw("same");
  wp15->Draw("same");
  wp20->Draw("same");
#endif
  

  TLegend * lg = new TLegend(0.59, 0.67, 0.89, 0.89);
  lg->SetTextSize(0.03);
  lg->SetBorderSize(0);
  lg->SetFillColor(0);
  lg->AddEntry(bgd, "Total bgd", "F");
#if 0
  lg->AddEntry(wp10, "W ' (1.0 TeV)", "F");
  lg->AddEntry(wp15, "W ' (1.5 TeV)", "F");
  lg->AddEntry(wp20, "W ' (2.0 TeV)", "F");
#endif
  lg->AddEntry(data, data_ipb, "LP");
  lg->Draw();
  
  c1->SaveAs(file.c_str());

  c1 = new TCanvas();
  c1->SetLogy();
		       
  //  data->Sumw2();
  // bgd->Sumw2();
  ratio->Divide(data, bgd);

  data_cumu->Draw("e");
  bgd_cumu->Draw("same");
  TLatex * l1c = new TLatex(350 + x_offset, 350, cms_prelim.c_str());
  l1c->SetTextSize(0.04); 
  TLatex * l2c = new TLatex(330 + x_offset, 100, lumi_sqrts.c_str());
  l2c->SetTextSize(0.04); 
  l1c->Draw(); l2c->Draw();

  TLegend * lgc = new TLegend(0.59, 0.67, 0.89, 0.89);
  lgc->SetTextSize(0.03);
  lgc->SetBorderSize(0);
  lgc->SetFillColor(0);
  lgc->AddEntry(bgd_cumu, "Total bgd", "F");
  lgc->AddEntry(data_cumu, data_ipb, "LP");
  lgc->Draw();

  c1->SaveAs(file_cumu.c_str());


  c1 = new TCanvas();

  ratio->SetMaximum(6);
  ratio->SetMinimum(-4);
  ratio->SetMarkerStyle(4);
  ratio->SetMarkerSize(1.3);
  ratio->Draw();

  TLatex * l1 = new TLatex(330 + x_offset, -1.9, cms_prelim.c_str());
  l1->SetTextSize(0.04); 
  TLatex * l2 = new TLatex(300 + x_offset, -3.0, lumi_sqrts.c_str());
  l2->SetTextSize(0.04); 
  l1->Draw(); l2->Draw();

  TLine * line = new TLine(xmin, 1, xmax_ratio, 1);
  line->Draw();

  c1->SaveAs(file_ratio.c_str());

  float lower_limit = 400;
  float upper_limit = 800;
  int bin_min = data->FindBin(lower_limit);
  int bin_max = data->FindBin(upper_limit);
  for(int bin = bin_min; bin <= bin_max; ++bin)
    {
      cout << " # of events above " << data->GetBinLowEdge(bin)
	   << " expected = " << bgd->Integral(bin, Nbins+1) 
	   << " observed = " << data->Integral(bin, Nbins+1) << endl;
    }
  

}
