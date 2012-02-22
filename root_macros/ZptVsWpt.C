#include <iostream>
#include <iomanip>

#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TError.h"
#include "TStyle.h"
#include "TLine.h"

#include "common.h"
#include "CMSStyle.C"


using namespace std;

void ZptVsWpt(){
  gErrorIgnoreLevel = kWarning;
  CMSstyle();
  gStyle->SetOptStat(0);

  TFile *fin = TFile::Open("../../../WprimeWZ.root", "read"); assert(fin);
//  TFile *fin = TFile::Open("../../../WprimeWZ.root", "read"); assert(fin);
//  TFile *fin = TFile::Open("PreApproval/Wprime_DoubleLepton.root", "read"); assert(fin);
  
  vector<string> vBkg, samples;
  vBkg.push_back("DYJetsToLL");
  vBkg.push_back("TTJets");
  vBkg.push_back("ZZ");
  vBkg.push_back("GVJets");  
  vBkg.push_back("WJetsToLNu");
  vBkg.push_back("WWTo2L2Nu");
  vBkg.push_back("WZJetsTo3LNu");

//  samples = vBkg;

  samples.push_back("Bkg");
  samples.push_back("Sig");
  samples.push_back("data");


  float minWindow = 0;
  float maxWindow = 20000;
  float minHt = 300;
  float minZpt = 0;
  float minWpt = 0;

  const string SignalName = "WprimeToWZTo3LNu_M-600";
  const string cuts = Form("(WZMass > %.0f && WZMass < %.0f && Ht > %.0f && Zpt > %.0f && Wpt > %.0f)*weight",
                           minWindow, maxWindow, minHt, minZpt, minWpt); 
  cout<<"Cuts are "<<cuts<<endl;
  //const string NoHtcuts = Form("(WZMass > %.0f && WZMass < %.0f && Ht > %.0f && Zpt > %.0f && Wpt > %.0f)",
  const string NoHtcuts = Form("(WZMass > %.0f && WZMass < %.0f && Ht > %.0f && Zpt > %.0f && Wpt > %.0f)*weight",
                               minWindow, maxWindow, 0., minZpt, minWpt);
//  const string NoHtcuts = Form("(Ht < 300 && Ht > WZMass/2 && WZMass > 500 && WZMass < 700)*weight"); 
//  const string NoHtcuts = Form("(WZMass > 500 && WZMass < 700)*weight"); 
  cout<<"No Ht Cuts are "<<NoHtcuts<<endl;

  vector<TH2D> hZptVsWpt  (samples.size(),TH2D("hZptVsWpt", ";p_{T}^{W} (GeV);p_{T}^{Z} (GeV)", 2500, 0, 500, 250, 0, 500));
  vector<TH2D> hWZMassVsHt(samples.size(),TH2D("hWZMassVsHt", ";H_{T} (GeV);M_{WZ} (GeV)", 5000, 0, 1000, 5000, 0, 1500));

  TCanvas c1,c2;
  c1.Divide(2,1);
  c2.Divide(2,1);
  TLegend legend(0.0,0.0,0.40,0.13,"");
	legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetNColumns(3);
  legend.SetTextSize(0.05);
  
  for(unsigned i=0; i<samples.size(); ++i){
    int color = kGray;
    string legendName;
    vector<string> names;
    if(samples[i] == "Bkg"){
      names = vBkg;
      legendName = "Bkg";
      color = kRed;
    }else if(samples[i] == "Sig"){
      names.push_back(SignalName);
      legendName = "Sig";
      color = kGreen;
    }else{
      names.push_back(samples[i]);
      legendName = samples[i];
      color = kBlack;
    }
//        cout<<"---------------------------------\n";
    cout<<"sample now is "<<samples[i]<<endl;

    get_sum_of_hists(fin, names, "tEvts_ValidWZCand", "Wpt:Zpt", cuts, hZptVsWpt[i]);
    get_sum_of_hists(fin, names, "tEvts_ValidWZCand", "Ht:WZMass", NoHtcuts, hWZMassVsHt[i]);
    
    cout<<" and # of X bins is "<<hWZMassVsHt[i].GetNbinsX()<<endl; 
    cout<<" and # of Y bins is "<<hWZMassVsHt[i].GetNbinsY()<<endl; 
    cout<<" and # of entries is "<<hWZMassVsHt[i].Integral()<<endl;
    cout<<" and max bin is "<<hWZMassVsHt[i].GetMaximumBin()<<endl;
    cout<<" and max # of entries is "<<hWZMassVsHt[i].GetBinContent(hWZMassVsHt[i].GetMaximumBin())<<endl;

    hZptVsWpt[i].SetMarkerColor(color);
    hZptVsWpt[i].SetFillColor(color);
    hZptVsWpt[i].SetMarkerSize(0.5);

    hWZMassVsHt[i].SetMarkerColor(color);
    hWZMassVsHt[i].SetFillColor(color);
    hWZMassVsHt[i].SetMarkerSize(0.5);
    
    if(samples[i] == "data"){
      hZptVsWpt[i].SetMarkerSize(0.4);
      hWZMassVsHt[i].SetMarkerSize(0.4);
    }

    legend.AddEntry(&hZptVsWpt[i],legendName.c_str(), "F");
//    legend.AddEntry(&hWZMassVsHt[i],legendName.c_str(), "F");

    cout<<"Drawing\n";

    string opt = samples[i] == "data" ? "scat" : "box";
    opt += i==0 ? "" : "same";
    
    c1.cd(1); hZptVsWpt[i].Draw(opt.c_str());
    if(samples[i] != "Sig"){
      c1.cd(2); hZptVsWpt[i].Draw(opt.c_str());
    }

    c2.cd(1); hWZMassVsHt[i].Draw(opt.c_str());
    if(samples[i] != "Sig"){
      c2.cd(2); hWZMassVsHt[i].Draw(opt.c_str());
    }             
  }
  c1.cd();
  legend.Draw();

  cout<<"Saving 1\n";
  c1.SaveAs("ZptVsWpt.png");
//  c1.SaveAs("ZptVsWpt_Wprime600.png");

  cout<<"Pad 2\n";
  c2.cd();
  legend.Draw();
/*
  cout<<"Adding Lines\n";
  TLine windLow (0, 500, 1000, 500); windLow.SetLineColor(kViolet);
  TLine windHigh(0, 700, 1000, 700); windHigh.SetLineColor(kViolet);

  TLine cutOld(minHt, 0, minHt, 1000); cutOld.SetLineColor(kBlue);
  TLine cutNew(0, 0, 500, 1000); cutNew.SetLineColor(kBlue);

  for(int i=0; i<2; ++i){
    c2.cd(i+1);

    windLow.Draw();
    windHigh.Draw();
    
    cutOld.Draw();
//  cutNew.Draw();
  }
*/
  float maxHt = 1000;
  TBox wprime200(  0,190,maxHt,210); 
  TBox wprime250(150,230,maxHt,270);
  TBox wprime300(160,280,maxHt,320);
  TBox wprime400(220,360,maxHt,440);
  TBox wprime500(230,450,maxHt,550);
  TBox wprime600(290,540,maxHt,660);
  TBox wprime700(360,620,maxHt,780);
  TBox wprime800(400,710,maxHt,890);
  
  wprime200.SetFillStyle(0);
  wprime250.SetFillStyle(0);
  wprime300.SetFillStyle(0);
  wprime400.SetFillStyle(0);
  wprime500.SetFillStyle(0);
  wprime600.SetFillStyle(0);
  wprime700.SetFillStyle(0);
  wprime800.SetFillStyle(0);
  /*
  wprime200.SetLineWidth(10);
  wprime250.SetLineWidth(10);
  wprime300.SetLineWidth(10);
  wprime400.SetLineWidth(10);
  wprime500.SetLineWidth(10);
  wprime600.SetLineWidth(10);
  wprime700.SetLineWidth(10);
  wprime800.SetLineWidth(10);
  */
  wprime200.SetLineColor(kBlue);
  wprime250.SetLineColor(kBlue);
  wprime300.SetLineColor(kBlue);
  wprime400.SetLineColor(kBlue);
  wprime500.SetLineColor(kBlue);
  wprime600.SetLineColor(kBlue);
  wprime700.SetLineColor(kBlue);
  wprime800.SetLineColor(kBlue);

  for(int i=0; i<2; ++i){
    c2.cd(i+1);
    wprime200.Draw("l");
    wprime250.Draw("l");
    wprime300.Draw("l");
    wprime400.Draw("l");
    wprime500.Draw("l");
    wprime600.Draw("l");
    wprime700.Draw("l");
    wprime800.Draw("l");
  }

  cout<<"Saving 2\n";
  c2.SaveAs("WZMassVsHt.png");
//  c2.SaveAs("WZMassVsHt_Wprime600.png");

}
