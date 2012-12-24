//Usage: root -b -l -q ZptVsWpt.C+

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
  
  vector<string> vBkg, samples;
  //vBkg.push_back("DYJetsToLL");
  //vBkg.push_back("TTJets");
  //vBkg.push_back("ZZ");
  //vBkg.push_back("GVJets");  
  //vBkg.push_back("WJetsToLNu");
  //vBkg.push_back("WWTo2L2Nu");
  vBkg.push_back("WZJetsTo3LNu");

//  samples = vBkg;

  samples.push_back("Bkg");
  samples.push_back("Sig");
  samples.push_back("data");


  float minWindow = 0;
  float maxWindow = 20000;
  float minLt = 300;
  float minZpt = 0;
  float minWpt = 0;

  const string SignalName = "WprimeToWZTo3LNu_M-1000-MyGen";
  const string cuts = Form("(WZMass > %.0f && WZMass < %.0f && Lt > %.0f && Zpt > %.0f && Wpt > %.0f)*weight",
                           minWindow, maxWindow, minLt, minZpt, minWpt); 
  cout<<"Cuts are "<<cuts<<endl;
  //const string NoLtcuts = Form("(WZMass > %.0f && WZMass < %.0f && Lt > %.0f && Zpt > %.0f && Wpt > %.0f)",
  const string NoLtcuts = Form("(WZMass > %.0f && WZMass < %.0f && Lt > %.0f && Zpt > %.0f && Wpt > %.0f)*weight",
                               minWindow, maxWindow, 0., minZpt, minWpt);
//  const string NoLtcuts = Form("(Lt < 300 && Lt > WZMass/2 && WZMass > 500 && WZMass < 700)*weight"); 
//  const string NoLtcuts = Form("(WZMass > 500 && WZMass < 700)*weight"); 
  cout<<"No Lt Cuts are "<<NoLtcuts<<endl;

  vector<TH2D> hZptVsWpt  (samples.size(),TH2D("hZptVsWpt", ";p_{T}^{W} (GeV);p_{T}^{Z} (GeV)", 250, 0, 750, 250, 0, 750));
  vector<TH2D> hWZMassVsLt(samples.size(),TH2D("hWZMassVsLt", ";L_{T} (GeV);M_{WZ} (GeV)", 5000, 0, 1000, 5000, 0, 1500));
  vector<TH2D> hWZMassVsMeff(samples.size(),TH2D("hWZMassVsMeff", ";M_{eff} (GeV);M_{WZ} (GeV)", 5000, 0, 1500, 5000, 0, 1500));
  vector<TH2D> hWZMassVsZDr(samples.size(),TH2D("hWZMassVsZDr", ";#Delta_{R}^{ll} (GeV);M_{WZ} (GeV)", 500, 0, 5, 5000, 0, 1500));
  vector<TH2D> hZmassVsWTMass(samples.size(),TH2D("hZmassVsWTMass", ";M_{Z} (GeV);M_{W}^{T} (GeV)", 60, 60, 120, 100, 0, 100));

  TCanvas c[5]; for(int i=0; i<5;  ++i) c[i].Divide(2,1);
  TLegend legend(0.0,0.02,0.40,0.11,"");
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
    }else if(samples[i] == "data"){
      names.push_back(samples[i]);
      legendName = "Data";
      color = kBlack;
    }else{
      names.push_back(samples[i]);
      legendName = samples[i];
      color = kYellow;
    }

//        cout<<"---------------------------------\n";
    cout<<"sample now is "<<samples[i]<<endl;

    get_sum_of_hists(fin, names, "tEvts_MET", "Wpt:Zpt", cuts, hZptVsWpt[i]);
    get_sum_of_hists(fin, names, "tEvts_MET", "Lt:WZMass", NoLtcuts, hWZMassVsLt[i]);
    get_sum_of_hists(fin, names, "tEvts_MET", "Lt+MET:WZMass", NoLtcuts, hWZMassVsMeff[i]);
    get_sum_of_hists(fin, names, "tEvts_MET", "ZDr:WZMass", NoLtcuts, hWZMassVsZDr[i]);
    get_sum_of_hists(fin, names, "tEvts_MET", "ZMass:WTransMass", NoLtcuts, hZmassVsWTMass[i]);
    
    cout<<" and # of X bins is "<<hWZMassVsLt[i].GetNbinsX()<<endl; 
    cout<<" and # of Y bins is "<<hWZMassVsLt[i].GetNbinsY()<<endl; 
    cout<<" and # of entries is "<<hWZMassVsLt[i].Integral()<<endl;
    cout<<" and max bin is "<<hWZMassVsLt[i].GetMaximumBin()<<endl;
    cout<<" and max # of entries is "<<hWZMassVsLt[i].GetBinContent(hWZMassVsLt[i].GetMaximumBin())<<endl;

    hZptVsWpt[i].SetMarkerColor(color);
    hZptVsWpt[i].SetFillColor(color);
    hZptVsWpt[i].SetMarkerSize(0.5);

    hWZMassVsLt[i].SetMarkerColor(color);
    hWZMassVsLt[i].SetFillColor(color);
    hWZMassVsLt[i].SetMarkerSize(0.5);

    hWZMassVsMeff[i].SetMarkerColor(color);
    hWZMassVsMeff[i].SetFillColor(color);
    hWZMassVsMeff[i].SetMarkerSize(0.5);

    hWZMassVsZDr[i].SetMarkerColor(color);
    hWZMassVsZDr[i].SetFillColor(color);
    hWZMassVsZDr[i].SetMarkerSize(0.5);

    hZmassVsWTMass[i].SetMarkerColor(color);
    hZmassVsWTMass[i].SetFillColor(color);
    hZmassVsWTMass[i].SetMarkerSize(0.5);
    
    if(samples[i] == "Sig"){
      hZptVsWpt[i].SetMarkerSize(0.6);
      hWZMassVsLt[i].SetMarkerSize(0.6);
      hWZMassVsMeff[i].SetMarkerSize(0.6);
      hWZMassVsZDr[i].SetMarkerSize(0.6);
      hZmassVsWTMass[i].SetMarkerSize(0.6);
    }

    if(samples[i] == "data"){
      hZptVsWpt[i].SetMarkerSize(0.4);
      hWZMassVsLt[i].SetMarkerSize(0.4);
      hWZMassVsMeff[i].SetMarkerSize(0.4);
      hWZMassVsZDr[i].SetMarkerSize(0.4);
      hZmassVsWTMass[i].SetMarkerSize(0.4);
    }

    legend.AddEntry(&hZptVsWpt[i],legendName.c_str(), "F");
//    legend.AddEntry(&hWZMassVsLt[i],legendName.c_str(), "F");

    cout<<"Drawing\n";

    string opt = samples[i] == "data" ? "scat" : "box";
    opt += i==0 ? "" : "same";
    
    c[0].cd(1); hZptVsWpt[i].Draw(opt.c_str());
    if(samples[i] != "Sig"){
      c[0].cd(2); hZptVsWpt[i].Draw(opt.c_str());
    }

    c[1].cd(1); hWZMassVsLt[i].Draw(opt.c_str());
    if(samples[i] != "Sig"){
      c[1].cd(2); hWZMassVsLt[i].Draw(opt.c_str());
    }             

    c[2].cd(1); hWZMassVsMeff[i].Draw(opt.c_str());
    if(samples[i] != "Sig"){
      c[2].cd(2); hWZMassVsMeff[i].Draw(opt.c_str());
    }             

    c[3].cd(1); hWZMassVsZDr[i].Draw(opt.c_str());
    if(samples[i] != "Sig"){
      c[3].cd(2); hWZMassVsZDr[i].Draw(opt.c_str());
    }             

    c[4].cd(1); hZmassVsWTMass[i].Draw(opt.c_str());
    if(samples[i] != "Sig"){
      c[4].cd(2); hZmassVsWTMass[i].Draw(opt.c_str());
    }             
  }
  c[0].cd();
  legend.Draw();

  cout<<"Saving 1\n";
  c[0].SaveAs("ZptVsWpt.png");
//  c1.SaveAs("ZptVsWpt_Wprime600.png");

  cout<<"Pad 2\n";
  c[1].cd();
  legend.Draw();

  cout<<"Pad 3\n";
  c[2].cd();
  legend.Draw();

  cout<<"Pad 4\n";
  c[3].cd();
  legend.Draw();

  cout<<"Pad 5\n";
  c[4].cd();
  legend.Draw();
/*
  cout<<"Adding Lines\n";
  TLine windLow (0, 500, 1000, 500); windLow.SetLineColor(kViolet);
  TLine windHigh(0, 700, 1000, 700); windHigh.SetLineColor(kViolet);

  TLine cutOld(minLt, 0, minLt, 1000); cutOld.SetLineColor(kBlue);
  TLine cutNew(0, 0, 500, 1000); cutNew.SetLineColor(kBlue);

  for(int i=0; i<2; ++i){
    c2.cd(i+1);

    windLow.Draw();
    windHigh.Draw();
    
    cutOld.Draw();
//  cutNew.Draw();
  }
*/
  float maxLt = 1000;
  TBox wprime200(  0,190,maxLt,210); 
  TBox wprime250(150,230,maxLt,270);
  TBox wprime300(160,280,maxLt,320);
  TBox wprime400(220,360,maxLt,440);
  TBox wprime500(230,450,maxLt,550);
  TBox wprime600(290,540,maxLt,660);
  TBox wprime700(360,620,maxLt,780);
  TBox wprime800(400,710,maxLt,890);
  
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
    c[1].cd(i+1);
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
  c[1].SaveAs("WZMassVsLt.png");
//  c2.SaveAs("WZMassVsLt_Wprime600.png");

  cout<<"Saving 3\n";
  c[2].SaveAs("WZMassVsMeff.png");

  cout<<"Saving 3\n";
  c[3].SaveAs("WZMassVsZDr.png");

  cout<<"Saving 4\n";
  c[4].SaveAs("ZmassVsWTMass.png");
}

//  LocalWords:  str
