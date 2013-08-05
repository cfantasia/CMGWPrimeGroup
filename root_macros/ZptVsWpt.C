//Author: Cory Fantasia 2012
//Purpose: Plot 2D plots e.g. Zpt vs Wpt
//Usage: root -b -l -q 'ZptVsWpt.C+(string inName, string outSuffix="")'
//Usage: root -b -l -q 'ZptVsWpt.C+("../../../WprimeWZ.root", "")'

#include <iostream>
#include <iomanip>

#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TError.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"

#include "common.h"
#include "CMSStyle.C"
#include "setTDR_modified.C"


using namespace std;

string outSuffix_;

void
Make2DPlot(const string & name, const string & title, 
           const int & nx, const float & minx, const float & maxx, 
           const int & ny, const float & miny, const float & maxy, 
           TFile* fin, const string & plot, const string & cuts);

void DrawBoxes(TCanvas & c);


void ZptVsWpt(string inName="../../../WprimeWZ.root", string outSuffix=""){
  gErrorIgnoreLevel = kWarning;
  setTDRStyle();
  //CMSstyle();
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  outSuffix_ = outSuffix;

  TFile *fin = TFile::Open(inName.c_str(), "read"); assert(fin);
  
  TFile *fout = TFile::Open("Wprime-2DPlots.root", "recreate"); assert(fout);

  float minWindow = 0;
  float maxWindow = 20000;
  float minLt = 400;
  float minZpt = 0;
  float minWpt = 0;

  const string cuts = Form("(WZMass > %.0f && WZMass < %.0f && Lt > %.0f && Zpt > %.0f && Wpt > %.0f)*weight",
                           minWindow, maxWindow, minLt, minZpt, minWpt); 
  cout<<"Cuts are "<<cuts<<endl;
  //const string NoLtcuts = Form("(WZMass > %.0f && WZMass < %.0f && Lt > %.0f && Zpt > %.0f && Wpt > %.0f)",
  const string NoLtcuts = Form("(WZMass > %.0f && WZMass < %.0f && Lt > %.0f && Zpt > %.0f && Wpt > %.0f)*weight",
                               minWindow, maxWindow, 0., minZpt, minWpt);
//  const string NoLtcuts = Form("(Lt < 300 && Lt > WZMass/2 && WZMass > 500 && WZMass < 700)*weight"); 
//  const string NoLtcuts = Form("(WZMass > 500 && WZMass < 700)*weight"); 
  cout<<"No Lt Cuts are "<<NoLtcuts<<endl;

  
  //Make2DPlot(name, title, nx, minx, maxx, ny, miny, maxy, fin, plots, cuts, outname);
  //Make2DPlot("ZptVsWpt", ";p_{T}^{W} (GeV);p_{T}^{Z} (GeV)", 150, 0, 750, 150, 0, 750, fin, "Wpt:Zpt", cuts);

  //Make2DPlot("ZptVsWpt", ";p_{T}^{W} (GeV);p_{T}^{Z} (GeV)", 150, 0, 750, 150, 0, 750, fin, "Wpt:Zpt", NoLtcuts);
  Make2DPlot("WZMassVsLt", ";L_{T} (GeV);M_{WZ} (GeV)", 1000, 0, 1000, 1500, 0, 1500, fin, "Lt:WZMass", NoLtcuts);
  //Make2DPlot("WZMassVsMeff", ";M_{eff} (GeV);M_{WZ} (GeV)", 5000, 0, 1500, 5000, 0, 1500, fin, "Lt+MET:WZMass", NoLtcuts);
  //Make2DPlot("WZMassVsZpt", ";p_{T}^{Z} (GeV);M_{WZ} (GeV)", 150, 0, 750, 5000, 0, 1500, fin, "Zpt:WZMass", NoLtcuts);
  //Make2DPlot("WZMassVsWpt", ";p_{T}^{W} (GeV);M_{WZ} (GeV)", 150, 0, 750, 5000, 0, 1500, fin, "Wpt:WZMass", NoLtcuts);
  //Make2DPlot("WZMassVsZDr", ";#Delta_{R}^{ll} (GeV);M_{WZ} (GeV)", 500, 0, 5, 5000, 0, 1500, fin, "ZDr:WZMass", NoLtcuts);
  //Make2DPlot("MeffVsLt", ";L_{T} (GeV);M_{eff} (GeV)", 5000, 0, 1500, 5000, 0, 1500, fin, "Lt:Lt+MET", NoLtcuts);
  //Make2DPlot("ZdrVsLt", ";L_{T} (GeV);#Delta_{R}^{ll} (GeV)", 5000, 0, 1000, 500, 0, 5, fin, "Lt:ZDr", NoLtcuts);
  //Make2DPlot("ZptVsZDr", ";#Delta_{R}^{ll} (GeV);p_{T}^{Z} (GeV)", 500, 0, 5, 150, 0, 750, fin, "ZDr:Zpt", NoLtcuts);
  //Make2DPlot("ZmassVsWTMass", ";M_{Z} (GeV);M_{W}^{T} (GeV)", 60, 60, 120, 50, 0, 100, fin, "ZMass:WTransMass", NoLtcuts);



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

  fout->Close();

}

void
Make2DPlot(const string & name, const string & title, 
           const int & nx, const float & minx, const float & maxx, 
           const int & ny, const float & miny, const float & maxy, 
           TFile* fin, const string & plot, const string & cuts){

  int mass = 1000;//signal sample mass
  const string SignalName = Form("WprimeToWZTo3LNu_M-%i",mass);

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

  vector<TH2D> h2D(samples.size(),TH2D(name.c_str(), title.c_str(), nx, minx, maxx, ny, miny, maxy));

  TCanvas c(name.c_str());
  //c.Divide(2,1);
  //TLegend legend(0.0,0.02,0.40,0.11,"");
  TLegend legend(0.17,0.7,0.37,0.93,"");
  prepLegend(&legend);
  //legend.SetNColumns(3);

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
    
    cout<<"sample now is "<<samples[i]<<endl;
    h2D[i].SetStats(0);
    get_sum_of_hists(fin, names, "tEvts_MET", plot.c_str(), cuts, h2D[i]);
    
    h2D[i].SetMarkerColor(color);
    h2D[i].SetFillColor(color);
    h2D[i].SetMarkerSize(0.5);
    
    if(samples[i] == "Sig") h2D[i].SetMarkerSize(0.6);
    if(samples[i] == "data") h2D[i].SetMarkerSize(0.4);
    
    legend.AddEntry(&h2D[i],legendName.c_str(), samples[i] == "data" ? "P" : "F");
    
    string opt = samples[i] == "data" ? "scat" : "box";
    opt += i==0 ? "" : "same";
    
    //c.cd(1); 
    h2D[i].Draw(opt.c_str());
    if(0 && samples[i] != "Sig"){
      c.cd(2); 
      h2D[i].Draw(opt.c_str());
    }
    h2D[i].GetXaxis()->SetNdivisions(505);
    
    //c.cd();
    if(i == (samples.size()-1)){
      if(name.find("WZMassVsLt") != string::npos) DrawBoxes(c);

      legend.Draw();

      TLatex latexLabel;
      latexLabel.SetNDC();
      latexLabel.SetTextSize(0.05);
      latexLabel.SetTextFont(42);
      latexLabel.DrawLatex(0.33, 0.96, "CMS Preliminary 2012");
      latexLabel.DrawLatex(0.37, 0.84, "#sqrt{s} = 8 TeV");
      latexLabel.DrawLatex(0.60, 0.84, Form("#intL dt = 19.6 fb^{-1}"));
    }

  }//loop over samples
  c.SaveAs(Form("%s%s.png", name.c_str(), outSuffix_.c_str()));
  c.Write();
}

void
DrawBox(TCanvas & c, int mass){
  float maxLt = 1000;//upper end of the box
  TBox* b = new TBox(LtCut(mass),mass - WindowWidth(mass)/2.,maxLt,mass + WindowWidth(mass)/2.); 
  
  b->SetFillStyle(0);
  b->SetLineColor(kBlue);
  for(int i=0; i<2; ++i){
    c.cd(i+1);
    b->Draw("l");
  }
}
void
DrawBoxes(TCanvas & c){
  float maxLt = 1000;//upper end of the box
  for(int mass=200; mass<=1000; mass+=100){
    TBox* b = new TBox(LtCut(mass),mass - WindowWidth(mass)/2.,maxLt,mass + WindowWidth(mass)/2.); 
    
    b->SetFillStyle(0);
    b->SetLineColor(kBlue);
    for(int i=0; i<2; ++i){
      c.cd(i+1);
      b->Draw("l");
    }
  }
  DrawBox(c, 170);
  DrawBox(c, 180);
  DrawBox(c, 190);
  DrawBox(c, 210);
  DrawBox(c, 220);
  DrawBox(c, 230);
  DrawBox(c, 240);
  DrawBox(c, 250);
  DrawBox(c, 275);
  DrawBox(c, 325);
  DrawBox(c, 350);
  DrawBox(c, 450);
  DrawBox(c, 550);
  
}
