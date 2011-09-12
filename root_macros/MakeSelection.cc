#include <stdio.h>
#include <math.h>
#include <string>
#include <iostream>
#include "TFile.h"
#include "TTree.h" 
#include "TH1F.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TMath.h"
#include "TLegend.h"
#include "THStack.h"
#include <vector>
#include "TROOT.h"
#include "TError.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"

#include "common.h"
#include "CMSStyle.C"
using namespace std;

struct Cut{
  string name;
  bool isMin;
  float cutVal;

  Cut(){name=""; isMin = true; cutVal = 0;}
  Cut(string n,bool i=true,float c=0){name = n; isMin = i; cutVal = c;}
};

void MakeSelection(string inName, string opt="");
void DrawSelection(TFile* fin, Cut& thisCut);
bool Straddles(float a, float b, float num);

vector<string> SigSample;
vector<string> BkgSamples;
vector<Cut> Cuts;

float lumiUsed_;

bool debug_;

void
MakeSelection(string inName, string opt){
  gErrorIgnoreLevel = kWarning;
  CMSstyle();

  TFile *fin = TFile::Open(inName.c_str(), "read");
  lumiUsed_ = GetLumiUsed(fin);
  cout<<"Lumi Used is "<<lumiUsed_<<" inv pb\n";

  if(opt.find("debug") != string::npos) debug_ = true;

  vector<string> SigSamples;
  
  if(inName.find("Wprime") != string::npos || inName.find("EWKWZ") != string::npos){
    BkgSamples.push_back("WZJetsTo3LNu");
    BkgSamples.push_back("WJetsToLNu");
    BkgSamples.push_back("ZZ");
    BkgSamples.push_back("GVJets");
    BkgSamples.push_back("WWTo2L2Nu");
    BkgSamples.push_back("TTJets");
    BkgSamples.push_back("DYJetsToLL");
/*
  BkgSamples.push_back("ZBB0JetsToLNu");
  BkgSamples.push_back("ZBB1JetsToLNu");
  BkgSamples.push_back("ZBB2JetsToLNu");
  BkgSamples.push_back("ZBB3JetsToLNu");
  BkgSamples.push_back("ZCC0JetsToLNu");
  BkgSamples.push_back("ZCC1JetsToLNu");
  BkgSamples.push_back("ZCC2JetsToLNu");
  BkgSamples.push_back("ZCC3JetsToLNu");
*/
    SigSamples.push_back("WprimeToWZTo3LNu_M-300");
    SigSamples.push_back("WprimeToWZTo3LNu_M-400");
    //SigSamples.push_back("WprimeToWZTo3LNu_M-500");
    SigSamples.push_back("WprimeToWZTo3LNu_M-600");
    SigSamples.push_back("WprimeToWZTo3LNu_M-700");
    SigSamples.push_back("WprimeToWZTo3LNu_M-800");
    SigSamples.push_back("WprimeToWZTo3LNu_M-900");
    SigSamples.push_back("TC_WZ_300");
    SigSamples.push_back("TC_WZ_400");
    SigSamples.push_back("TC_WZ_500");
  
    Cuts.push_back(Cut("Ht", true));
    Cuts.push_back(Cut("Zpt", true));
    Cuts.push_back(Cut("Wpt", true));
  }else if(inName.find("HadVZ") != string::npos){
    BkgSamples.push_back("Summer11_ZZ");
    BkgSamples.push_back("Summer11_ZZJets_2l2q");
    BkgSamples.push_back("Summer11_VGamma");
    BkgSamples.push_back("Summer11_WW");
    BkgSamples.push_back("Summer11_WZ");
    BkgSamples.push_back("Summer11_TTJets");
    BkgSamples.push_back("Summer11_DYJetsToLL");

    SigSamples.push_back("Summer11_RSZZmmjj_750");
    SigSamples.push_back("Summer11_RSZZmmjj_1000");
    SigSamples.push_back("Summer11_RSZZmmjj_1250");
    SigSamples.push_back("Summer11_RSZZmmjj_1500");
    SigSamples.push_back("Summer11_RSZZmmjj_1750");
    SigSamples.push_back("Summer11_RSZZmmjj_2000");
    
    Cuts.push_back(Cut("Zpt", true));
    Cuts.push_back(Cut("Vpt", true));
  }

  for(uint i=0; i<SigSamples.size(); ++i){
    //Reset cuts for a new signal sample
    for(uint j=0; j<Cuts.size(); ++j){ 
      Cuts[j].cutVal = Cuts[j].isMin ? 0 : 9e9; 
    }

    SigSample.clear();
    SigSample.push_back(SigSamples[i]);
    TCanvas c1;
    c1.Print(Form("Selection_%s.pdf[", SigSample[0].c_str()), "pdf"); 
    
    for(uint j=0; j<Cuts.size(); ++j){
      DrawSelection(fin, Cuts[j]);
    }

    cout<<"For sample: "<<SigSample[0]<<" final cut is: ";
    for(uint j=0; j<Cuts.size(); ++j){
      cout<<" "<<Cuts[j].name<<" "<<Cuts[j].cutVal;
    }
    cout<<endl;

    c1.Print(Form("Selection_%s.pdf]", SigSample[0].c_str()), "pdf"); 
  }
}

void
DrawSelection(TFile* fin, Cut& thisCut){
  //printf("  Draw Selection\n");
  assert(SigSample.size() == 1);

  vector<TCut> cuts(Cuts.size());
  for(uint j=0; j<Cuts.size(); ++j){
    cuts[j] = Form("%s %s %.0f", Cuts[j].name.c_str(), Cuts[j].isMin ? ">" : "<", Cuts[j].cutVal);
  }
  
  TCut cutString("cutString");
  for(uint j=0; j<Cuts.size(); ++j){
    if(j==0) cutString = cuts[j];
    else     cutString = cutString && cuts[j];
  }
  cutString *= "weight";

  if(debug_) cout<<" cut is now "<<cutString<<endl;

  TH1F hSig("hSig","hSig",90,0,900);
  TH1F hBkg("hBkg","hBkg",90,0,900);
  get_sum_of_hists(fin, SigSample, "tWZCand", thisCut.name, cutString.GetTitle(), hSig);
  get_sum_of_hists(fin, BkgSamples, "tWZCand", thisCut.name, cutString.GetTitle(), hBkg);

  int first = 0;//Include under flow
  int last  = hSig.GetNbinsX() + 1;//Include over flow
  int nbins = hSig.GetNbinsX() + 2;

  const float Nsig_tot  = hSig.Integral(first, last);
  const float Nbkg_tot  = hBkg.Integral(first, last);

  //printf("Tot sig:%f   Tot Bkg:%f\n", Nsig_tot, Nbkg_tot);
  //printf("first:%i, last:%i, nbins:%i or %i\n", first, last, nbins, hists[0]->GetNbinsX());

  double* Nsig = new double[nbins];
  double* Nbkg = new double[nbins];
  double* fsig = new double[nbins];
  double* fbkg = new double[nbins];

  double* cut = new double[nbins];
  double* Signif = new double[nbins];

  for(int bin=first; bin<=last; ++bin){//for min cut
    //printf("bin:%i low:%f\n",bin, hists[0]->GetBinCenter(bin));
    Nsig[bin] = 0;
    Nbkg[bin] = 0;
    fsig[bin] = 0;
    fbkg[bin] = 0;
    cut[bin] = 0;
    Signif[bin] = 0;

    if(thisCut.isMin){
      Nsig[bin] = hSig.Integral(bin, last);
      Nbkg[bin] = hBkg.Integral(bin, last);
    }else{
      Nsig[bin] = hSig.Integral(first, bin);
      Nbkg[bin] = hBkg.Integral(first, bin);
    }
    fsig[bin] = Nsig[bin] / Nsig_tot;
    fbkg[bin] = Nbkg[bin] / Nbkg_tot;
    
    cut[bin] = hSig.GetBinCenter(bin);
    if(Nsig[bin] + Nbkg[bin] > 0) Signif[bin] = Nsig[bin] / sqrt(Nsig[bin] + Nbkg[bin]);
  }
  //Find the optimial cut
  int cutbin = 0;
  for(int bin=first; bin<last; ++bin){
    if(Signif[bin] > Signif[cutbin]) cutbin = bin;
  }


  TCanvas c1;
  hBkg.Draw();
  c1.Print(Form("Selection_%s.pdf",SigSample[0].c_str()), "pdf");
  c1.Clear();
  hSig.Draw();
  c1.Print(Form("Selection_%s.pdf",SigSample[0].c_str()), "pdf");
  c1.Clear();

  //c1.Divide(1,2);
  
  //TVirtualPad* p1 = c1.cd(1);
  string graph_name = "g" + thisCut.name;
  string graph_title = ";Signal Eff;Background Eff";

  TMultiGraph mg(graph_name.c_str(), graph_title.c_str());
  TGraph graph(nbins, fsig, fbkg);
  TGraph point(1, &fsig[cutbin], &fbkg[cutbin]);
  point.SetMarkerColor(kRed);

  mg.Add(&graph,"*");
  mg.Add(&point,"*");
    
  mg.Draw("A");

  TText tCut;
  tCut.SetNDC();
  tCut.SetTextSize(0.07);
  tCut.SetTextColor(2);
  tCut.DrawText(0.2, 0.75, Form("%s Cut=%.0f GeV", thisCut.name.c_str(), hSig.GetBinLowEdge(cutbin)));

  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextSize(0.045);
  latexLabel.DrawLatex(0.18, 0.96, "#font[132]{CMS Preliminary 2011}");
  latexLabel.DrawLatex(0.5, 0.955, "#font[132]{#sqrt{s} = 7 TeV}");
  latexLabel.DrawLatex(0.25, 0.87, Form("#font[132]{#intL dt = %.2f fb^{-1}}",lumiUsed_/1000.));

  c1.SaveAs(Form("plots/%s_%sEff.pdf", SigSample[0].c_str(), thisCut.name.c_str()));
  c1.Print(Form("Selection_%s.pdf",SigSample[0].c_str()), "pdf");
  c1.Clear();

  //////////////////////////////
  //TVirtualPad* p2 = c1.cd(2);
  graph_title = "; Cut Value (GeV);#frac{S}{#sqrt{S+B}}";
  TMultiGraph mgSignif(graph_name.c_str(), graph_title.c_str());
  TGraph gSignif(nbins, cut, Signif);
  TGraph gSignifPoint(1, &cut[cutbin], &Signif[cutbin]);
  gSignifPoint.SetMarkerColor(kRed);

  mgSignif.Add(&gSignif,"*");
  mgSignif.Add(&gSignifPoint,"*");
    
  mgSignif.Draw("A");
  

  TText tSignifCut;
  tSignifCut.SetNDC();
  tSignifCut.SetTextSize(0.07);
  tSignifCut.SetTextColor(2);
  //tSignif.SetTextAlign(22);
  tSignifCut.DrawText(0.60, 0.75, Form("%s Cut=%.0f GeV", thisCut.name.c_str(), hSig.GetBinLowEdge(cutbin)));

  latexLabel.DrawLatex(0.18, 0.96, "#font[132]{CMS Preliminary 2011}");
  latexLabel.DrawLatex(0.5, 0.955, "#font[132]{#sqrt{s} = 7 TeV}");
  latexLabel.DrawLatex(0.73, 0.87, Form("#font[132]{#intL dt = %.2f fb^{-1}}",lumiUsed_/1000.));

  c1.SaveAs(Form("plots/%s_%sSignificance.pdf", SigSample[0].c_str(), thisCut.name.c_str()));
  c1.Print(Form("Selection_%s.pdf",SigSample[0].c_str()), "pdf");


  thisCut.cutVal = hSig.GetBinLowEdge(cutbin);
  //if (thisCut.name == "Ht")  minHt_  = hSig.GetBinLowEdge(cutbin);
  //if (thisCut.name == "Zpt") minZpt_ = hSig.GetBinLowEdge(cutbin);
  //if (thisCut.name == "Wpt") minWpt_ = hSig.GetBinLowEdge(cutbin);

}

bool 
Straddles(float a, float b, float num){
  if(a >= num && num >= b) return true;
  if(b >= num && num >= a) return true;
    
  return false;
}
