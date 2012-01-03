/* 
   Usage: root -b -q 'MakeSelection.cc+("input.root", "options")'
   Options include:
   "debug" = show debuging statements
   "hist" = use histograms instead of default trees to make cuts
   "eff" = choose cut which gives 98% signal eff instead of peak significance
   
 */

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

bool debug_ = false;
bool useHists_ = false;
bool useSig_ = true;
string treeName = "tWZCand";

void
MakeSelection(string inName, string opt){
  gErrorIgnoreLevel = kWarning;
  CMSstyle();

  if(opt.find("debug") != string::npos) debug_ = true;
  if(opt.find("hist") != string::npos) useHists_ = true;
  if(opt.find("eff") != string::npos) useSig_ = false;

  TFile *fin = TFile::Open(inName.c_str(), "read"); assert(fin);

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
    SigSamples.push_back("WprimeToWZTo3LNu_M-200");
    SigSamples.push_back("WprimeToWZTo3LNu_M-250");
    SigSamples.push_back("WprimeToWZTo3LNu_M-300");
    SigSamples.push_back("WprimeToWZTo3LNu_M-400");
    SigSamples.push_back("WprimeToWZTo3LNu_M-500");
    SigSamples.push_back("WprimeToWZTo3LNu_M-600");
    SigSamples.push_back("WprimeToWZTo3LNu_M-700");
    SigSamples.push_back("WprimeToWZTo3LNu_M-800");
    SigSamples.push_back("WprimeToWZTo3LNu_M-900");
    SigSamples.push_back("WprimeToWZTo3LNu_M-1000");
    SigSamples.push_back("WprimeToWZTo3LNu_M-1100");
    SigSamples.push_back("WprimeToWZTo3LNu_M-1200");
    SigSamples.push_back("WprimeToWZTo3LNu_M-1300");
    SigSamples.push_back("WprimeToWZTo3LNu_M-1400");
    SigSamples.push_back("WprimeToWZTo3LNu_M-1500");
    SigSamples.push_back("TC_WZ_300");
    SigSamples.push_back("TC_WZ_400");
    SigSamples.push_back("TC_WZ_500");
    SigSamples.push_back("TC_WZ_600");
    SigSamples.push_back("TC_WZ_700");
    SigSamples.push_back("TC_WZ_800");
    SigSamples.push_back("TC_WZ_900");
  
    Cuts.push_back(Cut("Ht", true));
    //Cuts.push_back(Cut("Zpt", true));
    //Cuts.push_back(Cut("Wpt", true));

    treeName = "tEvts_ValidWZCand";

  }else if(inName.find("HadVZ") != string::npos){
    BkgSamples.push_back("Summer11_ZZ");
    BkgSamples.push_back("Summer11_VGamma");
    BkgSamples.push_back("Summer11_WW");
    BkgSamples.push_back("Summer11_WZ");
    BkgSamples.push_back("Summer11_TTJets");
    BkgSamples.push_back("Summer11_DYJetsToLL_PtZ100");

    SigSamples.push_back("Summer11_RSZZmmjj_750");
    SigSamples.push_back("Summer11_RSZZmmjj_1000");
    SigSamples.push_back("Summer11_RSZZmmjj_1250");
    SigSamples.push_back("Summer11_RSZZmmjj_1500");
    SigSamples.push_back("Summer11_RSZZmmjj_1750");
    SigSamples.push_back("Summer11_RSZZmmjj_2000");

    SigSamples.push_back("Summer11_WPrimeZZlljj_300");
    SigSamples.push_back("Summer11_WPrimeZZlljj_400");
    //SigSamples.push_back("Summer11_WPrimeZZlljj_500");
    SigSamples.push_back("Summer11_WPrimeZZlljj_600");
    SigSamples.push_back("Summer11_WPrimeZZlljj_700");
    SigSamples.push_back("Summer11_WPrimeZZlljj_800");
    SigSamples.push_back("Summer11_WPrimeZZlljj_900");
    SigSamples.push_back("Summer11_WPrimeZZlljj_1000");
    //SigSamples.push_back("Summer11_WPrimeZZlljj_1100");
    SigSamples.push_back("Summer11_WPrimeZZlljj_1200");
    //SigSamples.push_back("Summer11_WPrimeZZlljj_1300");
    SigSamples.push_back("Summer11_WPrimeZZlljj_1400");
    SigSamples.push_back("Summer11_WPrimeZZlljj_1500");
    
    Cuts.push_back(Cut("Zpt", true));
    Cuts.push_back(Cut("Vpt", true));

    treeName = "tVZCand";
  }

  if(debug_) cout<<"Using "<<SigSamples.size()<<" signal samples\n"
                 <<"Using "<<BkgSamples.size()<<" background samples\n"
                 <<"Using "<<Cuts.size()<<" cuts\n"<<endl;

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
    cout<<endl<<endl;

    c1.Print(Form("Selection_%s.pdf]", SigSample[0].c_str()), "pdf"); 
  }
}

void
DrawSelection(TFile* fin, Cut& thisCut){
  if(debug_) printf("  Draw Selection\n");
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
  if(!useHists_){
    get_sum_of_hists(fin, SigSample, treeName, thisCut.name, cutString.GetTitle(), hSig);
    get_sum_of_hists(fin, BkgSamples, treeName, thisCut.name, cutString.GetTitle(), hBkg);
  }else{
    string histName;
    if(thisCut.name == "Zpt") histName = "hZpt_ZMass";
    if(thisCut.name == "Vpt") histName = "hVpt_VMass";
    if(debug_) cout<<"Using histName: "<<histName<<endl;
    TH1F* pSig = get_sum_of_hists(fin, SigSample, histName);
    TH1F* pBkg = get_sum_of_hists(fin, BkgSamples, histName);
    hSig = *pSig;
    hBkg = *pBkg;
  }
  
  int first = 0;//Include under flow
  int last  = hSig.GetNbinsX() + 1;//Include over flow
  int nbins = hSig.GetNbinsX() + 2;

  const float Nsig_tot  = hSig.Integral(first, last);
  const float Nbkg_tot  = hBkg.Integral(first, last);

  if(debug_) printf("Tot sig:%f   Tot Bkg:%f\n", Nsig_tot, Nbkg_tot);
  //if(debug_) printf("first:%i, last:%i, nbins:%i or %i\n", first, last, nbins, hists[0]->GetNbinsX());

  double* Nsig = new double[nbins];
  double* Nbkg = new double[nbins];
  double* fsig = new double[nbins];
  double* fbkg = new double[nbins];

  double* cut = new double[nbins];
  double* Signif = new double[nbins];

  for(int bin=first; bin<=last; ++bin){
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
  int cutBinSig = 0;
  int cutBinEff = 0;
  for(int bin=first; bin<last; ++bin){
    if(Signif[bin] > Signif[cutBinSig]) cutBinSig = bin;
    if(Straddles(fsig[bin],fsig[bin+1], 0.98)) 
      cutBinEff = fsig[bin+1] > fsig[bin] ? bin+1 : bin;
  }

  if(useSig_) thisCut.cutVal = hSig.GetBinLowEdge(cutBinSig);
  else        thisCut.cutVal = hSig.GetBinLowEdge(cutBinEff);

  if(debug_){
    cout<<" For the significance cut, the sig eff is "
        <<fsig[cutBinSig]*100<<"\% and bkg eff is "<<fbkg[cutBinSig]*100
        <<"\%"<<endl;
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
  TGraph pointSig(1, &fsig[cutBinSig], &fbkg[cutBinSig]);
  pointSig.SetMarkerColor(kRed);
  TGraph pointEff(1, &fsig[cutBinEff], &fbkg[cutBinEff]);
  pointEff.SetMarkerColor(kGreen);

  mg.Add(&graph,"*");
  mg.Add(&pointSig,"*");
  mg.Add(&pointEff,"*");
    
  mg.Draw("A");

  string cutResult = thisCut.name + (thisCut.isMin ? " >" : " <");

  TText tCutSig;
  tCutSig.SetNDC();
  tCutSig.SetTextSize(0.07);
  tCutSig.SetTextColor(kRed);
  tCutSig.DrawText(0.2, 0.75, Form("%s %.0f GeV (Sig)", cutResult.c_str(), hSig.GetBinLowEdge(cutBinSig)));

  TText tCutEff;
  tCutEff.SetNDC();
  tCutEff.SetTextSize(0.07);
  tCutEff.SetTextColor(kGreen);
  tCutEff.DrawText(0.2, 0.65, Form("%s %.0f GeV (Eff)", cutResult.c_str(), hSig.GetBinLowEdge(cutBinEff)));

  c1.SaveAs(Form("plots/%s_%sEff.pdf", SigSample[0].c_str(), thisCut.name.c_str()));
  c1.Print(Form("Selection_%s.pdf",SigSample[0].c_str()), "pdf");
  c1.Clear();

  //////////////////////////////
  //TVirtualPad* p2 = c1.cd(2);
  graph_title = "; Cut Value (GeV);#frac{N_{S}}{#sqrt{N_{S}+N_{B}}}";
  TMultiGraph mgSignif(graph_name.c_str(), graph_title.c_str());
  TGraph gSignif(nbins, cut, Signif);
  TGraph gSignifPointSig(1, &cut[cutBinSig], &Signif[cutBinSig]);
  TGraph gSignifPointEff(1, &cut[cutBinEff], &Signif[cutBinEff]);
  gSignifPointSig.SetMarkerColor(kRed);
  gSignifPointEff.SetMarkerColor(kGreen);

  mgSignif.Add(&gSignif,"*");
  mgSignif.Add(&gSignifPointSig,"*");
  mgSignif.Add(&gSignifPointEff,"*");
    
  mgSignif.Draw("A");
  
  tCutSig.DrawText(0.2, 0.75, Form("%s %.0f GeV (Sig)", cutResult.c_str(), hSig.GetBinLowEdge(cutBinSig)));
  tCutEff.DrawText(0.2, 0.65, Form("%s %.0f GeV (Eff)", cutResult.c_str(), hSig.GetBinLowEdge(cutBinEff)));

  c1.SaveAs(Form("plots/%s_%sSignificance.pdf", SigSample[0].c_str(), thisCut.name.c_str()));
  c1.Print(Form("Selection_%s.pdf",SigSample[0].c_str()), "pdf");


}

bool 
Straddles(float a, float b, float num){
  if(a >= num && num >= b) return true;
  if(b >= num && num >= a) return true;
    
  return false;
}
