#include <stdio.h>
#include <math.h>
#include <string>
#include <iostream>
#include "TFile.h"
#include "TTree.h" 
#include "TH1F.h"
#include "TMultiGraph.h"
#include "TGraph.h"
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

struct Sample{
  std::string name;
  std::vector<std::string> names;
  int line;
  int style;
  int fill;
  TH1F* hist;
  Sample(){
    line = 0; style = 0; fill=0; hist=NULL;
  }
  Sample(std::string n){
    name = n; names.push_back(n); line = 1; style = 0; fill = 0; hist=NULL;
  }
  Sample(std::string n, int l, int s, int f){
    name = n; names.push_back(n); line = l; style = s; fill = f; hist=NULL;
  }
  Sample(std::string n, std::vector<std::string> ns, int l, int s, int f){
    name = n; names = ns; line = l; style = s; fill = f; hist=NULL;
  }
};

std::vector< std::vector<Sample>* > samples_;
std::vector<Sample> Data;
std::vector<Sample> Bkg;
std::vector<Sample> Sig;
std::map<std::string, std::string> SampleNames;

const bool debug = false;
float lumiUsed_ = 0.;

int main(int argc, char ** argv);
void MakePlots(string inName, string outName, string opt="");
void DrawandSave(TFile* fin, std::string pdfName, std::string title, std::string bookmark, bool logy, bool eff, bool cum, TLine* line=NULL);
TH1F* FillCum(TH1F* h);
void GetHistograms(TFile* fin, std::string title, bool eff, bool cum);
void Draw(std::string filename, std::string pdfName, std::string bookmark, bool logy, bool eff, TLine* line);
void CheckSamples(TFile* fin, std::vector<Sample> & sample);
TH1F* GetValidHist(TFile* f);
vector<string> GetListofCuts(TFile* f, TH1F* hist=NULL);

void
MakePlots(string inName, string outName, string opt){  
  gErrorIgnoreLevel = kWarning;
  CMSstyle();
  
  TFile *fin = TFile::Open(inName.c_str(), "read");//Cory:
  lumiUsed_ = GetLumiUsed(fin);
  cout<<"Cory: Using lumiUsed Hack bc merging histos!!!\n";

  cout<<"Lumi Used is "<<lumiUsed_<<" inv pb\n";

  SampleNames["data"]="Data";
  SampleNames["WJetsToLNu"]="W+Jets";
  SampleNames["WWTo2L2Nu"]="WW\\rightarrow2l2\\nu";
  SampleNames["TTJets"]="t\\bar{t}";
  SampleNames["ZZTo4L_TuneZ2"]="ZZ\\rightarrow4l";
  SampleNames["PhotonVJets"]="V\\gamma";
  SampleNames["VV"]="VV";
  SampleNames["DYJetsToLL"]="DY+Jets\\rightarrow2l";
  SampleNames["ZBB0JetsToLNu"]="Z+0Jets\\rightarrowbb";
  SampleNames["ZBB1JetsToLNu"]="Z+1Jets\\rightarrowbb";
  SampleNames["ZBB2JetsToLNu"]="Z+2Jets\\rightarrowbb";
  SampleNames["ZBB3JetsToLNu"]="Z+3Jets\\rightarrowbb";
  SampleNames["ZCC0JetsToLNu"]="Z+0Jets\\rightarrowcc";
  SampleNames["ZCC1JetsToLNu"]="Z+1Jets\\rightarrowcc";
  SampleNames["ZCC2JetsToLNu"]="Z+2Jets\\rightarrowcc";
  SampleNames["ZCC3JetsToLNu"]="Z+3Jets\\rightarrowcc";
  SampleNames["WZTo3LNu"]="WZ\\rightarrow3l\\nu";
  SampleNames["WprimeToWZTo3LNu_M-300"]="W' 300";
  SampleNames["WprimeToWZTo3LNu_M-400"]="W' 400";
  SampleNames["WprimeToWZTo3LNu_M-500"]="W' 500";
  SampleNames["WprimeToWZTo3LNu_M-600"]="W' 600";
  SampleNames["WprimeToWZTo3LNu_M-700"]="W' 700";
  SampleNames["WprimeToWZTo3LNu_M-800"]="W' 800";
  SampleNames["WprimeToWZTo3LNu_M-900"]="W' 900";
  SampleNames["TC225"]="\\rho_{TC} 225";
  SampleNames["TC_WZ_300"]="\\rho_{TC} 300";
  SampleNames["TC_WZ_400"]="\\rho_{TC} 400";
  SampleNames["TC_WZ_500"]="\\rho_{TC} 500";

  Data.push_back(Sample("data"));
  CheckSamples(fin,Data);
  
  Bkg.push_back(Sample("WJetsToLNu", kOrange+3, 1, kOrange+10));
  vector<string> VV;
  VV.push_back("ZZTo4e");
  VV.push_back("ZZTo4mu");
  VV.push_back("ZZTo2e2mu");
  VV.push_back("GVJets");
  VV.push_back("WWTo2L2Nu");
  Bkg.push_back(Sample("VV", VV, kOrange+3, 1, kOrange+3));

  Bkg.push_back(Sample("TTJets"  , kRed+4, 1, kRed+2));
  
  vector<string> ZJets; 
  ZJets.push_back("DYJetsToLL");

  ZJets.push_back("ZBB0JetsToLNu");
  ZJets.push_back("ZBB1JetsToLNu");
  ZJets.push_back("ZBB2JetsToLNu");
  ZJets.push_back("ZBB3JetsToLNu");
  ZJets.push_back("ZCC0JetsToLNu");
  ZJets.push_back("ZCC1JetsToLNu");
  ZJets.push_back("ZCC2JetsToLNu");
  ZJets.push_back("ZCC3JetsToLNu");

  Bkg.push_back(Sample("ZJets", ZJets, kOrange+3, 1, kOrange+7));

  Bkg.push_back(Sample("WZTo3LNu"       , kOrange+3, 1, kOrange-2));
  
  CheckSamples(fin,Bkg);

  if(inName.find("Wprime") != string::npos){
    //Sig.push_back(Sample("WprimeToWZTo3LNu_M-300", 1, 1, kGreen));
    //Sig.push_back(Sample("WprimeToWZTo3LNu_M-400", 1, 1, 10));
    //Sig.push_back(Sample("WprimeToWZTo3LNu_M-500", 1, 1, 10));
    Sig.push_back(Sample("WprimeToWZTo3LNu_M-600", 1, 1, 10));
    //Sig.push_back(Sample("WprimeToWZTo3LNu_M-700", 1, 1, 10));
    //Sig.push_back(Sample("WprimeToWZTo3LNu_M-800", 1, 1, 10));
    //Sig.push_back(Sample("WprimeToWZTo3LNu_M-900", 1, 1, 10));

    //Sig.push_back(Sample("TC225",     2, 1, 10));
    //Sig.push_back(Sample("TC_WZ_300",     1, 1, kBlue));
    //Sig.push_back(Sample("TC_WZ_400",     1, 1, kBlue));
    //Sig.push_back(Sample("TC_WZ_500",     1, 1, kRed));
  }
  CheckSamples(fin,Sig);

  samples_.push_back(&Data);
  samples_.push_back(&Sig);
  samples_.push_back(&Bkg);
  if(debug) cout<<"Size of samples: "<<samples_.size()
                <<" Size of Data: "<<Data.size()
                <<" Size of Sig: "<<Sig.size()
                <<" Size of Bkg: "<<Bkg.size()<<endl;

  string efftitle[] = {"hNumEvts", "hEffAbs", "hEffRel"};

  vector<string> Cuts = GetListofCuts(fin);
  vector<string> variable; 
  
  if(inName.find("Wprime") != string::npos){
    variable.push_back("hWZMass");
    variable.push_back("hHt");
    variable.push_back("hWpt");
    variable.push_back("hZpt");
  }else if(inName.find("EWKWZ") != string::npos){
      variable.push_back("hZMass");
      variable.push_back("hWTransMass");
      variable.push_back("hMET");
  }else if(inName.find("HadVZ") != string::npos){
    variable.push_back("hVZMass");
    variable.push_back("hZMass");
    variable.push_back("hZpt");
    variable.push_back("hVMass");
    variable.push_back("hVpt");
  }

  if(opt.find("show") == string::npos){//Extra Plots
    if(!inName.find("Wprime") != string::npos){
      variable.push_back("hWZ3e0muMass");
      variable.push_back("hWZ2e1muMass");
      variable.push_back("hWZ1e2muMass");
      variable.push_back("hWZ0e3muMass");
      
      //variable.push_back("hWZTransMass");
      variable.push_back("hWZpt");
    
      variable.push_back("hQ");       
    }else if(inName.find("EWKWZ") != string::npos){
      variable.push_back("hEvtType");          
      variable.push_back("hEvtTypeP");          
      variable.push_back("hEvtTypeM");          
      
      variable.push_back("hZMass");      
      variable.push_back("hZeeMass");      
      variable.push_back("hZmmMass");  

      variable.push_back("hNJets");
      variable.push_back("hNVtxs");
      variable.push_back("hNLLeps");
    }else if(inName.find("HadVZ") != string::npos){
      variable.push_back("hVZpt");
      variable.push_back("hNLJets");
    }

  TCanvas c1;
  
  c1.Print(Form("%s[", outName.c_str()), "pdf"); 

  if(opt.find("show") == string::npos){
    for(int i=0;i<1;++i) DrawandSave(fin,outName,efftitle[i],"Title: "+efftitle[i], i==0, 1, 0);
  }
  uint begin = 0;
  if(opt.find("final") != string::npos || opt.find("show") != string::npos){
    begin = Cuts.size()-3;
  }else{
    for(uint i=0; i<Cuts.size(); ++i) if(!Cuts[i].compare("ZMass")) begin = i;
  }
    
  for(uint i=0;i<variable.size();++i){
    for(uint j=begin;j<Cuts.size();++j){
      string title = variable[i] + "_" + Cuts[j];
      string bkmark = (j==begin) ? "Title: " + variable[i] : "";

      bool log = 1;
      DrawandSave(fin,outName,title,bkmark,log,0,0);

      if(opt.find("show") != string::npos) 
        DrawandSave(fin,outName,title,bkmark,log,0,1);

    }
  }

  c1.Print(Form("%s]", outName.c_str()), "pdf"); 
}

void
DrawandSave(TFile* fin, string pdfName, string title, string bookmark, bool logy, bool eff, bool cum, TLine* line){
  GetHistograms(fin, title, eff, cum);
  string filename = logy ? "plots/" + title + ".pdf" : "plots/" + title + "_Linear.pdf";
  Draw(filename, pdfName, bookmark, logy, eff, line);
}

void
GetHistograms(TFile* fin, string title, bool eff, bool cum){
  if(debug) printf("  GetHisto %s\n", title.c_str());
  string hist_name;
  
  int rebin = (title.find("WZMass") != string::npos || 
               title.find("WZ3e0muMass") != string::npos ||
               title.find("WZ2e1muMass") != string::npos ||
               title.find("WZ1e2muMass") != string::npos ||
               title.find("WZ0e3muMass") != string::npos) ? 10 : 0;

  for(unsigned int i=0; i<samples_.size(); ++i){
    for(unsigned int j=0; j<samples_[i]->size(); ++j){
      if(debug) cout<<"i: "<<i<<" j:"<<j<<endl;
      samples_[i]->at(j).hist = get_sum_of_hists(fin, samples_[i]->at(j).names, title, rebin);
      samples_[i]->at(j).hist->SetLineStyle(samples_[i]->at(j).style);
      samples_[i]->at(j).hist->SetLineColor(samples_[i]->at(j).line); 
      if(!eff){
        samples_[i]->at(j).hist->SetFillColor(samples_[i]->at(j).fill);
        if(cum){
          samples_[i]->at(j).hist = FillCum(samples_[i]->at(j).hist);
        }
      }
    }
  }
}

void
Draw(string filename, string pdfName, string bookmark, bool logy, bool eff, TLine* line){
  if(debug) printf("  Draw %s\n", filename.c_str());

  TCanvas c1;
  if(logy) c1.SetLogy();

  string title = "";
  if(Data.size()){
    if(0) title = Data[0].hist->GetTitle();
    title += ";"; 
    title += Data[0].hist->GetXaxis()->GetTitle();
    title += ";";
    title += Data[0].hist->GetYaxis()->GetTitle();
  }else if(Bkg.size()){
    if(0) title = Bkg[0].hist->GetTitle();
    title += ";"; 
    title += Bkg[0].hist->GetXaxis()->GetTitle();
    title += ";";
    title += Bkg[0].hist->GetYaxis()->GetTitle();
  }
  TH1F* hData = NULL;
  for(uint i=0; i<Data.size(); i++){
    if(i==0) hData = (TH1F*)Data[i].hist->Clone("hData");
    else     hData->Add(Data[i].hist);
  }
  if(hData) hData->SetMarkerStyle(20);

  if(!eff){
    THStack* sBkg  = new THStack("sBkg",title.c_str());
    std::vector<THStack*> sSigs(Sig.size(),NULL);

    for(uint i=0; i<Bkg.size(); i++){
      sBkg->Add(Bkg[i].hist);
    } 
    
    for(uint i=0; i<Sig.size(); i++){
      sSigs[i] = (THStack*) sBkg->Clone();
      sSigs[i]->Add(Sig[i].hist);
    }
    
    Double_t max = sBkg->GetMaximum();
    sBkg->Draw("HIST");
    for(unsigned int i=0; i<Sig.size(); i++){
      max = TMath::Max(max, sSigs[i]->GetMaximum());
      sSigs[i]->Draw("HIST SAME");
    }
    if(hData){
      max = TMath::Max(max, hData->GetMaximum());
      hData->Draw("E1 SAME");
    }
    if(logy){
      sBkg->SetMaximum(50*max);
      sBkg->SetMinimum(0.01);
    }else{
      sBkg->SetMaximum(1.5*max);
      sBkg->SetMinimum(0.);
    }
    sBkg->GetXaxis()->SetTitleFont(132);
    sBkg->GetYaxis()->SetTitleFont(132);
    sBkg->GetXaxis()->SetTitleSize(0.06);
    sBkg->GetYaxis()->SetTitleSize(0.06);

  }else{
    THStack* hs = new THStack("hs",title.c_str());
    for(unsigned int i=0; i<Bkg.size(); ++i){
      hs->Add(Bkg[i].hist);
    }
    for(unsigned int i=0; i<Sig.size(); ++i){
      hs->Add(Sig[i].hist);
    }

    hs->Add(hData, "E1");
    
    hs->Draw("nostack");
    hs->GetXaxis()->SetTitleFont(132);
    hs->GetYaxis()->SetTitleFont(132);
    hs->GetXaxis()->SetTitleSize(0.06);
    hs->GetYaxis()->SetTitleSize(0.06);
  }
  if(debug) cout<<"Title: "<<title<<endl;

  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextSize(0.045);
  latexLabel.DrawLatex(0.18, 0.96, "#font[132]{CMS Preliminary 2011}");
  latexLabel.DrawLatex(0.5, 0.955, "#font[132]{#sqrt{s} = 7 TeV}");
  latexLabel.DrawLatex(0.73, 0.87, Form("#font[132]{#intL dt = %.2f fb^{-1}}",lumiUsed_/1000.));

  if(debug) cout<<"Creating Legend\n";
  TLegend *legend = new TLegend(0.62,0.69,0.92,0.85,"");
  
  if(Data.size()) legend->AddEntry(hData, "Data", "PLE");
  for(unsigned int i=0; i<Bkg.size(); ++i){
    legend->AddEntry(Bkg[i].hist,SampleNames[Bkg[i].name].c_str(), "F");
  }
  for(unsigned int i=0; i<Sig.size(); ++i){
    legend->AddEntry(Sig[i].hist,SampleNames[Sig[i].name].c_str(), "F");
  }
 
  legend->SetTextSize(0.05);

	legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetNColumns(2);
  legend->Draw();

  if(line) line->Draw();
  c1.RedrawAxis(); 

  c1.SaveAs(filename.c_str());
  //string bookmark = string Form("Title: %s",bookmark.c_str());
  c1.Print(pdfName.c_str(), bookmark.c_str());
}

void
CheckSamples(TFile* fin, vector<Sample> & sample){
  if(debug) cout<<"Checking sample "<<endl;
  for(size_t i=0; i<sample.size(); ++i){
    for(size_t j=0; j<sample[i].names.size(); ++j){
      if(debug) cout<<"Checking key "<<sample[i].names[j]<<endl;
      if(!fin->GetKey((sample[i].names[j]).c_str() )){
        cout<<"Didn't find "<<sample[i].names[j]<<". Removing."<<endl;
        sample.erase(sample.begin()+i);
        i--;
      }
    }
  }
}

TH1F* GetValidHist(TFile* f){
  string name;
  for(unsigned int i=0; i<samples_.size(); ++i){
    for(unsigned int j=0; j<samples_[i]->size(); ++j){
      name = samples_[i]->at(j).names[0] + "/hNumEvts";
      TH1F* hist = (TH1F*) f->Get(name.c_str());
      if(hist) return hist;
    }
  }
  cout<<"Didn't find a valid histogram\n";
  abort();
}

vector<string> GetListofCuts(TFile* f, TH1F* hist){
  if(!hist) hist = GetValidHist(f);
  if(!hist) return vector<string>();
  vector<string> Cuts(hist->GetXaxis()->GetNbins());
  for(uint i=0; i<Cuts.size(); ++i){
    Cuts[i] = hist->GetXaxis()->GetBinLabel(i+1);
  }
  return Cuts;
}

TH1F*
FillCum(TH1F* h){
  TH1F* cumhist = (TH1F*)h->Clone();
  string title = h->GetTitle();
  title += "(Cumlative)";

  int size = h->GetXaxis()->GetNbins();
  float sum = h->GetBinContent(size+1);//Overflow
  cumhist->SetTitle(title.c_str());
  for(int i=size; i>0; --i){
    sum += h->GetBinContent(i);
    cumhist->SetBinContent(i,sum);
    cumhist->SetBinError(i,sqrt(sum));
  }
  return cumhist;
}

int 
main(int argc, char ** argv){
  if(argc > 4){
    fprintf(stderr,"%s usage: %s inFile outFile <opt>\n",argv[0],argv[0]);
    exit( 1 );
  }
  MakePlots(argv[1], argv[2], argc > 3 ? argv[3] : "");

  return 0;
}

