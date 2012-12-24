/* 
   Usage: root -b -l -q 'MakePlots.cc+("input.root", "output.pdf", "options", <lumiWanted>)'
   Options include:
   ""    = every plot after every cut, note some may be blank
   "final" = show me all plots but only after the last cut
   "show" = show me only the main plots(eg for weekly update), only after the last cut
   "debug" = show debugging statements
   "lowSig" = show low xsec signal unstacked
   "paper" = Change style to official cms style
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
#include "setTDR_modified.C"
using namespace std;

struct Sample{
  std::string name;
  std::vector<std::string> names;
  std::vector<float> weights;
  std::vector<float> lumiUsed;
  int line;
  int style;
  int fill;
  TH1F* hist;
  Sample(){
    line = 0; style = 0; fill=0; hist=NULL;
  }
  Sample(std::string n, float scale=1.){
    name = n; names.push_back(n); 
    weights.push_back(scale); lumiUsed.push_back(0.);
    line = 1; style = 0; fill = 0; hist=NULL; 
  }
  Sample(std::string n, int l, int s, int f, float scale=1.){
    name = n; names.push_back(n); 
    weights.push_back(scale); lumiUsed.push_back(0.); 
    line = l; style = s; fill = f; hist=NULL;
  }
  Sample(std::string n, std::vector<std::string> ns, int l, int s, int f, float scale=1.){
    name = n; names = ns; 
    weights.assign(names.size(), scale); lumiUsed.assign(names.size(), 0.); 
    line = l; style = s; fill = f; hist=NULL;
  }
  void init(){
    line = 0; style = 0; fill=0; hist=NULL;
  }
};

enum Mode {kEWKWZ, kWprimeWZ, kHadVZ, kWprimeTB, kWprimeVW, kTTbar, kWZFakeRate};
Mode mode_;

std::vector< std::vector<Sample>* > samples_;
std::vector<Sample> Data;
std::vector<Sample> Bkg;
std::vector<Sample> Sig;
std::map<std::string, std::string> SampleNames;

bool debug_ = false;
bool drawLowSig_ = false;
float lumiUsed_ = 0.;
float lumiWanted_ = -1;
bool scaleLumi_ = false;
bool paperMode_ = false;
bool identifyPlot_ = false;

typedef std::vector<std::string> vstring;

//int main(int argc, char ** argv);
void MakePlots(string inName, string outName="", string opt="", float lumiWanted=-1);
void DrawandSave(TFile* fin, std::string pdfName, std::string title, std::string bookmark, bool logy=1, bool eff=0, bool cum=0, bool pull=0, TLine* line=NULL);
void DrawandSave(TFile* fin, std::string pdfName, vstring title, std::string bookmark, bool logy=1, bool eff=0, bool cum=0, bool pull=1, TLine* line=NULL);
TH1F* FillCum(TH1F* h);
bool GetHistograms(TFile* fin, std::string title, bool eff=0, bool cum=0);
void CheckSamples(TFile* fin, std::vector<Sample> & sample);
void ScaleSampleLumi(std::vector<Sample> & sample);
TH1F* GetValidHist(TFile* f);
vector<string> GetListofCuts(TFile* f, TH1F* hist=NULL);
void ChangeAxisRange(const string & filename, TAxis* xaxis);
TH1F* GetPullPlot(TH1F* hData, THStack * h);
string GetTitle();
TH1F* GetDataStack(string title);
THStack* GetBkgStack(string title);
vector<THStack*> GetSigStacks(string title, THStack* sBkg);
void DrawLabels();
void DrawLegend(TH1F* hData, bool eff=0);
int GetRebin(string title);
void FillNameMap();
 
void
MakePlots(string inName, string outName, string opt, float lumiWanted){  
  gErrorIgnoreLevel = kWarning;

  if(opt.find("debug") != string::npos) debug_ = true;
  if(opt.find("lowSig") != string::npos) drawLowSig_ = true;
  if(opt.find("paper") != string::npos) paperMode_ = true;
  if(opt.find("identify") != string::npos) identifyPlot_ = true;

  if(paperMode_) setTDRStyle();
  else{
    CMSstyle();
    gStyle->SetOptStat(0);
  }
  gROOT->ForceStyle();
  
  if(outName.empty()){
    size_t found = inName.find_last_of("/\\");
    outName = inName.substr(found+1);
    outName.replace(outName.find(".root"), 5, ".pdf"); 
    cout<<"Using output file "<<outName<<endl;
  }

  TFile *fin = TFile::Open(inName.c_str(), "read"); assert(fin);
  //Determine what analysis this is
  if     (inName.find("EWKWZ") != string::npos) mode_ = kEWKWZ;
  else if(inName.find("WprimeWZ") != string::npos) mode_ = kWprimeWZ;
  else if(inName.find("WprimeVW") != string::npos) mode_ = kWprimeVW;
  else if(inName.find("WprimeTB") != string::npos) mode_ = kWprimeTB;
  else if(inName.find("HadVZ") != string::npos) mode_ = kHadVZ;
  else if(inName.find("TTbar") != string::npos) mode_ = kTTbar;
  else if(inName.find("WZFakeRate") != string::npos) mode_ = kWZFakeRate;
  else{
    cerr<<" Don't know what you're trying to plot with input file. "
        <<inName<<endl;
    abort();
  }
  if(debug_) cout<<"Mode is determined to be "<<mode_<<endl;

  /////Data Samples
  if(mode_ == kWprimeWZ || mode_ == kEWKWZ ||
     mode_ == kWprimeVW || mode_ == kWprimeTB ||
     mode_ == kWZFakeRate){
    Data.push_back(Sample("data"));
  }else if(mode_ == kHadVZ){
    vector<string> vData;

    vData.push_back("data_DoubleMu-Run2011A-May10ReReco-v1");
    vData.push_back("data_SingleMu-Run2011A-May10ReReco-v1");
    vData.push_back("data_DoubleElectron-Run2011A-May10ReReco-v1");

    vData.push_back("data_DoubleMu-Run2011A-PromptReco-v4");
    vData.push_back("data_SingleMu-Run2011A-PromptReco-v4");
    vData.push_back("data_DoubleElectron-Run2011A-PromptReco-v4");

    vData.push_back("data_DoubleMu-Run2011A-05Aug2011-v1");
    vData.push_back("data_SingleMu-Run2011A-05Aug2011-v1");
    vData.push_back("data_DoubleElectron-Run2011A-05Aug2011-v1");

    vData.push_back("data_DoubleMu-Run2011A-PromptReco-v6");
    vData.push_back("data_SingleMu-Run2011A-PromptReco-v6");
    vData.push_back("data_DoubleElectron-Run2011A-PromptReco-v6");
    
    vData.push_back("data_SingleMu-Run2011B-PromptReco-v1");
    vData.push_back("data_DoubleMu-Run2011B-PromptReco-v1");
    vData.push_back("data_DoubleElectron-Run2011B-PromptReco-v1");
    
    Data.push_back(Sample("data", vData, 1, 0, 0));     
  }else if(mode_ == kTTbar){
    vector<string> vData;

    vData.push_back("data-DiLeptonJet-V360B-SingleMu-Run2011A-May10ReReco-v1.txt");
    vData.push_back("data-DiLeptonJet-V360B-SingleElectron-Run2011A-May10ReReco-v1.txt");

    vData.push_back("data-DiLeptonJet-V360B-SingleMu-Run2011A-PromptReco-v4.txt");
    vData.push_back("data-DiLeptonJet-V360B-SingleElectron-Run2011A-PromptReco-v4.txt");

    vData.push_back("data-DiLeptonJet-V360B-SingleMu-Run2011A-05Aug2011-v1.txt");
    vData.push_back("data-DiLeptonJet-V360B-SingleElectron-Run2011A-05Aug2011-v1.txt");

    vData.push_back("data-DiLeptonJet-V360B-SingleMu-Run2011A-PromptReco-v6.txt");
    vData.push_back("data-DiLeptonJet-V360B-SingleElectron-Run2011A-PromptReco-v6.txt");
    
    vData.push_back("data-DiLeptonJet-V360B-SingleMu-Run2011B-PromptReco-v1.txt");
    vData.push_back("data-DiLeptonJet-V360B-SingleElectron-Run2011B-PromptReco-v1.txt");
    
    Data.push_back(Sample("data", vData, 1, 0, 0));
  } 
  CheckSamples(fin,Data);
    
  /////Background Samples

  if(mode_ == kWprimeWZ || mode_ == kEWKWZ || mode_ == kWZFakeRate){
    vector<string> ZJets; 
    ZJets.push_back("DYJetsToLL");
    Bkg.push_back(Sample("ZJets", ZJets, kOrange+3, 1, kOrange+7));
    
    if(mode_ == kWZFakeRate || mode_ == kEWKWZ){
      //Bkg.push_back(Sample("WJetsToLNu", kOrange+6, 1, kOrange+0));
      Bkg.push_back(Sample("ZZ", kGreen+3, 1, kGreen+2));
      Bkg.push_back(Sample("GVJets", kOrange+3, 1, kOrange+5));
      Bkg.push_back(Sample("WWTo2L2Nu", kOrange+3, 1, kOrange+10));
    }else{
      vector<string> VV;
      VV.push_back("ZZ");
      VV.push_back("GVJets");
      //VV.push_back("WWTo2L2Nu");
      Bkg.push_back(Sample("VV", VV, kOrange+3, 1, kOrange+3));
    }
    Bkg.push_back(Sample("TTJets"  , kViolet+4, 1, kViolet+2));
    
    Bkg.push_back(Sample("WZJetsTo3LNu"       , kOrange+3, 1, kOrange-2));
  }else if(mode_ == kHadVZ || mode_ == kWprimeVW || mode_ == kWprimeTB){
    vector<string> VV;
    VV.push_back("Fall11-ZZ");
    VV.push_back("Fall11-VGamma");
    VV.push_back("Fall11-WW");
    VV.push_back("Fall11-WZ");
    Bkg.push_back(Sample("VV", VV, kGreen-7, 1, kGreen-7));
    
    Bkg.push_back(Sample("Fall11-TTJets"  , kBlue-7, 1, kBlue-7));
    
    vector<string> ZJets; 
    ZJets.push_back("Fall11-DYJetsToLL_PtZ100");
    Bkg.push_back(Sample("ZJets", ZJets, kRed-7, 1, kRed-7));
    
  }else if(mode_ == kTTbar){
    Bkg.push_back(Sample("TTJets_TuneZ2_7TeV-madgraph-tauola"  , kBlue-7, 1, kBlue-7));
    Bkg.push_back(Sample("WJetsToLNu_TuneZ2_7TeV-madgraph-tauola", kRed-7, 1, kRed-7));

  }
  CheckSamples(fin,Bkg);

  /////Signal Samples

  if(mode_ == kWprimeWZ){
    Sig.push_back(Sample("WprimeToWZTo3LNu_M-200", kBlue, 1, 0));
    Sig.push_back(Sample("WprimeToWZTo3LNu_M-500", kRed, 1, 0));
    Sig.push_back(Sample("WprimeToWZTo3LNu_M-1000", kGreen, 1, 0));
    Sig.push_back(Sample("WprimeToWZTo3LNu_M-1500", kOrange, 1, 0));
    Sig.push_back(Sample("WprimeToWZTo3LNu_M-2000", kCyan, 1, 0));
    //Sig.push_back(Sample("WprimeToWZTo3LNu_M-1900-MyGen", kGreen, 1, 0));
    //Sig.push_back(Sample("WprimeToWZTo3LNu_M-250", kRed, 1, 0));
    //Sig.push_back(Sample("WprimeToWZTo3LNu_M-300", kGreen, 1, 0));
    //Sig.push_back(Sample("WprimeToWZTo3LNu_M-400", 1, 1, 0));
    //Sig.push_back(Sample("WprimeToWZTo3LNu_M-500", 1, 1, 0));
    //////Sig.push_back(Sample("WprimeToWZTo3LNu_M-600", 1, 1, 0));
    //Sig.push_back(Sample("WprimeToWZTo3LNu_M-700", 1, 1, 0));
    //Sig.push_back(Sample("WprimeToWZTo3LNu_M-800", 1, 1, 0));
    //Sig.push_back(Sample("WprimeToWZTo3LNu_M-900", 1, kDashed, 0));

    //Sig.push_back(Sample("TC225",     2, 1, 0));
    //Sig.push_back(Sample("TC_WZ_300",     1, 1, kBlue));
    //Sig.push_back(Sample("TC_WZ_400",     1, 1, kBlue));
    //Sig.push_back(Sample("TC_WZ_500",     1, 1, kRed));
  }else if(mode_ == kHadVZ){
    vector<string> RS750;
    RS750.push_back("Summer11_RSZZeejj_750");
    RS750.push_back("Summer11_RSZZmmjj_750");
    //Bkg.push_back(Sample("RSZZ_750", RS750, kMagenta, 1, 0));
    //Sig.push_back(Sample("RSZZ_750", RS750, kMagenta+3, 1, 0)); 
    vector<string> RS1000;
    RS1000.push_back("Summer11_RSZZeejj_1000");
    RS1000.push_back("Summer11_RSZZmmjj_1000");
    //Bkg.push_back(Sample("RSZZ_1000", RS1000, kViolet, 1, 0));
    //    Sig.push_back(Sample("RSZZ_1000", RS1000, kViolet+4, 2, 0));

    vector<string> WPrime1000;
    WPrime1000.push_back("Summer11_WPrimeZZlljj_1000");
    //Sig.push_back(Sample("WPrimeWZ_1000", WPrime1000, kViolet+4, 2, 0));
 
    vector<string> RS1250;
    RS1250.push_back("Summer11_RSZZeejj_1250");
    RS1250.push_back("Summer11_RSZZmmjj_1250");
    //Bkg.push_back(Sample("RSZZ_1250", RS1250, kMagenta, 1, 0));
    vector<string> RS1500;
    RS1500.push_back("Summer11_RSZZeejj_1500");
    RS1500.push_back("Summer11_RSZZmmjj_1500");
    //Bkg.push_back(Sample("RSZZ_1500", RS1500, kMagenta, 1, 0));
    vector<string> RS1750;
    RS1750.push_back("Summer11_RSZZeejj_1750");
    RS1750.push_back("Summer11_RSZZmmjj_1750");
    //Bkg.push_back(Sample("RSZZ_1750", RS1750, kMagenta, 1, 0));
    vector<string> RS2000;
    RS2000.push_back("Summer11_RSZZeejj_2000");
    RS2000.push_back("Summer11_RSZZmmjj_2000");
    //Bkg.push_back(Sample("RSZZ_2000", RS2000, kCyan, 1, 0));
    // Sig.push_back(Sample("Summer11_RSZZmmjj_750",     kMagenta, 1, 0));
    // Sig.push_back(Sample("Summer11_RSZZmmjj_1000",     kMagenta, 2, 0));
    //Sig.push_back(Sample("Summer11_RSZZmmjj_1250",     1, 1, 0));
    //Sig.push_back(Sample("Summer11_RSZZmmjj_1500",     1, 1, 0));
    //Sig.push_back(Sample("Summer11_RSZZmmjj_1750",     1, 1, 0));
    //Sig.push_back(Sample("Summer11_RSZZmmjj_2000",     1, 1, 0));
    //Sig.push_back(Sample("Summer11_WprimeToWZTo2Q2L_M-500", kGray, 1, 0));
    //Sig.push_back(Sample("Summer11_WprimeToWZTo2Q2L_M-1000", kGray, 2, 0));
  }else if(mode_ == kWprimeVW){
  }else if(mode_ == kWprimeTB){
    Sig.push_back(Sample("WprimeTB_M-800", 1, 1, 0));
    Sig.push_back(Sample("WprimeTB_M-1000", kCyan, 1, 0));
    Sig.push_back(Sample("WprimeTB_M-1200", kPink, 1, 0));
    //Sig.push_back(Sample("WprimeTB_M-1500", 1, 1, 0));
    Sig.push_back(Sample("WprimeTB_M-2000", kMagenta, 1, 0));
  }
  CheckSamples(fin,Sig);

  //cout<<"New method for lumi is "<<lumiUsed_<<endl;
  //lumiUsed_ = GetLumiUsed(fin);//Replace by per sample lumi
  lumiWanted_ = lumiWanted > 0 ? lumiWanted : lumiUsed_;
  ScaleSampleLumi(Sig);
  ScaleSampleLumi(Bkg);

  cout<<"Lumi Used was "<<lumiUsed_<<" inv pb\n";
  cout<<"Lumi Wanted is "<<lumiWanted_<<" inv pb\n";

  FillNameMap();

  //Put all the samples into a single structure
  samples_.push_back(&Data);
  samples_.push_back(&Sig);
  samples_.push_back(&Bkg);
  if(debug_) cout<<"Size of samples: "<<samples_.size()
                 <<" Size of Data: "<<Data.size()
                 <<" Size of Sig: "<<Sig.size()
                 <<" Size of Bkg: "<<Bkg.size()
                 <<endl;


  //Cuts plotted after every cut defined here
  vector<string> Cuts = GetListofCuts(fin);
  vector<string> variable; 

  //These variables will be plotted after each cut always
  if(mode_ == kWprimeWZ){
    variable.push_back("hWZMass");
  }else if(mode_ == kEWKWZ){
    variable.push_back("hZMass");
    variable.push_back("hWTransMass");
    variable.push_back("hMET");
  }else if(mode_ == kHadVZ){
    variable.push_back("hVZMass");
    variable.push_back("hVZeeMass");
    variable.push_back("hVZmmMass");
  }else if(mode_ == kWZFakeRate){
    variable.push_back("hNJets");
    variable.push_back("hWTransMass");
    variable.push_back("hDeltaWMT");
    variable.push_back("hDeltaPhiWJet");
    variable.push_back("hDeltaPhiLepJet");
    variable.push_back("hDeltaRLepJet");
    variable.push_back("hMET");
    variable.push_back("hLeadElecPt");
    variable.push_back("hLeadMuonPt");
  }

  //These variables will be plotted after each cut unless you say "show"
  //They are cross checks but not critical to be shown
  if(opt.find("show") == string::npos){
    if(mode_ == kWprimeWZ){
      variable.push_back("hWZ3e0mMass");
      variable.push_back("hWZ2e1mMass");
      variable.push_back("hWZ1e2mMass");
      variable.push_back("hWZ0e3mMass");
      
      variable.push_back("hWZTransMass");
      variable.push_back("hWZpt");
    
      variable.push_back("hLt");
      variable.push_back("hQ");       

      variable.push_back("hZMass");
      variable.push_back("hZeeMass");
      variable.push_back("hZmmMass");

      variable.push_back("hZpt");
      variable.push_back("hZeept");
      variable.push_back("hZmmpt");

      variable.push_back("hWTransMass");
      variable.push_back("hWenuTransMass");
      variable.push_back("hWmnuTransMass");

      variable.push_back("hWpt");
      variable.push_back("hWenupt");
      variable.push_back("hWmnupt");

      variable.push_back("hMET");
      variable.push_back("hMETSig");
      variable.push_back("hMETee");
      variable.push_back("hMETmm");

      variable.push_back("hZeeDr");
      variable.push_back("hZmmDr");

      variable.push_back("hEvtType");
      variable.push_back("hEvtTypeP");
      variable.push_back("hEvtTypeM");

      variable.push_back("hWenuCombRelIso");
      variable.push_back("hWmnuCombRelIso");

      variable.push_back("hLeadPt");
      variable.push_back("hLeadElecPt");
      variable.push_back("hLeadMuonPt");

      variable.push_back("hNJets");

      variable.push_back("hNVtxs");
      variable.push_back("hWeight");

      variable.push_back("hNLLeps");

      variable.push_back("hLeadPt");
      variable.push_back("hLeadPtZee");
      variable.push_back("hLeadPtZmm");

      variable.push_back("hTriLepMass");
      variable.push_back("hL1FastJet");

    }else if(mode_ == kEWKWZ){
      variable.push_back("hEvtType");          
      variable.push_back("hEvtTypeP");          
      variable.push_back("hEvtTypeM");          
      
      variable.push_back("hZMass");      
      variable.push_back("hZeeMass");      
      variable.push_back("hZmmMass");  

      variable.push_back("hWTransMass");
      variable.push_back("hWenuTransMass");
      variable.push_back("hWmnuTransMass");

      variable.push_back("hNJets");
      variable.push_back("hNVtxs");
      variable.push_back("hNLLeps");
    }else if(mode_ == kHadVZ){
      variable.push_back("hVZpt");

      variable.push_back("hZMass");
      variable.push_back("hZpt");

      variable.push_back("heeVMass");
      variable.push_back("heeVpt");

      variable.push_back("hmmVMass");
      variable.push_back("hmmVpt");


      variable.push_back("hNLJets");
      variable.push_back("hNLLeps");
    }else if(mode_ == kWprimeVW){
      variable.push_back("hVWMass");
      //Cory:variable.push_back("hVWenuMass");
      //Cory:variable.push_back("hVWmnuMass");
      variable.push_back("hQ");
      variable.push_back("hVMass");
      variable.push_back("hVpt");
      variable.push_back("hMET");
      //variable.push_back("hMETSig");//Cory
      variable.push_back("hWTransMass");
      variable.push_back("hWenuTransMass");
      variable.push_back("hWmnuTransMass");
      variable.push_back("hWpt");
      variable.push_back("hWenupt");
      variable.push_back("hWmupt");//Cory: typo, fix
      variable.push_back("hNLLeps");
      variable.push_back("hNJets");
      variable.push_back("hNVtxs");
      //variable.push_back("hWeight");
    }else if(mode_ == kWprimeTB){
      variable.push_back("hTBMass");
      variable.push_back("hTBenuMass");
      variable.push_back("hTBmnuMass");
      variable.push_back("hQ");
      variable.push_back("hTMass");
      variable.push_back("hTenuMass");
      variable.push_back("hTmnuMass");
      variable.push_back("hEvtType");
      variable.push_back("hTpt");
      variable.push_back("hMET");
      variable.push_back("hMETSig");
      variable.push_back("hWTransMass");
      variable.push_back("hWenuTransMass");
      variable.push_back("hWmnuTransMass");
      variable.push_back("hWpt");
      //variable.push_back("hWenupt");
      //variable.push_back("hWmnupt");
      variable.push_back("hB1Mass");
      variable.push_back("hB1Disc");
      variable.push_back("hB1pt");
      variable.push_back("hB2Mass");
      variable.push_back("hB2Disc");
      variable.push_back("hB2pt");
      variable.push_back("hMbl");
      variable.push_back("hNLLeps");
      variable.push_back("hNLJets");
      variable.push_back("hNLBJets");
      variable.push_back("hNTBJets");
      variable.push_back("hWeight");
    }else if(mode_ == kTTbar){
      variable.push_back("hTMass");
      variable.push_back("hT2Mass");

      variable.push_back("hTpt");
      variable.push_back("hT2pt");

      variable.push_back("hMET");
      variable.push_back("hWTransMass");
      variable.push_back("hW2Mass");

      variable.push_back("hWpt");
      variable.push_back("hW2pt");
    }
  }

  TCanvas c1;
  c1.Print(Form("%s[", outName.c_str()), "pdf"); //Open .pdf

  //Plot the # of Evts plot
  string efftitle[] = {"hNumEvts"};
  if(opt.find("show") == string::npos){
    //This is the # of evts after each cut (stacked)
    DrawandSave(fin,outName,efftitle[0],"Title: "+efftitle[0], 1, 0);
    //This is the # of evts after each cut (buggy)
    //DrawandSave(fin,outName,efftitle[0],"Title: "+efftitle[0], 1, 1);
  }

  //Determine where to start showing plots
  uint begin = 0;
  if(opt.find("final") != string::npos || opt.find("show") != string::npos){
    begin = Cuts.size()-1;
  }else{
    for(uint i=0; i<Cuts.size(); ++i) if(!Cuts[i].compare("ZMass")) begin = i;
  }
    
  //Make plots for multiple cut stages
  for(uint i=0;i<variable.size();++i){
    for(uint j=begin;j<Cuts.size();++j){
      string title = variable[i] + "_" + Cuts[j];
      string bkmark = "Title: " + title;

      bool log = 1;
      DrawandSave(fin,outName,title,bkmark,log,0,0);

      if(opt.find("show") != string::npos) 
        DrawandSave(fin,outName,title,bkmark,log,0,1);

    }
  }

  //Make plots for a single cut
  if(mode_ == kWprimeWZ){
    if(opt.find("show") == string::npos) {
      DrawandSave(fin,outName,"hZeeMass_Lt","Title: Z Mass After Lt Zee",0);
      DrawandSave(fin,outName,"hZmmMass_Lt","Title: Z Mass After Lt Zmm",0);

      DrawandSave(fin,outName,"hZMass_ValidZ","Title: Z Mass After Valid Z",1);
      DrawandSave(fin,outName,"hZeeMass_ValidZ","Title: Z Mass After Valid Zee",1);
      DrawandSave(fin,outName,"hZmmMass_ValidZ","Title: Z Mass After Valid Zmm",1);

      DrawandSave(fin,outName,"hZMass_ValidZ","Title: Z Mass After Valid Z",0);
      DrawandSave(fin,outName,"hZeeMass_ValidZ","Title: Z Mass After Valid Zee",0);
      DrawandSave(fin,outName,"hZmmMass_ValidZ","Title: Z Mass After Valid Zmm",0);

      DrawandSave(fin,outName,"hWTransMass_ValidW","Title: W TMass After Valid W",1);
      DrawandSave(fin,outName,"hWenuTransMass_ValidW","Title: W TMass After Valid Wen",1);
      DrawandSave(fin,outName,"hWmnuTransMass_ValidW","Title: W TMass After Valid Wmu",1);

      //Adding bunch of MET plots
      DrawandSave(fin,outName,"hMET_ValidW","Title: MET After Valid W",1);
      DrawandSave(fin,outName,"hMETee_ValidW","Title: MET After Valid W (zee)",1);
      DrawandSave(fin,outName,"hMETmm_ValidW","Title: MET After Valid W (Zmm)",1);

      DrawandSave(fin,outName,"hMET_ValidW","Title: MET After Valid W",1,0,1);
      DrawandSave(fin,outName,"hMETee_ValidW","Title: MET After Valid W (zee)",1,0,1);
      DrawandSave(fin,outName,"hMETmm_ValidW","Title: MET After Valid W (Zmm)",1,0,1);

      DrawandSave(fin,outName,"hMET_ValidW","Title: MET After Valid W",0);
      DrawandSave(fin,outName,"hMETee_ValidW","Title: MET After Valid W (zee)",0);
      DrawandSave(fin,outName,"hMETmm_ValidW","Title: MET After Valid W (Zmm)",0);

      DrawandSave(fin,outName,"hMET_ValidW","Title: MET After Valid W",0,0,1);
      DrawandSave(fin,outName,"hMETee_ValidW","Title: MET After Valid W (zee)",0,0,1);
      DrawandSave(fin,outName,"hMETmm_ValidW","Title: MET After Valid W (Zmm)",0,0,1);
      //End MET

      DrawandSave(fin,outName,"hLt_MET","Title: Cumlative Lt before Lt Cut",1,0,1);

      DrawandSave(fin,outName,"hWZMass_MET","Title: WZ Mass before Lt Cut",1);
      DrawandSave(fin,outName,"hWZMass_MET","Title: WZ Mass before Lt Cut",1,0,1);
      DrawandSave(fin,outName,"hWZMass_Lt","Title: WZ Mass After Lt Cut ",1,0,0);
      DrawandSave(fin,outName,"hWZMass_Lt","Title: WZ Mass After Lt Cut (Cumlative) ",1,0,1);

      ////////////////
      DrawandSave(fin,outName,"hNLLeps_MET","Title: NLepAfter MET",1);
      DrawandSave(fin,outName,"hNJets_MET","Title: NJet After MET",1);
      DrawandSave(fin,outName,"hZeept_MET","Title: ZeePt After MET",1);
      DrawandSave(fin,outName,"hZmmpt_MET","Title: ZmmPt After MET",1);
      DrawandSave(fin,outName,"hNVtxs_MET","Title: NVtx After MET",1);
      DrawandSave(fin,outName,"hLeadPt_MET","Title: Lead Pt After MET",1);

      DrawandSave(fin,outName,"hLt_MET","Title: Lt After MET",1);
      DrawandSave(fin,outName,"hZpt_MET","Title: Zpt After MET",1);
      DrawandSave(fin,outName,"hWpt_MET","Title: Wpt After MET",1);
      
      DrawandSave(fin,outName,"hWZ3e0mMass_MET","Title: WZ 3e0mu Mass before Lt Cut",1,0,0);
      DrawandSave(fin,outName,"hWZ2e1mMass_MET","Title: WZ 2e1mu Mass before Lt Cut",1,0,0);
      DrawandSave(fin,outName,"hWZ1e2mMass_MET","Title: WZ 1e2mu Mass before Lt Cut",1,0,0);
      DrawandSave(fin,outName,"hWZ0e3mMass_MET","Title: WZ 0e3mu Mass before Lt Cut",1,0,0);
/*
      DrawandSave(fin,outName,"hDiscriminant","Title: Disc ",1,0,0);
      DrawandSave(fin,outName,"hDiscriminantFrac","Title: Disc Frac",1,0,0);
      DrawandSave(fin,outName,"hDiscriminantAngle","Title: Disc Angle",1,0,0);
      DrawandSave(fin,outName,"hDiscriminantReal","Title: Disc Real",1,0,0);
      DrawandSave(fin,outName,"hDiscriminantImag","Title: Disc Imag",1,0,0);
*/     
      vector<string> WZMassChannels;
      WZMassChannels.push_back("hWZ3e0mMass_MET");
      WZMassChannels.push_back("hWZ2e1mMass_MET");
      WZMassChannels.push_back("hWZ1e2mMass_MET");
      WZMassChannels.push_back("hWZ0e3mMass_MET");
      //DrawandSave(fin,outName,WZMassChannels,"Title: WZ Mass By Channel",1);
      DrawandSave(fin,outName,WZMassChannels,"Title: WZ Mass By Channel",1,0,0,0);
      //DrawandSave(fin,outName,WZMassChannels,"Title: WZ Mass By Channel",0);
      //DrawandSave(fin,outName,WZMassChannels,"Title: WZ Mass By Channel",0,0,0,0);
    }
  }else if(mode_ == kEWKWZ){
      vector<string> NVtxs2Channels;
      NVtxs2Channels.push_back("hNVtxsZee_ValidZ");
      NVtxs2Channels.push_back("hNVtxsZmm_ValidZ");
      DrawandSave(fin,outName,NVtxs2Channels,"Title: NVtx By Channel",0);

      vector<string> ZMass2Channels_ValidZ;
      ZMass2Channels_ValidZ.push_back("hZeeMass_ValidZ");
      ZMass2Channels_ValidZ.push_back("hZmmMass_ValidZ");
      DrawandSave(fin,outName,ZMass2Channels_ValidZ,"Title: Z Mass By Channel After Valid Z",0);

      vector<string> ZMass2Channels;
      ZMass2Channels.push_back("hZeeMass_ValidW");
      ZMass2Channels.push_back("hZmmMass_ValidW");
      DrawandSave(fin,outName,ZMass2Channels,"Title: Z Mass 2 By Channel",0);

      vector<string> ZMass4Channels;
      ZMass4Channels.push_back("hZ3e0mMass_ValidW");
      ZMass4Channels.push_back("hZ2e1mMass_ValidW");
      ZMass4Channels.push_back("hZ1e2mMass_ValidW");
      ZMass4Channels.push_back("hZ0e3mMass_ValidW");
      DrawandSave(fin,outName,ZMass4Channels,"Title: Z Mass By 4 Channels",0);

      vector<string> MET2Channels;
      MET2Channels.push_back("hMETee_ValidW");
      MET2Channels.push_back("hMETmm_ValidW");
      DrawandSave(fin,outName,MET2Channels,"Title: MET By 2 Channel",0);

      vector<string> MET4Channels;
      MET4Channels.push_back("hMET3e0m_ValidW");
      MET4Channels.push_back("hMET2e1m_ValidW");
      MET4Channels.push_back("hMET1e2m_ValidW");
      MET4Channels.push_back("hMET0e3m_ValidW");
      DrawandSave(fin,outName,MET4Channels,"Title: MET By 4 Channel",0, 0, 0, 0);

      vector<string> WTMass2Channels;
      WTMass2Channels.push_back("hWenuTransMass_ValidW");
      WTMass2Channels.push_back("hWmnuTransMass_ValidW");
      DrawandSave(fin,outName,WTMass2Channels,"Title: W TransMass By 2 Channel",0);

      vector<string> WTMass4Channels;
      WTMass4Channels.push_back("hW3e0mTransMass_ValidW");
      WTMass4Channels.push_back("hW2e1mTransMass_ValidW");
      WTMass4Channels.push_back("hW1e2mTransMass_ValidW");
      WTMass4Channels.push_back("hW0e3mTransMass_ValidW");
      DrawandSave(fin,outName,WTMass4Channels,"Title: W TransMass By 4 Channel",0);

      //DrawandSave(fin,outName,"TREEWLepPtCUTSweight*(ZTightCode==3)*(MET<20)*(WTransMass<20.)*(EvtType==2)*(abs(WLepEta)<=1.5)","Title: Obj Faking Barrel Elec Pt Loose",0);
      //DrawandSave(fin,outName,"TREEWLepPtCUTSweight*(ZTightCode==3)*(WTightCode==1)*(MET<20)*(WTransMass<20.)*(EvtType==2)*(abs(WLepEta)<=1.5)","Title: Obj Faking Barrel Elec Pt Tight",0);
      //DrawandSave(fin,outName,"TREEWLepPtCUTSweight*(ZTightCode==3)*(MET<20)*(WTransMass<20.)*(EvtType==2)*(abs(WLepEta)>1.5)","Title: Obj Faking Endcap Elec Pt Loose",0);
      //DrawandSave(fin,outName,"TREEWLepPtCUTSweight*(ZTightCode==3)*(WTightCode==1)*(MET<20)*(WTransMass<20.)*(EvtType==2)*(abs(WLepEta)>1.5)","Title: Obj Faking Endcap Elec Pt Tight",0);
      //DrawandSave(fin,outName,"TREEWLepPtCUTSweight*(ZTightCode==3)*(MET<20)*(WTransMass<20.)*(EvtType==1)","Title: Obj Faking Muon Pt Loose",0);
      //DrawandSave(fin,outName,"TREEWLepPtCUTSweight*(ZTightCode==3)*(WTightCode==1)*(MET<20)*(WTransMass<20.)*(EvtType==1)","Title: Obj Faking Muon Pt Tight",0);
      

  }else if(mode_ == kHadVZ){
    if(opt.find("show") == string::npos) {
      DrawandSave(fin,outName,"h_bestmass","Title: Best Mass",1);
      DrawandSave(fin,outName,"heeVMass_ValidV","Title: Leading Jet Mass",1);
      DrawandSave(fin,outName,"hmmVMass_ValidV","Title: Leading Jet Mass",1);
      vector<string> VZMassChannels;
      VZMassChannels.push_back("hVZeeMass_ValidVZ");
      VZMassChannels.push_back("hVZmmMass_ValidVZ");
      DrawandSave(fin,outName,VZMassChannels,"Title: VZ Mass By Channel",1);

      DrawandSave(fin,outName,"h_deltaR_elec1elec2","Title: DeltaR ee",1);
      DrawandSave(fin,outName,"h_deltaR_muon1muon2","Title: DeltaR mm",1);

      DrawandSave(fin,outName,"h_deltaR_HadVelec1","Title: DeltaR Ve1",1);
      DrawandSave(fin,outName,"h_deltaR_HadVelec2","Title: DeltaR Ve2",1);
      DrawandSave(fin,outName,"h_deltaR_HadVmuon1","Title: DeltaR Vm1",1);
      DrawandSave(fin,outName,"h_deltaR_HadVmuon2","Title: DeltaR Vm2",1);

      DrawandSave(fin,outName,"h_deltaR_jet1Z_R1","Title: DeltaR Zj1",1);
      DrawandSave(fin,outName,"h_deltaR_jet1Z_R2","Title: DeltaR Zj2",1);
      
    }
    
  }

  c1.Print(Form("%s]", outName.c_str()), "pdf"); //Close .pdf

}

void
DrawandSave(TFile* fin, string pdfName, string title, string bookmark, bool logy, bool eff, bool cum, bool pull, TLine* line){
  vstring titles(1,title);
  DrawandSave(fin, pdfName, titles, bookmark, logy, eff, cum, pull, line);
}

void
DrawandSave(TFile* fin, string pdfName, vstring title, string bookmark, bool logy, bool eff, bool cum, bool pull, TLine* line){
  string filename = "plots/" + title[0];
  if(title.size()>1) filename += "_Channels";
  if(!logy) filename += "_Linear"; 
  if(cum)  filename += "_Cumlative";
  filename += ".pdf";
  if(debug_) cout<<"In DrawandSave with filename "<<filename
                 <<" and "<<title.size()<<" sub plots"<<endl;
  
  TCanvas* canvas = new TCanvas();
  int nsubplots = title.size();
  vector<TH1F*> vhData(nsubplots);
  vector<THStack*> vsBkg(nsubplots);
  vector<vector<THStack*> >vsSig(nsubplots);
  for(int i=0; i<nsubplots; ++i) vsSig[i].assign(Sig.size(),NULL);

  if(!eff){
    //Cory: do the subplots here
    int ny = sqrt(nsubplots);//truncate to int
    int nx = nsubplots / ny;
    canvas->Divide(nx,ny);
    
    for(int i=0; i<nsubplots; ++i){
      if(!GetHistograms(fin, title[i], eff, cum)) return;
      string histname = GetTitle();
      vhData[i] = GetDataStack(histname);
      vsBkg[i] = GetBkgStack(histname);
      vsSig[i] = GetSigStacks(histname, vsBkg[i]);
      TH1F* hData = vhData[i];
      THStack* sBkg = vsBkg[i];
      vector<THStack*> & sSigs = vsSig[i];

      TPad* subpad = (TPad*) canvas->cd(i+1);//0 is whole page
      TPad* pad = subpad;
      if(pull){ //Cory: do all the pull plot stuff here
        subpad->Divide(1,2);
        pad = (TPad*) subpad->cd(2);//Change to lower plot
        pad->SetPad(0.,0., 1.,0.20);
        TH1F* hpull = GetPullPlot(hData, sBkg);
        hpull->Draw("e");

        //This changes the range of the X axis
        TAxis* axis = hpull->GetXaxis();
        ChangeAxisRange(title[i], hpull->GetXaxis());
        float xmin = hpull->GetBinLowEdge(axis->GetFirst());
        float xmax = hpull->GetBinLowEdge(axis->GetLast()+1);
        TLine* baseline = new TLine(xmin,0,xmax,0);
        baseline->SetLineColor(kRed);
        baseline->Draw();
        pad->RedrawAxis(); 

        pad = (TPad*) subpad->cd(1);//Change to upper plot 
        pad->SetPad(0.,0.20, 1.,1.);
      }//end pull plots
      
      //Now Draw!!!!!
      pad->SetLogy(logy);

      //Find Maximum
      Double_t max = sBkg->GetMaximum();
      sBkg->Draw("HIST");

      for(uint isig=0; isig<Sig.size(); ++isig){
        max = TMath::Max(max, sSigs[isig]->GetMaximum());
        sSigs[isig]->Draw("HIST SAME");
      }
           
      if(drawLowSig_){//This is if you want to see the shape of a low xsec signal
        for(uint isig=0; isig<Sig.size(); ++isig){
          Sig[isig].hist->Draw("HIST SAME");
        }
      }
      
      if(hData){//Draw Data if it exists
        max = TMath::Max(max, hData->GetMaximum());
        hData->Draw("E SAME");

        DrawLabels();
        DrawLegend(hData, 0);

      }

      if(debug_) cout<<" max is "<<max<<" for logy = "<<logy<<endl;
      sBkg->SetMaximum(logy ? 100*max : 1.5*max);
      sBkg->SetMinimum(logy ? 0.5 : 0.);

      //This changes the range of the X axis
      ChangeAxisRange(title[i], sBkg->GetXaxis());

      if(line) line->Draw();
      pad->RedrawAxis(); 
      
    }//end of all subplots
    canvas->cd();
    //End subplots
    
    
  }else{//Eff
    if(!GetHistograms(fin, title[0], eff)) return;
    string histname = GetTitle();
    vhData[0] = GetDataStack(histname);
    vsBkg[0] = GetBkgStack(histname);
    vsSig[0] = GetSigStacks(histname, vsBkg[0]);
    TH1F* hData = vhData[0];
    THStack* sBkg = vsBkg[0];
    vector<THStack*> & sSigs = vsSig[0];

    Double_t max = sBkg->GetMaximum();
    sBkg->Draw("nostack"); 

    for(uint isig=0; isig<Sig.size(); ++isig){
      max = TMath::Max(max, sSigs[isig]->GetMaximum());
      sSigs[isig]->Draw("nostack same");
    }
    max = TMath::Max(max, hData->GetMaximum());
    hData->Draw("E same");
    sBkg->SetMaximum(logy ? 100*max : 1.5*max);
    sBkg->SetMinimum(logy ? 0.5 : 0.);
    
  }

  canvas->SaveAs(filename.c_str());
  if(paperMode_){
    filename.replace(filename.find(".pdf"), 4, ".png"); 
    canvas->SaveAs(filename.c_str());
  }

  //string bookmark = string Form("Title: %s",bookmark.c_str());
  canvas->Print(pdfName.c_str(), bookmark.c_str());
  //Draw(filename, pdfName, bookmark, logy, eff, pull, subplots, line);

  //Clean Up
  for(int i=0; i<nsubplots; ++i){
    delete vhData[i], vsBkg[i];
    vector<THStack*> & sSigs = vsSig[i];
    for(uint isig=0; isig<Sig.size(); ++isig){ 
      delete sSigs[isig];
    }
  }
  delete canvas;
}

bool
GetHistograms(TFile* fin, string title, bool eff, bool cum){
  if(debug_) printf("  GetHisto %s\n", title.c_str());
  //How many bins should I combine
  int rebin = GetRebin(title);

  bool validHist = false;//If none of the histograms are filled ,don't draw it
  for(unsigned int i=0; i<samples_.size(); ++i){//Loop over Data, Bkg, Sig
    vector<Sample> & samples = *samples_[i];
    for(unsigned int j=0; j<samples.size(); ++j){
      if(debug_) cout<<"i: "<<i<<" j:"<<j<<endl;
      Sample& curSample = samples.at(j);//Let's make this easy

      ///Add ability to use Data Driven Methods instead of MC//////
      vector<string> names = curSample.names;
      if(1 && mode_ == kEWKWZ){
        if(title.find("hMET3e0m_ValidW") != string::npos ||//met
           title.find("hMET2e1m_ValidW") != string::npos ||
           title.find("hMET1e2m_ValidW") != string::npos ||
           title.find("hMET0e3m_ValidW") != string::npos ||

           title.find("hWpt3e0m_ValidW") != string::npos ||//Wpt
           title.find("hWpt2e1m_ValidW") != string::npos ||
           title.find("hWpt1e2m_ValidW") != string::npos ||
           title.find("hWpt0e3m_ValidW") != string::npos ||

           title.find("hZMass3e0m_ValidW") != string::npos ||//zmass
           title.find("hZMass2e1m_ValidW") != string::npos ||
           title.find("hZMass1e2m_ValidW") != string::npos ||
           title.find("hZMass0e3m_ValidW") != string::npos ||

           title.find("hWTransMass3e0m_ValidW") != string::npos ||//wmt
           title.find("hWTransMass2e1m_ValidW") != string::npos ||
           title.find("hWTransMass1e2m_ValidW") != string::npos ||
           title.find("hWTransMass0e3m_ValidW") != string::npos ||

           title.find("hZpt3e0m_ValidW") != string::npos ||//zpt
           title.find("hZpt2e1m_ValidW") != string::npos ||
           title.find("hZpt1e2m_ValidW") != string::npos ||
           title.find("hZpt0e3m_ValidW") != string::npos ){
          //loop over names and replace MC dir with Data driven
          replace (names.begin(), names.end(), (string)"DYJetsToLL", (string)"DYJetsToLL-DataDriven");
          //curSample.name = "DataDrivenBkg";

          //remove ttbar
//           for(int p=0; p<(int)names.size(); ++p){
//             cout<<" p="<<p<<endl;
//             cout<<" names is "<<names[p]<<endl;
//             if(names[p].compare("TTJets") == 0){
//               curSample.names.erase(curSample.names.begin()+p);
//               p--;
//               if(curSample.names.size() == 0 ){//Sample is now empty
//                 samples.erase(samples.begin()+j);
//                 j--;
//               }
//               continue;
//             }
//           }//p loop
          //don't change title for other MC/data
        }//j loop
      }//i loop

      ///////////
      if(title.find("TREE") == string::npos){
        //cout<<"title is "<<title<<endl;
        curSample.hist = get_sum_of_hists(fin, names, title, rebin, curSample.weights);
      }else{//Get Histograms from Trees///
        curSample.hist = new TH1F("a", "", 10, 0, 50);
        string newtitle = title;
        newtitle.replace(newtitle.find("TREE"), 4, "");
        size_t pos = newtitle.find("CUTS");
        string var = newtitle.substr(0, pos);
        string cuts = newtitle.substr(pos);
        cuts.replace(cuts.find("CUTS"), 4, "");
        //cout<<" pos: "<<pos
        //<<" var: "<<var
        //<<" cuts: "<<cuts
        //<<endl;
        
        string obj = "tEvts_ValidW";
        TTree* t = getTree(fin, names, "tEvts_ValidW"); assert(t); 
        get_sum_of_hists(fin, names, obj, var, cuts, *curSample.hist, curSample.weights);
        assert(curSample.hist);
      }
      /////
      if(!validHist && curSample.hist->Integral() > 0)//Only count filled histos
        validHist = true;
      curSample.hist->SetLineStyle(curSample.style);
      curSample.hist->SetLineColor(curSample.line); 

      //if(samples_[i] != &Data) curSample.hist->Scale(lumiWanted_/lumiUsed_); //Don't scale data


      
      if(!eff){
        curSample.hist->SetFillColor(curSample.fill);
        if(cum){//Do you want a cumlative distrubution
          curSample.hist = FillCum(curSample.hist);
          string newtitle = "Cumlative ";
          newtitle += curSample.hist->GetYaxis()->GetTitle();
          curSample.hist->SetYTitle(newtitle.c_str());
        }
      }else{//Eff Histos
        if     (samples_[i] == &Sig) curSample.hist->SetMarkerStyle(kOpenCircle);//Sig
        else if(samples_[i] == &Bkg) curSample.hist->SetMarkerColor(curSample.fill); //Bkg
        
        if(debug_){
          cout<<" Sample: "<<curSample.name<<endl;
          TH1F* hist = curSample.hist;
          for(int iBin=1; iBin<=hist->GetNbinsX(); ++iBin){
            cout<<"\tCut: "<<hist->GetXaxis()->GetBinLabel(iBin)<<": "
                <<hist->GetBinContent(iBin)
                <<" +/- "<<hist->GetBinError(iBin)<<endl;
          }
        }
      }
    }
  }
  if(debug_) cout<<" validHist = "<<validHist<<endl;
  return validHist;
}

void
CheckSamples(TFile* fin, vector<Sample> & samples){
  if(debug_) cout<<"Checking sample "<<endl;
  for(size_t i=0; i<samples.size(); ++i){
    Sample & sample = samples[i];
    for(size_t j=0; j<sample.names.size(); ++j){
      string & name = sample.names[j];
      float & sampleLumiUsed = sample.lumiUsed[j];
      cout<<" looking at sample "<<name<<endl;

      if(debug_) cout<<"Checking key "<<name<<endl;
      TH1F* hInfo = (TH1F*) fin->Get(Form("%s/hFileInfo", name.c_str()));
      //Check if folder exists
      if(!fin->GetKey((name).c_str() ) || !hInfo){
        if(!fin->GetKey((name).c_str() )) 
          cout<<"Didn't find folder: "<<name<<". Removing."<<endl;
        else
          cout<<"Didn't find hFileInfo for sample "<<name<<". Removing."<<endl;
        sample.names.erase(sample.names.begin()+j);
        j--;
        if(sample.names.size() == 0 ){//Sample is now empty
          samples.erase(samples.begin()+i);
          i--;
        }
        continue;
      }

      sampleLumiUsed = GetSampleInfo(hInfo, "#intL dt (pb^{-1})");
      lumiUsed_ = max(lumiUsed_, sampleLumiUsed);//Cory: can't use lumi wanted here, not determined yet

      //Check if all subjobs finished
      int nJobsTotal = GetSampleInfo(hInfo, "Number of SubSamples");
      int nJobsDone  = GetSampleInfo(hInfo, "Number of Files Merged");
      if(nJobsTotal == 0){
        printf("Not trying to scale sample %s because I couldn't find bin \"Number of SubSamples\"\n", name.c_str());
        continue;//Cory: Didn't record this so skip for now.  Delete soon (2012-05-31)
      }
      if(nJobsDone != nJobsTotal){
        if(nJobsDone < nJobsTotal ) printf(" Only %i of %i jobs finished for %s.  Scaling to compensate\n",
                                          nJobsDone, nJobsTotal, name.c_str());
        else printf("Job totals for %s don't match (%i expected, see %i) and you're combined different job types so the results are junk\n",
                    name.c_str(), nJobsTotal, nJobsDone);
      }
      //Scale sample to compensate for missing jobs (except data)
      if(&samples != &Data) sample.weights[j] *= (float) nJobsTotal / nJobsDone;
      if(debug_) cout<<"Scaling set to "<<sample.weights[j]<<endl;

    }//subsample loop
  }//sample loop
}

void
ScaleSampleLumi(vector<Sample> & samples){
  if(debug_) cout<<"Scaling sample lumi "<<endl;
  for(size_t i=0; i<samples.size(); ++i){
    Sample & sample = samples[i];
    for(size_t j=0; j<sample.names.size(); ++j){
      string & name = sample.names[j];

      //Scale by lumi
      sample.weights[j] *= lumiWanted_ / sample.lumiUsed[j];//Cory: should i do something diff for data
      if(debug_) printf("Scaling sample %s by %.2f to equalize lumi\n", name.c_str(), lumiWanted_ / sample.lumiUsed[j]);
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

  int size = h->GetNbinsX();
  float sum = h->GetBinContent(size+1);//Overflow
  cumhist->SetTitle(title.c_str());
  for(int i=size; i>0; --i){
    sum += h->GetBinContent(i);
    cumhist->SetBinContent(i,sum);
    cumhist->SetBinError(i,sqrt(sum));
  }
  return cumhist;
}

void
ChangeAxisRange(const string & filename, TAxis* xaxis){
  if(debug_) cout<<"Changing axis range\n";
  //This fn changes the range of the X axis
  if (filename.find("hWZ3e0muMass_") != string::npos ||
      filename.find("hWZ2e1muMass_") != string::npos ||
      filename.find("hWZ1e2muMass_") != string::npos ||
      filename.find("hWZ0e3muMass_") != string::npos ||
      filename.find("hWZMass_") != string::npos){
    xaxis->SetRangeUser(0,1500);
    //hpull->SetAxisRange( -3., 3., "Y");
  }else if (filename.find("hMET_") != string::npos ||
            filename.find("hMETee_") != string::npos ||
            filename.find("hMETmm_") != string::npos ||
            filename.find("hMET3e0m_") != string::npos ||
            filename.find("hMET2e1m_") != string::npos ||
            filename.find("hMET1e2m_") != string::npos ||
            filename.find("hMET0e3m_") != string::npos){
    xaxis->SetRangeUser(0,300);
    //hpull->SetAxisRange( -3., 3., "Y");
  }else if (filename.find("hLeadMuonPt_") != string::npos ||
            filename.find("hLeadElecPt_") != string::npos){
    xaxis->SetRangeUser(0,500);
  }else if (filename.find("hZMass_") != string::npos ||
            filename.find("hZMass_") != string::npos){
    xaxis->SetRangeUser(70,110);
  }

  

}

TH1F*
GetPullPlot(TH1F* hData, THStack * h){//Actually a residue plot for now
  if(debug_) cout<<"Creating pull plot\n";
  string title = hData->GetTitle();
  title += " (Pull)";
  TH1F* hpull = (TH1F*)hData->Clone(title.c_str());
  
  int size = hData->GetNbinsX();
  hpull->SetTitle(";;#sigma(Data-MC)");
  //hpull->SetTitle(title.c_str());
  TH1F* hsum = (TH1F*) h->GetStack()->Last();
  for(int i=1; i<=size; ++i){
    float diff = hData->GetBinContent(i) - hsum->GetBinContent(i);
    //diff = (hsum->GetBinContent(i) > 0) ? diff/hData->GetBinError(i) : 0.;
    if(debug_) cout<<"diff for bin "<<i<<" is "<<diff<<endl;
    hpull->SetBinContent(i,diff);
    hpull->SetBinError(i,hData->GetBinError(i));//Cory: include mc error
  }
  //hpull->SetAxisRange( -3., 3., "Y");
  return hpull;
}

TH1F*
GetDataStack(string title){
  if(debug_) cout<<"Creating Data histo\n";
  TH1F* hData = NULL;
  for(uint i=0; i<Data.size(); i++){
    if(i==0) hData = (TH1F*)Data[i].hist->Clone("hData");
    else     hData->Add(Data[i].hist);
  }
  hData->SetTitle(title.c_str());
  return hData;
}

THStack*
GetBkgStack(string title){
  if(debug_) cout<<"Creating Bkg Stack with title: "<<title<<endl;
  THStack* sBkg = new THStack("sBkg",title.c_str());
  
  for(uint i=0; i<Bkg.size(); i++){
      sBkg->Add(Bkg[i].hist);
  } 
  return sBkg;
}

vector<THStack*>
GetSigStacks(string title, THStack* sBkg){
  if(debug_) cout<<"Creating Sig Stack with title: "<<title<<endl;
  vector<THStack*> sSigs(Sig.size(),NULL);
  for(uint isig=0; isig<Sig.size(); ++isig){
    sSigs[isig] = (THStack*) sBkg->Clone();
    sSigs[isig]->Add(Sig[isig].hist);
  }
  return sSigs;
}

string
GetTitle(){
  string title;
  vector<Sample> & samples = Data.size() ? Data : Bkg;
  if(identifyPlot_) title = samples[0].hist->GetTitle();
  title += ";"; 
  title += samples[0].hist->GetXaxis()->GetTitle();
  title += ";";
  title += samples[0].hist->GetYaxis()->GetTitle();

  if(debug_) cout<<"Title: "<<title<<endl;
  return title;
}

void
DrawLabels(){
  if(debug_) cout<<"Creating Labels\n";
  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextSize(0.05);
  latexLabel.SetTextFont(42);
  latexLabel.DrawLatex(0.33, 0.96, "CMS Preliminary 2012");
  latexLabel.DrawLatex(paperMode_ ? 0.20 : 0.25, 0.77, "#sqrt{s} = 8 TeV");
  latexLabel.DrawLatex(paperMode_ ? 0.16 : 0.20, 0.85, Form("#intL dt = %.1f fb^{-1}",lumiWanted_/1000.));
}


void
DrawLegend(TH1F* hData, bool eff){//Cory: is this a problem?
  if(debug_) cout<<"Creating Legend\n";
  float legMinX = paperMode_ ? 0.55 : 0.43;
  float legMinY = paperMode_ ? 0.65 : 0.68;
  TLegend *legend = new TLegend(legMinX,legMinY,0.91,0.92,"");
  
  if(Data.size()) legend->AddEntry(hData, "Data", "PE");
  for(unsigned int i=0; i<Bkg.size(); ++i){
    string legName = SampleNames.find(Bkg[i].name) != SampleNames.end() ? SampleNames[Bkg[i].name] : Bkg[i].name;
    legend->AddEntry(Bkg[i].hist,legName.c_str(), eff ? "P" : "F");
  }
  for(unsigned int i=0; i<Sig.size(); ++i){
    string legName = SampleNames.find(Sig[i].name) != SampleNames.end() ? SampleNames[Sig[i].name] : Sig[i].name;
    legend->AddEntry(Sig[i].hist,legName.c_str(), eff ? "P" : "FL");
  }
 
  legend->SetTextSize(0.05);
  legend->SetTextFont(42);
	legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetNColumns(paperMode_ ? 1 : 2);
  legend->SetColumnSeparation(0.05);
  legend->Draw();
}

int GetRebin(string title){
  int rebin = 0;
  if(title.find("WZMass") != string::npos || 
     title.find("WZ3e0mMass") != string::npos ||
     title.find("WZ2e1mMass") != string::npos ||
     title.find("WZ1e2mMass") != string::npos ||
     title.find("WZ0e3mMass") != string::npos) rebin = 5;
  if(title.find("VZMass") != string::npos ||
     title.find("VZeeMass") != string::npos ||
     title.find("VZmmMass") != string::npos) rebin = 2;

  if(title.find("eeVMass") != string::npos) rebin = 3;
  if(title.find("mmVMass") != string::npos) rebin = 3;

  if(title.find("hMET_") != string::npos) rebin = 1;
  if(title.find("hLt_") != string::npos) rebin = 1;
  if(title.find("hWZTransMass_") != string::npos) rebin = 2;
  if(title.find("hWZpt_") != string::npos) rebin = 2;
  
  if(title.find("hTBMass_") != string::npos ||
     title.find("hTBenuMass_") != string::npos ||
     title.find("hTBmnuMass_") != string::npos) rebin = 5;
  if(title.find("hVWMass_") != string::npos ||
     title.find("hVWenuMass_") != string::npos ||
     title.find("hVWmnuMass_") != string::npos) rebin = 5;
  if(title.find("hDeltaPhiWJet_") != string::npos ||
     title.find("hDeltaPhiLepJet_") != string::npos) rebin = 4;
  if(title.find("hDeltaRLepJet_") != string::npos) rebin = 5;

  return rebin;
}

void
FillNameMap(){
  SampleNames["data"]="Data";
  SampleNames["WJetsToLNu"]="W+Jets";
  SampleNames["WWTo2L2Nu"]="WW\\rightarrow2l2\\nu";
  SampleNames["TTJets"]="t\\bar{t}";
  SampleNames["Fall11-TTJets"]="t\\bar{t}";
  SampleNames["ZZTo4L_TuneZ2"]="ZZ\\rightarrow4l";
  SampleNames["PhotonVJets"]="V\\gamma";
  SampleNames["VV"]="ZZ/Z#gamma";//"VV";
  SampleNames["GVJets"]="V#gamma";//"VV";
  SampleNames["ZZ"]="ZZ";//"VV";
  SampleNames["DYJetsToLL"]="DY+Jets\\rightarrow2l";
  SampleNames["Fall11-DYJetsToLL"]="DY+Jets\\rightarrow2l";
  SampleNames["ZJets"]="Z+Jets";
  SampleNames["ZBB0JetsToLNu"]="Z+0Jets\\rightarrowbb";
  SampleNames["ZBB1JetsToLNu"]="Z+1Jets\\rightarrowbb";
  SampleNames["ZBB2JetsToLNu"]="Z+2Jets\\rightarrowbb";
  SampleNames["ZBB3JetsToLNu"]="Z+3Jets\\rightarrowbb";
  SampleNames["ZCC0JetsToLNu"]="Z+0Jets\\rightarrowcc";
  SampleNames["ZCC1JetsToLNu"]="Z+1Jets\\rightarrowcc";
  SampleNames["ZCC2JetsToLNu"]="Z+2Jets\\rightarrowcc";
  SampleNames["ZCC3JetsToLNu"]="Z+3Jets\\rightarrowcc";
  SampleNames["WZ"]="WZ";
  SampleNames["WZJetsTo3LNu"]="WZ";//\\rightarrow3l\\nu";
  SampleNames["WprimeToWZTo3LNu_M-200"]="W' (200 GeV)";
  SampleNames["WprimeToWZTo3LNu_M-250"]="W' (250 GeV)";
  SampleNames["WprimeToWZTo3LNu_M-300"]="W' (300 GeV)";
  SampleNames["WprimeToWZTo3LNu_M-400"]="W' (400 GeV)";
  SampleNames["WprimeToWZTo3LNu_M-500"]="W' (500 GeV)";
  SampleNames["WprimeToWZTo3LNu_M-600"]="W' (600 GeV)";
  SampleNames["WprimeToWZTo3LNu_M-700"]="W' (700 GeV)";
  SampleNames["WprimeToWZTo3LNu_M-800"]="W' (800 GeV)";
  SampleNames["WprimeToWZTo3LNu_M-900"]="W' (900 GeV)";
  SampleNames["WprimeToWZTo3LNu_M-1000"]="W' (1000 GeV)";
  SampleNames["WprimeToWZTo3LNu_M-1100"]="W' (1100 GeV)";
  SampleNames["WprimeToWZTo3LNu_M-1200"]="W' (1200 GeV)";
  SampleNames["WprimeToWZTo3LNu_M-1300"]="W' (1300 GeV)";
  SampleNames["WprimeToWZTo3LNu_M-1400"]="W' (1400 GeV)";
  SampleNames["WprimeToWZTo3LNu_M-1500"]="W' (1500 GeV)";
  SampleNames["WprimeToWZTo3LNu_M-1600"]="W' (1600 GeV)";
  SampleNames["WprimeToWZTo3LNu_M-1700"]="W' (1700 GeV)";
  SampleNames["WprimeToWZTo3LNu_M-1800"]="W' (1800 GeV)";
  SampleNames["WprimeToWZTo3LNu_M-1900"]="W' (1900 GeV)";
  SampleNames["WprimeToWZTo3LNu_M-2000"]="W' (2000 GeV)";
  SampleNames["TC225"]="\\rho_{TC} 225";
  SampleNames["TC_WZ_300"]="\\rho_{TC} 300";
  SampleNames["TC_WZ_400"]="\\rho_{TC} 400";
  SampleNames["TC_WZ_500"]="\\rho_{TC} 500";
  SampleNames["RSZZ_750"]="RS 750";
  SampleNames["RSZZ_1000"]="RS 1000";
  SampleNames["RSZZ_1250"]="RS 1250";
  SampleNames["RSZZ_1500"]="RS 1500";
  SampleNames["RSZZ_1750"]="RS 1750";
  SampleNames["RSZZ_2000"]="RS 2000";
  SampleNames["WPrimeWZ_1000"]="W'1000";
  SampleNames["Summer11_WprimeToWZTo2Q2L_M-500"]="W' 500";
  SampleNames["Summer11_WprimeToWZTo2Q2L_M-1000"]="W' 1000";
}
/*
int 
main(int argc, char ** argv){
  if(argc > 4){
    fprintf(stderr,"%s usage: %s inFile outFile <opt>\n",argv[0],argv[0]);
    exit( 1 );
  }
  MakePlots(argv[1], argv[2], argc > 3 ? argv[3] : "");

  return 0;
}
*/
