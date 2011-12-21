/* 
   Usage: root -b -q 'MakePlots.cc+("input.root", "output.pdf, "options")'
   Options include:
   ""    = every plot after every cut, note some may be blank
   "final" = show me all plots but only after the last cut
   "show" = show me only the main plots(eg for weekly update), only after the last cut
   "debug" = show debugging statements
   "lowSig" = show low xsec signal unstacked
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

bool debug_ = false;
bool drawLowSig_ = false;
float lumiUsed_ = 0.;

int main(int argc, char ** argv);
void MakePlots(string inName, string outName, string opt="");
void DrawandSave(TFile* fin, std::string pdfName, std::string title, std::string bookmark, bool logy, bool eff, bool cum, TLine* line=NULL);
TH1F* FillCum(TH1F* h);
bool GetHistograms(TFile* fin, std::string title, bool eff, bool cum);
void Draw(std::string filename, std::string pdfName, std::string bookmark, bool logy, bool eff, TLine* line);
void CheckSamples(TFile* fin, std::vector<Sample> & sample);
TH1F* GetValidHist(TFile* f);
vector<string> GetListofCuts(TFile* f, TH1F* hist=NULL);

void
MakePlots(string inName, string outName, string opt){  
  gErrorIgnoreLevel = kWarning;
  CMSstyle();

  if(opt.find("debug") != string::npos) debug_ = true;
  if(opt.find("lowSig") != string::npos) drawLowSig_ = true;
  
  TFile *fin = TFile::Open(inName.c_str(), "read"); assert(fin);
  lumiUsed_ = GetLumiUsed(fin);
  cout<<"Lumi Used is "<<lumiUsed_<<" inv pb\n";

  SampleNames["data"]="Data";
  SampleNames["WJetsToLNu"]="W+Jets";
  SampleNames["WWTo2L2Nu"]="WW\\rightarrow2l2\\nu";
  SampleNames["TTJets"]="t\\bar{t}";
  SampleNames["Summer11_TTJets"]="t\\bar{t}";
  SampleNames["ZZTo4L_TuneZ2"]="ZZ\\rightarrow4l";
  SampleNames["PhotonVJets"]="V\\gamma";
  SampleNames["VV"]="VV";
  SampleNames["DYJetsToLL"]="DY+Jets\\rightarrow2l";
  SampleNames["Summer11_DYJetsToLL"]="DY+Jets\\rightarrow2l";
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
  SampleNames["WZJetsTo3LNu"]="WZ\\rightarrow3l\\nu";
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
  SampleNames["RSZZ_750"]="RS 750";
  SampleNames["RSZZ_1000"]="RS 1000";
  SampleNames["RSZZ_1250"]="RS 1250";
  SampleNames["RSZZ_1500"]="RS 1500";
  SampleNames["RSZZ_1750"]="RS 1750";
  SampleNames["RSZZ_2000"]="RS 2000";
  SampleNames["WPrimeWZ_1000"]="W'1000";
  SampleNames["Summer11_WprimeToWZTo2Q2L_M-500"]="W' 500";
  SampleNames["Summer11_WprimeToWZTo2Q2L_M-1000"]="W' 1000";

  /////Data Samples
  
  if(inName.find("WprimeWZ") != string::npos || inName.find("EWKWZ") != string::npos){
  Data.push_back(Sample("data"));
  }else if(inName.find("HadVZ") != string::npos){
    vector<string> vData;

    vData.push_back("data_DoubleMu-Run2011A-May10ReReco-v1");
    vData.push_back("data_SingleMu-Run2011A-May10ReReco-v1");
    vData.push_back("data_DoubleElectron-Run2011A-May10ReReco-v1");
    vData.push_back("data_SingleElectron-Run2011A-May10ReReco-v1");

    vData.push_back("data_DoubleMu-Run2011A-PromptReco-v4");
    vData.push_back("data_SingleMu-Run2011A-PromptReco-v4");
    vData.push_back("data_DoubleElectron-Run2011A-PromptReco-v4");
    vData.push_back("data_SingleElectron-Run2011A-PromptReco-v4");

    vData.push_back("data_DoubleMu-Run2011A-05Aug2011-v1");
    vData.push_back("data_SingleMu-Run2011A-05Aug2011-v1");
    vData.push_back("data_DoubleElectron-Run2011A-05Aug2011-v1");
    vData.push_back("data_SingleElectron-Run2011A-05Aug2011-v1");

    vData.push_back("data_DoubleMu-Run2011A-PromptReco-v6");
    vData.push_back("data_SingleMu-Run2011A-PromptReco-v6");
    vData.push_back("data_DoubleElectron-Run2011A-PromptReco-v6");
    vData.push_back("data_SingleElectron-Run2011A-PromptReco-v6");
    
    vData.push_back("data_DoubleMu-Run2011B-PromptReco-v1");
    vData.push_back("data_DoubleElectron-Run2011B-PromptReco-v1");
    
    Data.push_back(Sample("data", vData, 1, 0, 0));
      
} 
  CheckSamples(fin,Data);
    
  /////Background Samples

  if(inName.find("WprimeWZ") != string::npos || inName.find("EWKWZ") != string::npos){
    Bkg.push_back(Sample("WJetsToLNu", kOrange+6, 1, kOrange+10));

    //Bkg.push_back(Sample("ZZ", kOrange+3, 1, kOrange+10));
    //Bkg.push_back(Sample("GVJets", kOrange+3, 1, kOrange+10));
    //Bkg.push_back(Sample("WWTo2L2Nu", kOrange+3, 1, kOrange+10));
    vector<string> VV;
    VV.push_back("ZZ");
    VV.push_back("GVJets");
    VV.push_back("WWTo2L2Nu");
    Bkg.push_back(Sample("VV", VV, kOrange+3, 1, kOrange+3));
    
    Bkg.push_back(Sample("TTJets"  , kViolet+4, 1, kViolet+2));
    
    vector<string> ZJets; 
    ZJets.push_back("DYJetsToLL");
    Bkg.push_back(Sample("ZJets", ZJets, kOrange+3, 1, kOrange+7));
 
    Bkg.push_back(Sample("WZJetsTo3LNu"       , kOrange+3, 1, kOrange-2));
  }else if(inName.find("HadVZ") != string::npos){
    vector<string> VV;
    VV.push_back("Summer11_ZZ");
    VV.push_back("Summer11_VGamma");
    VV.push_back("Summer11_WW");
    VV.push_back("Summer11_WZ");
    Bkg.push_back(Sample("VV", VV, kGreen-7, 1, kGreen-7));
    
    Bkg.push_back(Sample("Summer11_TTJets"  , kBlue-7, 1, kBlue-7));
    
    vector<string> ZJets; 
    ZJets.push_back("Summer11_DYJetsToLL_PtZ100");
    Bkg.push_back(Sample("ZJets", ZJets, kRed-7, 1, kRed-7));
    
  }
  CheckSamples(fin,Bkg);

  /////Signal Samples

  if(inName.find("WprimeWZ") != string::npos){
    //Sig.push_back(Sample("WprimeToWZTo3LNu_M-300", 1, 1, kGreen));
    //Sig.push_back(Sample("WprimeToWZTo3LNu_M-400", 1, 1, 0));
    //Sig.push_back(Sample("WprimeToWZTo3LNu_M-500", 1, 1, 0));
    Sig.push_back(Sample("WprimeToWZTo3LNu_M-600", 1, 1, 0));
    //Sig.push_back(Sample("WprimeToWZTo3LNu_M-700", 1, 1, 0));
    //Sig.push_back(Sample("WprimeToWZTo3LNu_M-800", 1, 1, 0));
    //Sig.push_back(Sample("WprimeToWZTo3LNu_M-900", 1, 1, 0));

    //Sig.push_back(Sample("TC225",     2, 1, 0));
    //Sig.push_back(Sample("TC_WZ_300",     1, 1, kBlue));
    //Sig.push_back(Sample("TC_WZ_400",     1, 1, kBlue));
    //Sig.push_back(Sample("TC_WZ_500",     1, 1, kRed));
  }else if(inName.find("HadVZ") != string::npos){
    vector<string> RS750;
    RS750.push_back("Summer11_RSZZeejj_750");
    RS750.push_back("Summer11_RSZZmmjj_750");
    //Bkg.push_back(Sample("RSZZ_750", RS750, kMagenta, 1, 0));
    Sig.push_back(Sample("RSZZ_750", RS750, kMagenta+3, 1, 0));                                                                                            
    vector<string> RS1000;
    RS1000.push_back("Summer11_RSZZeejj_1000");
    RS1000.push_back("Summer11_RSZZmmjj_1000");
    //Bkg.push_back(Sample("RSZZ_1000", RS1000, kViolet, 1, 0));
    //    Sig.push_back(Sample("RSZZ_1000", RS1000, kViolet+4, 2, 0));

    vector<string> WPrime1000;
    WPrime1000.push_back("Summer11_WPrimeZZlljj_1000");
    Sig.push_back(Sample("WPrimeWZ_1000", WPrime1000, kViolet+4, 2, 0));
 
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
  }else{
    cerr<<" Don't know what you're trying to plot with input file. "
        <<inName<<endl;
    abort();
  }
  CheckSamples(fin,Sig);

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
  if(inName.find("WprimeWZ") != string::npos){
    variable.push_back("hWZMass");
  }else if(inName.find("EWKWZ") != string::npos){
      variable.push_back("hZMass");
      variable.push_back("hWTransMass");
      variable.push_back("hMET");
  }else if(inName.find("HadVZ") != string::npos){
    variable.push_back("hVZMass");
    variable.push_back("hVZeeMass");
    variable.push_back("hVZmmMass");
  }

  //These variables will be plotted after each cut unless you say "show"
  //They are cross checks but not critical to be shown
  if(opt.find("show") == string::npos){
    if(inName.find("WprimeWZ") != string::npos){
      variable.push_back("hHt");
      variable.push_back("hMET");
      variable.push_back("hWTransMass");
      variable.push_back("hWpt");
      variable.push_back("hZMass");
      variable.push_back("hZpt");
      variable.push_back("hZeept");
      variable.push_back("hZmmpt");
      
      variable.push_back("hWZ3e0muMass");
      variable.push_back("hWZ2e1muMass");
      variable.push_back("hWZ1e2muMass");
      variable.push_back("hWZ0e3muMass");
      
      variable.push_back("hWZTransMass");
      variable.push_back("hWZpt");
    
      variable.push_back("hQ");       

      variable.push_back("hEvtType");
      variable.push_back("hEvtTypeP");
      variable.push_back("hEvtTypeM");

      variable.push_back("hZMass");
      variable.push_back("hZeeMass");
      variable.push_back("hZmmMass");

      variable.push_back("hWTransMass");
      variable.push_back("hWenuTransMass");
      variable.push_back("hWmnuTransMass");

      variable.push_back("hMET");
      variable.push_back("hMETSig");
      variable.push_back("hMETee");
      variable.push_back("hMETmm");

      variable.push_back("hWenuCombRelIso");
      variable.push_back("hWmnuCombRelIso");

      variable.push_back("hLeadPt");
      variable.push_back("hLeadElecPt");
      variable.push_back("hLeadMuonPt");

      variable.push_back("hNJets");

      variable.push_back("hNVtxs");

      variable.push_back("hNLLeps");

      variable.push_back("hLeadPt");
      variable.push_back("hLeadPtZee");
      variable.push_back("hLeadPtZmm");

      variable.push_back("hTriLepMass");
      variable.push_back("hL1FastJet");

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

      variable.push_back("hZMass");
      variable.push_back("hZpt");

      variable.push_back("hVMass");
      variable.push_back("hVpt");

      variable.push_back("hNLJets");
      variable.push_back("hNLLeps");
    }
  }

  TCanvas c1; 
  c1.Print(Form("%s[", outName.c_str()), "pdf"); //Open .pdf

  //Plot the # of Evts plot
  string efftitle[] = {"hNumEvts"};
  if(opt.find("show") == string::npos){
    //This is the # of evts after each cut
    DrawandSave(fin,outName,efftitle[0],"Title: "+efftitle[0], 1, 1, 0);
    //This is the # of evts after each cut (stacked)
    DrawandSave(fin,outName,efftitle[0],"Title: "+efftitle[0], 1, 0, 0);
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
  if(inName.find("WprimeWZ") != string::npos){
    if(opt.find("show") == string::npos) {
      DrawandSave(fin,outName,"hZeeMass_Ht","Title: Z Mass After Ht Zee",0,0,0);
      DrawandSave(fin,outName,"hZmmMass_Ht","Title: Z Mass After Ht Zmm",0,0,0);

      DrawandSave(fin,outName,"hZMass_ValidZ","Title: Z Mass After Valid Z",1,0,0);
      DrawandSave(fin,outName,"hZeeMass_ValidZ","Title: Z Mass After Valid Zee",1,0,0);
      DrawandSave(fin,outName,"hZmmMass_ValidZ","Title: Z Mass After Valid Zmm",1,0,0);

      DrawandSave(fin,outName,"hZMass_ValidZ","Title: Z Mass After Valid Z",0,0,0);
      DrawandSave(fin,outName,"hZeeMass_ValidZ","Title: Z Mass After Valid Zee",0,0,0);
      DrawandSave(fin,outName,"hZmmMass_ValidZ","Title: Z Mass After Valid Zmm",0,0,0);

      DrawandSave(fin,outName,"hWTransMass_ValidW","Title: W TMass After Valid W",1,0,0);
      DrawandSave(fin,outName,"hWenuTransMass_ValidW","Title: W TMass After Valid Wen",1,0,0);
      DrawandSave(fin,outName,"hWmnuTransMass_ValidW","Title: W TMass After Valid Wmu",1,0,0);

//Adding bunch of MET plots
      DrawandSave(fin,outName,"hMET_ValidW","Title: MET After Valid W",1,0,0);
      DrawandSave(fin,outName,"hMETee_ValidW","Title: MET After Valid W (zee)",1,0,0);
      DrawandSave(fin,outName,"hMETmm_ValidW","Title: MET After Valid W (Zmm)",1,0,0);

      DrawandSave(fin,outName,"hMET_ValidW","Title: MET After Valid W",1,0,1);
      DrawandSave(fin,outName,"hMETee_ValidW","Title: MET After Valid W (zee)",1,0,1);
      DrawandSave(fin,outName,"hMETmm_ValidW","Title: MET After Valid W (Zmm)",1,0,1);

      DrawandSave(fin,outName,"hMET_ValidW","Title: MET After Valid W",0,0,0);
      DrawandSave(fin,outName,"hMETee_ValidW","Title: MET After Valid W (zee)",0,0,0);
      DrawandSave(fin,outName,"hMETmm_ValidW","Title: MET After Valid W (Zmm)",0,0,0);

      DrawandSave(fin,outName,"hMET_ValidW","Title: MET After Valid W",0,0,1);
      DrawandSave(fin,outName,"hMETee_ValidW","Title: MET After Valid W (zee)",0,0,1);
      DrawandSave(fin,outName,"hMETmm_ValidW","Title: MET After Valid W (Zmm)",0,0,1);
//End MET

      DrawandSave(fin,outName,"hHt_ValidWZCand","Title: Ht before Ht Cut",1,0,0);

      DrawandSave(fin,outName,"hWZMass_ValidWZCand","Title: WZ Mass before Ht Cut",1,0,0);
      DrawandSave(fin,outName,"hWZMass_ValidWZCand","Title: WZ Mass before Ht Cut",1,0,1);
      DrawandSave(fin,outName,"hWZMass_Ht","Title: WZ Mass After Ht Cut ",1,0,0);
      DrawandSave(fin,outName,"hWZMass_Ht","Title: WZ Mass After Ht Cut (Cumlative) ",1,0,1);

      ////////////////
      DrawandSave(fin,outName,"hNLLeps_ValidWZCand","Title: NLepAfter Valid WZ",1,0,0);
      DrawandSave(fin,outName,"hNJets_ValidWZCand","Title: NJet After Valid WZ",1,0,0);
      DrawandSave(fin,outName,"hZeept_ValidWZCand","Title: ZeePt After Valid WZ",1,0,0);
      DrawandSave(fin,outName,"hZmmpt_ValidWZCand","Title: ZmmPt After Valid WZ",1,0,0);
      DrawandSave(fin,outName,"hNVtxs_ValidWZCand","Title: NVtx After Valid WZ",1,0,0);
      DrawandSave(fin,outName,"hLeadPt_ValidWZCand","Title: Lead Pt After Valid WZ",1,0,0);

      DrawandSave(fin,outName,"hHt_ValidWZCand","Title: Ht After Valid WZ",1,0,0);
      DrawandSave(fin,outName,"hZpt_ValidWZCand","Title: Zpt After Valid WZ",1,0,0);
      DrawandSave(fin,outName,"hWpt_ValidWZCand","Title: Wpt After Valid WZ",1,0,0);

      DrawandSave(fin,outName,"hWZ3e0muMass_ValidWZCand","Title: WZ 3e0mu Mass before Ht Cut",1,0,0);
      DrawandSave(fin,outName,"hWZ2e1muMass_ValidWZCand","Title: WZ 2e1mu Mass before Ht Cut",1,0,0);
      DrawandSave(fin,outName,"hWZ1e2muMass_ValidWZCand","Title: WZ 1e2mu Mass before Ht Cut",1,0,0);
      DrawandSave(fin,outName,"hWZ0e3muMass_ValidWZCand","Title: WZ 0e3mu Mass before Ht Cut",1,0,0);
     
     }
  }else if(inName.find("EWKWZ") != string::npos){
  }else if(inName.find("HadVZ") != string::npos){
    if(opt.find("show") == string::npos) {
      DrawandSave(fin,outName,"h_bestmass","Title: Best Mass",1,0,0);
      DrawandSave(fin,outName,"hVMass_ValidV","Title: Leading Jet Mass",1,0,0);

      DrawandSave(fin,outName,"h_deltaR_elec1elec2","Title: DeltaR ee",1,0,0);
      DrawandSave(fin,outName,"h_deltaR_muon1muon2","Title: DeltaR mm",1,0,0);

      DrawandSave(fin,outName,"h_deltaR_HadVelec1","Title: DeltaR Ve1",1,0,0);
      DrawandSave(fin,outName,"h_deltaR_HadVelec2","Title: DeltaR Ve2",1,0,0);
      DrawandSave(fin,outName,"h_deltaR_HadVmuon1","Title: DeltaR Vm1",1,0,0);
      DrawandSave(fin,outName,"h_deltaR_HadVmuon2","Title: DeltaR Vm2",1,0,0);

      DrawandSave(fin,outName,"h_deltaR_jet1Z_R1","Title: DeltaR Zj1",1,0,0);
      DrawandSave(fin,outName,"h_deltaR_jet1Z_R2","Title: DeltaR Zj2",1,0,0);
      
    }
    
  }

  c1.Print(Form("%s]", outName.c_str()), "pdf"); //Close .pdf

}

void
DrawandSave(TFile* fin, string pdfName, string title, string bookmark, bool logy, bool eff, bool cum, TLine* line){
  if(!GetHistograms(fin, title, eff, cum)) return;
  string filename = "plots/" + title;
  if(!logy) filename += "_Linear"; 
  if(cum)  filename += "_Cumlative";
  filename += ".pdf";
  Draw(filename, pdfName, bookmark, logy, eff, line);
}

bool
GetHistograms(TFile* fin, string title, bool eff, bool cum){
  if(debug_) printf("  GetHisto %s\n", title.c_str());
  string hist_name;

  //How many bins should I combine
  int rebin = (title.find("WZMass") != string::npos || 
               title.find("WZ3e0muMass") != string::npos ||
               title.find("WZ2e1muMass") != string::npos ||
               title.find("WZ1e2muMass") != string::npos ||
               title.find("WZ0e3muMass") != string::npos) ? 5 : 0;
  if(title.find("VZMass") != string::npos ||
     title.find("VZeeMass") != string::npos ||
     title.find("VZmmMass") != string::npos) rebin = 2;

  //  if(title.find("VMass") != string::npos) rebin = 3;

  if(title.find("hMET_") != string::npos) rebin = 1;
  if(title.find("hHt_") != string::npos) rebin = 1;
  if(title.find("hWZTransMass_") != string::npos) rebin = 2;
  if(title.find("hWZpt_") != string::npos) rebin = 2;
  
  bool validHist = false;//If none of the histograms are filled ,don't draw it
  for(unsigned int i=0; i<samples_.size(); ++i){//Loop over Data, Bkg, Sig
    for(unsigned int j=0; j<samples_[i]->size(); ++j){
      if(debug_) cout<<"i: "<<i<<" j:"<<j<<endl;
      Sample& curSample = samples_[i]->at(j);//Let's make this easy

      curSample.hist = get_sum_of_hists(fin, curSample.names, title, rebin);
      if(!validHist && curSample.hist->Integral() > 0)
        validHist = true;
      curSample.hist->SetLineStyle(curSample.style);
      curSample.hist->SetLineColor(curSample.line); 

      if (title.find("VMass") != string::npos){
	int min = curSample.hist->GetXaxis()->FindBin(70.);
	int max = curSample.hist->GetXaxis()->FindBin(120.);
	curSample.hist->GetXaxis()->SetRangeUser(70,120);
      }

      if(!eff){
        curSample.hist->SetFillColor(curSample.fill);
        if(cum){//Do you want a cumlative distrubution
          curSample.hist = FillCum(curSample.hist);
          string newtitle = "Cumlative ";
          newtitle += curSample.hist->GetYaxis()->GetTitle();
          curSample.hist->SetYTitle(newtitle.c_str());
        }
      }else{//Eff Histos
        if(i!=1) curSample.hist->SetLineColor(curSample.fill); //If not signal
        
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
  return validHist;
}

void
Draw(string filename, string pdfName, string bookmark, bool logy, bool eff, TLine* line){
  if(debug_) printf("  Draw %s\n", filename.c_str());

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
    
    //Find Maximum
    Double_t max = sBkg->GetMaximum();
    sBkg->Draw("HIST");
    for(unsigned int i=0; i<Sig.size(); i++){
      max = TMath::Max(max, sSigs[i]->GetMaximum());
      sSigs[i]->Draw("HIST SAME");
    }

    if(drawLowSig_){//This is if you want to see the shape of a low xsec signal
      for(unsigned int i=0; i<Sig.size(); i++){
        Sig[i].hist->Draw("HIST SAME");
      }
    }

    if(hData){//Draw Data if it exists
      max = TMath::Max(max, hData->GetMaximum());
      hData->Draw("E1 SAME");
    }

    if(logy){
      sBkg->SetMaximum(50*max);
      sBkg->SetMinimum(0.1);
    }else{
      sBkg->SetMaximum(1.5*max);
      sBkg->SetMinimum(0.);
    }

    sBkg->GetXaxis()->SetTitleFont(132);
    sBkg->GetYaxis()->SetTitleFont(132);
    sBkg->GetXaxis()->SetTitleSize(0.06);
    sBkg->GetYaxis()->SetTitleSize(0.06);

  }else{//Eff Histos
    THStack* hs = new THStack("hs",title.c_str());
    for(unsigned int i=0; i<Bkg.size(); ++i){
      Bkg[i].hist->SetLineWidth(3);
      hs->Add(Bkg[i].hist);
    }
    for(unsigned int i=0; i<Sig.size(); ++i){
      Sig[i].hist->SetLineWidth(3);
      hs->Add(Sig[i].hist);
    }

    hs->Add(hData, "E1");
    
    hs->Draw("nostack");
    hs->GetXaxis()->SetTitleFont(132);
    hs->GetYaxis()->SetTitleFont(132);
    hs->GetXaxis()->SetTitleSize(0.06);
    hs->GetYaxis()->SetTitleSize(0.06);
  }
  if(debug_) cout<<"Title: "<<title<<endl;

  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextSize(0.045);
  latexLabel.DrawLatex(0.18, 0.96, "#font[132]{CMS Preliminary 2011}");
  latexLabel.DrawLatex(0.5, 0.955, "#font[132]{#sqrt{s} = 7 TeV}");
  latexLabel.DrawLatex(0.73, 0.87, Form("#font[132]{#intL dt = %.2f fb^{-1}}",lumiUsed_/1000.));

  if(debug_) cout<<"Creating Legend\n";
  TLegend *legend = new TLegend(0.62,0.69,0.92,0.85,"");
  
  if(Data.size()) legend->AddEntry(hData, "Data", "PLE");
  for(unsigned int i=0; i<Bkg.size(); ++i){
    legend->AddEntry(Bkg[i].hist,SampleNames[Bkg[i].name].c_str(), "F");
  }
  for(unsigned int i=0; i<Sig.size(); ++i){
    legend->AddEntry(Sig[i].hist,SampleNames[Sig[i].name].c_str(), "FL");
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
  if(debug_) cout<<"Checking sample "<<endl;
  for(size_t i=0; i<sample.size(); ++i){
    for(size_t j=0; j<sample[i].names.size(); ++j){
      if(debug_) cout<<"Checking key "<<sample[i].names[j]<<endl;
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

