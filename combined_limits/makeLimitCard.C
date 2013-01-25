#include "../root_macros/common.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include <fstream>

TGraphErrors* MakeSigEffGraph   (TFile* fIn, const string & moreCuts);
TGraphErrors* MakeSigEffErrGraph(TFile* fIn, const string & moreCuts);
void PrintSysErrors(const string & name, const string & type, TGraph* gSig, TGraph* gBkg, const float mass, const int nbkg, ofstream & fout);
string bin(int ch){
  if(ch == 0) return "eee";
  if(ch == 1) return "eem";
  if(ch == 2) return "mme";
  if(ch == 3) return "mmm";
  else abort();
  return "";
}
const float nGen = 6688;//60192 * 4 / 9 / nch (4/9 to remove tau events, 4=number of channels
const int nch=4;

void
makeLimitCard(string inName, int mass){
  TFile *fIn = TFile::Open(inName.c_str(), "read"); assert(fIn);

  float lumi = GetLumiUsed(fIn);

  TGraph* gxsec = new TGraph("../Limits/xSec_WZ.dat", "%lg %lg");

  TGraph* gSigMETRes   = new TGraph("../Systematics/SysSigMETRes.dat"   , "%lg %lg");
  TGraph* gSigMETScale = new TGraph("../Systematics/SysSigMETScale.dat" , "%lg %lg");
  TGraph* gSigPU       = new TGraph("../Systematics/SysSigPU.dat"       , "%lg %lg");
  TGraph* gSigMuPtRes  = new TGraph("../Systematics/SysSigMuPtRes.dat"  , "%lg %lg");
  TGraph* gSigMuPtScale= new TGraph("../Systematics/SysSigMuPtScale.dat", "%lg %lg");
  TGraph* gSigElEnScale= new TGraph("../Systematics/SysSigElEnScale.dat", "%lg %lg");
  TGraph* gSigPDF      = new TGraph("../Systematics/SysSigPDF.dat"      , "%lg %lg");

  TGraph* gBkgMETRes   = new TGraph("../Systematics/SysBkgMETRes.dat"  , "%lg %lg");
  TGraph* gBkgMETScale = new TGraph("../Systematics/SysBkgMETScale.dat", "%lg %lg");
  TGraph* gBkgPU       = new TGraph("../Systematics/SysBkgPU.dat"      , "%lg %lg");
  TGraph* gBkgMuPtRes  = new TGraph("../Systematics/SysBkgMuPtRes.dat" , "%lg %lg");
  TGraph* gBkgMuPtScale= new TGraph("../Systematics/SysBkgMuPtScale.dat", "%lg %lg");
  TGraph* gBkgElEnScale= new TGraph("../Systematics/SysBkgElEnScale.dat", "%lg %lg");
  TGraph* gBkgPDF      = new TGraph("../Systematics/SysBkgPDF.dat"     , "%lg %lg");

  //Make Sig Eff Graph
  TGraphErrors* geff   [4];
  TGraphErrors* gefferr[4];
  for(int ch=0; ch<nch; ch++){
    geff[ch]    = MakeSigEffGraph   (fIn, Form("EvtType == %i", ch));
    gefferr[ch] = MakeSigEffErrGraph(fIn, Form("EvtType == %i", ch));
  }

  map<string, Value> yields[4];

  vector<string> samples;
  samples.push_back("WZJetsTo3LNu");
  samples.push_back("DYJetsToLL");
  samples.push_back("TTJets");
  samples.push_back("ZZ");
  samples.push_back("GVJets");

  string outName = Form("card_Core_WprimeWZ_M%i.txt", mass);
  string outfile(outName.c_str());
  ofstream fout(outfile.c_str());
  if(!fout) { 
    cout << "Cannot open file " << outfile << endl; 
    abort();
  } 
  //float winWidth = WindowWidth(mass);
  //float minMass = mass - winWidth/2.;
  //float maxMass = mass + winWidth/2.;
  //float minLt   = LtCut(mass);
    
  string analysisCuts = AnalysisCuts(mass);//Form("WZMass > %.0f && WZMass < %.0f && Lt > %.0f", minMass, maxMass, minLt);
  fout<<"# Cuts are "<<analysisCuts<<" for mass="<<mass<<endl;
    
  for(int ch=0; ch<nch; ch++){
    string chCuts = Form("EvtType == %i", ch);
    string cuts = "(" + analysisCuts + " && " + chCuts + ") * weight";

    yields[ch]["data"] = GetNEvtsAndError(fIn, "data", "tEvts_MET", cuts);
    //float nDataEvt = yields[ch]["data"].val;//GetNEvts(fIn, "data", "tEvts_MET", cuts);
    //cout<<"Done with data="<<nDataEvt<<endl;
    //cout<<nDataEvt<<"  ";
      
    //nsig
    float sigEff    = geff[ch]->Eval(mass);
    float sigEffErr = gefferr[ch]->Eval(mass);
    float sigxsec = gxsec->Eval(mass);
    float nSigEvt    = lumi * sigxsec * sigEff; 
    float nSigEvtErr = 1. + (sigEff > 0. ? sigEffErr / sigEff : 0.); 

    fout<<"# ch="<<ch
        <<"  lumi="<<lumi
        <<"  sigxsec="<<sigxsec
        <<"  sigEff="<<sigEff
        <<"  sigEffErr="<<sigEffErr
        <<"  sigxsec="<<sigxsec
        <<"  nSigEvt="<<nSigEvt
        <<"  nSigEvtErr="<<nSigEvtErr
        <<endl;
    //cout<<"Done with signal="<<nSigEvt<<endl;
    yields[ch]["sig"] = Value(nSigEvt, nSigEvtErr);
    //printf("%.2f  ", nSigEvt);

    for(int isample=0; isample<(int)samples.size(); isample++){
      Value v = GetNEvtsAndError(fIn, samples[isample], "tEvts_MET", cuts);
      yields[ch][samples[isample]] = Value(v.val, 1. + (v.val > 0. ? v.err / v.val : 0.));
      //float nevts = yields[ch][samples[isample]].val;//GetNEvts(fIn, samples[isample], "tEvts_MET", cuts);
      //cout<<"Done with "<<samples[isample]<<"="<<nevts<<endl;
      //cout<<nevts<<"  ";
      //printf("%.2f  ", nevts);
    }

  }//Ch loop

  ///////////////////////
  /////////Print Out Card
  ///////////////////////
  fout<<"bin           eee   eem   mme   mmm"<<endl;
  fout<<"observation   ";
  for(int ch=0; ch<nch; ch++){
    fout<<yields[ch]["data"].val<<"  ";
  }
  fout<<endl;
  fout<<"------------"<<endl;

  //bin          eee   eee   eee   eee   eee   eee   eem   eem   eem   eem   eem   eem   emm   emm   emm   emm   emm   emm   mmm   mmm   mmm   mmm   mmm  mmm
  fout<<"bin          "; 
  for(int ch=0; ch<nch; ch++){
    for(int isample=0; isample<(int)samples.size()+1; isample++){
      fout<<bin(ch)<<"  ";
    }
  }
  fout<<endl;

  //process      sig   wz    zjet  ttbar zz    zg    sig   wz    zjet  ttbar zz    zg    sig   wz    zjet  ttbar zz    zg    sig   wz    zjet  ttbar zz   zg
  fout<<"process      ";
  for(int ch=0; ch<nch; ch++){
    fout<<"sig   ";
    for(int isample=0; isample<(int)samples.size(); isample++){
      fout<<samples[isample].substr(0,3)<<"  ";
    }
  }
  fout<<endl;
    
  //process       0     1     2     3     4     5     0     1     2     3     4     5     0     1     2     3     4     5     0     1     2     3     4     5
  fout<<"process      ";
  for(int ch=0; ch<nch; ch++){
    for(int isample=0; isample<(int)samples.size()+1; isample++){
      fout<<isample<<"    ";
    }
  }
  fout<<endl;

  //rate         20.   90    1     1     3     1     20.   90    1     1     3     1     20.   90    1     1     3     1     20.   90    1     1     3     1
  fout<<"rate         ";
  for(int ch=0; ch<nch; ch++){
    fout<<yields[ch]["sig"].val<<"  ";
    for(int isample=0; isample<(int)samples.size(); isample++){
      fout<<yields[ch][samples[isample]].val<<"  ";
    }
  }
  fout<<endl;
  fout<<"------------"<<endl;

  //Print Stat Errors
  fout<<"MCStat   lnN ";
  for(int ch=0; ch<nch; ch++){
    fout<<yields[ch]["sig"].err<<"  ";
    for(int isample=0; isample<(int)samples.size(); isample++){
      fout<<yields[ch][samples[isample]].err<<"  ";
    }
  }
  fout<<endl;

  //Print Sys Errors
  PrintSysErrors("METRes"  , "lnN", gSigMETRes  , gBkgMETRes  , mass, samples.size(),fout);
  PrintSysErrors("METScale", "lnN", gSigMETScale, gBkgMETScale, mass, samples.size(),fout);
  PrintSysErrors("PU"      , "lnN", gSigPU      , gBkgPU      , mass, samples.size(),fout);
  PrintSysErrors("MuPtRes" , "lnN", gSigMuPtRes , gBkgMuPtRes , mass, samples.size(),fout);
  PrintSysErrors("MuPtScale", "lnN", gSigMuPtScale, gBkgMuPtScale, mass, samples.size(),fout);
  PrintSysErrors("ElEnScale", "lnN", gSigElEnScale, gBkgElEnScale, mass, samples.size(),fout);
  PrintSysErrors("PDF"     , "lnN", gSigPDF     , gBkgPDF     , mass, samples.size(),fout);

  fout.close(); 
}

void
PrintSysErrors(const string & name, const string & type, TGraph* gSig, TGraph* gBkg, const float mass, const int nbkg, ofstream & fout){
  fout<<"#";//Cory: take this out when its ready
  fout<<name<<"   "<<type<<" ";
  for(int ch=0; ch<nch; ch++){
    fout<<1. + gSig->Eval(mass)<<"  ";
    for(int isample=0; isample<nbkg; isample++){
      fout<<1. + gBkg->Eval(mass)<<"  ";
    }
  }
  fout<<endl;
}

TGraphErrors*
MakeSigEffGraph(TFile* fIn, const string & moreCuts){
  TGraphErrors* g = new TGraphErrors(19);//??
  for(int mass=200, i=0; mass<=2000; mass+=100, ++i){
    float winWidth = WindowWidth(mass);
    float minMass = mass - winWidth/2.;
    float maxMass = mass + winWidth/2.;
    float minLt   = LtCut(mass);
    
    string analysisCuts = Form("WZMass > %.0f && WZMass < %.0f && Lt > %.0f", minMass, maxMass, minLt);
    string cuts = "(" + analysisCuts + " && " + moreCuts + ") * weight ";
    //get tree 
    TTree* tree = getTree(fIn, Form("WprimeToWZTo3LNu_M-%i", mass), "tEvts_MET");
    float nSigEvtsUnweighted = tree->Draw("WZMass", cuts.c_str());//Cory: Check that this works.

    //float nSigEvtsWeighted = GetNEvts(tree, cuts);
    //float nGenWeighted     = 18258 * 1.5554e-04 / nch;//ONLY FOR 1700 testing  nGen * lumi * xsec / nGen;

    double     Eff = nSigEvtsUnweighted / nGen;
    double statEff = TMath::Sqrt(Eff * (1-Eff)/nGen); 

    //cout<<" cuts: "<<cuts<<" Unweighted:"<<nSigEvtsUnweighted<<" Weighted:"<<nSigEvtsWeighted<<" ratio:"<<nSigEvtsWeighted/nSigEvtsUnweighted
    //    <<" effW:"<<nSigEvtsWeighted/nGenWeighted<<" effU:"<<Eff<<endl;
            
    g->SetPoint     (  i, mass, Eff);
    g->SetPointError(  i, 0., statEff);
    
  }
  return g;
}

TGraphErrors*
MakeSigEffErrGraph(TFile* fIn, const string & moreCuts){
  TGraphErrors* g = new TGraphErrors(19);//??
  for(int mass=200, i=0; mass<=2000; mass+=100, ++i){
    float winWidth = WindowWidth(mass);
    float minMass = mass - winWidth/2.;
    float maxMass = mass + winWidth/2.;
    float minLt   = LtCut(mass);
    
    string analysisCuts = Form("WZMass > %.0f && WZMass < %.0f && Lt > %.0f", minMass, maxMass, minLt);
    string cuts = analysisCuts + " && " + moreCuts;
    //get tree 
    TTree* tree = getTree(fIn, Form("WprimeToWZTo3LNu_M-%i", mass), "tEvts_MET");
    float nSigEvtsUnweighted = tree->Draw("WZMass", cuts.c_str());//Cory: Check that this works.
            
    double     Eff = nSigEvtsUnweighted / nGen;
    double statEff = TMath::Sqrt(Eff * (1-Eff)/nGen); 

    g->SetPoint(  i, mass, statEff);
  }
  return g;
}

