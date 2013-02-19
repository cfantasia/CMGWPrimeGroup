#include "../root_macros/common.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include <fstream>

typedef pair<TGraphErrors*, TGraphErrors*> gPair;

TGraphErrors* MakeSigEffGraph   (TFile* fIn, const string & moreCuts, float LtOffset=0., float WindOffset=0.);
TGraphErrors* MakeSigEffErrGraph(TFile* fIn, const string & moreCuts, float LtOffset=0., float WindOffset=0.);
void PrintSysErrors(const string & name, const string & type, gPair* gSys, const float mass, const int nbkg, ofstream & fout);
gPair getGraphPair(string sigFile, string bkgFile, int ch);
map<string, Value>
getYields(TFile* fIn, const int & mass, const string & cuts, const vector<string> & bkgSamples,
          const TGraph* gxsec, const TGraphErrors* geff, const TGraphErrors* gefferr);
const float nGen = 6688;//60192 * 4 / 9 / nch (4/9 to remove tau events, 4=number of channels
const int nch=4;

void
makeLimitCard(string inName, string outName, int mass, float LtOffset=0., float WindOffset=0.){
  TFile *fIn = TFile::Open(inName.c_str(), "read"); assert(fIn);
  
  float lumi = GetLumiUsed(fIn);

  TGraph* gxsec = new TGraph("../Limits/xSec_WZ.dat", "%lg %lg");
  //Make Sig Eff Graph
  TGraphErrors* geff   [nch];
  TGraphErrors* gefferr[nch];
  for(int ch=0; ch<nch; ch++){
    geff[ch]    = MakeSigEffGraph   (fIn, Form("EvtType == %i", ch), LtOffset, WindOffset);
    gefferr[ch] = MakeSigEffErrGraph(fIn, Form("EvtType == %i", ch), LtOffset, WindOffset);
  }
  
  map<string, Value> yields[nch];
  vector<string> bkgSamples = BkgSamples();
  
  ofstream fout(outName.c_str());
  if(!fout) { 
    cout << "Cannot open file " << outName << endl; 
    abort();
  } 
  
  string analysisCuts = AnalysisCuts(mass, LtOffset, WindOffset);
  fout<<"# Cuts are "<<analysisCuts<<" for mass="<<mass<<endl;
  
  for(int ch=0; ch<nch; ch++){
    string chCuts = Form("EvtType == %i", ch);
    string cuts = "(" + analysisCuts + " && " + chCuts + ") * weight"; 
    yields[ch] = getYields(fIn, mass, cuts, bkgSamples, gxsec, geff[ch], gefferr[ch]);
  }//Ch loop
  
  gPair gMETRes[nch], gMETScale[nch], gPU[nch], gMuPtRes[nch], gMuPtScale[nch], gElEnScale[nch], gPDF[nch];
  for(int ch=0; ch<nch; ch++){
    gMETRes[ch]   = getGraphPair("../Systematics/SysSigMETRes.dat"   , "../Systematics/SysBkgMETRes.dat", ch);
    gMETScale[ch] = getGraphPair("../Systematics/SysSigMETScale.dat" , "../Systematics/SysBkgMETScale.dat", ch);
    gPU[ch]       = getGraphPair("../Systematics/SysSigPU.dat"       , "../Systematics/SysBkgPU.dat", ch);
    gMuPtRes[ch]  = getGraphPair("../Systematics/SysSigMuPtRes.dat"  , "../Systematics/SysBkgMuPtRes.dat", ch);
    gMuPtScale[ch]= getGraphPair("../Systematics/SysSigMuPtScale.dat", "../Systematics/SysBkgMuPtScale.dat", ch);
    gElEnScale[ch]= getGraphPair("../Systematics/SysSigElEnScale.dat", "../Systematics/SysBkgElEnScale.dat", ch);
    gPDF[ch]      = getGraphPair("../Systematics/SysSigPDF.dat"      , "../Systematics/SysBkgPDF.dat", ch);
  }

  ////////////////////////
  /////////Print Out Table
  ////////////////////////
  if(true && mass%100 == 0){
    cout<<"Channel & ";
    for(int isample=0; isample<(int)bkgSamples.size(); isample++){
      cout<<bkgSamples[isample]<<" & ";
    }
    cout<<" & Data & Signal"<<endl;

    for(int ch=0; ch<nch; ch++){
      cout<<bin(ch)<<" & ";
      for(int isample=0; isample<(int)bkgSamples.size(); isample++){
        cout<<yields[ch][bkgSamples[isample]].val<<" & ";
      }
      cout<<yields[ch]["bkg"].val<<" & ";
      cout<<yields[ch]["data"].val<<" & ";
      cout<<yields[ch]["sig"].val<<"  ";
      cout<<endl;
    }
  }

  ///////////////////////
  /////////Print Out Card
  ///////////////////////
  for(int ch=0; ch<nch; ch++){
    fout<<"# ch="<<ch
        <<"  lumi="<<lumi
        <<"  sigEff="<<geff[ch]->Eval(mass)
        <<"  sigEffErr="<<gefferr[ch]->Eval(mass)
        <<"  sigxsec="<<gxsec->Eval(mass)
        <<"  nSigEvt="<<lumi * gxsec->Eval(mass) * geff[ch]->Eval(mass)
        <<"  nSigEvtErr="<<1. + (geff[ch]->Eval(mass) > 0. ? gefferr[ch]->Eval(mass) / geff[ch]->Eval(mass) : 0.)
        <<endl;
  }

  fout<<"bin           ";
  for(int ch=0; ch<nch; ch++) fout<<bin(ch)<<"   ";
  fout<<endl;
  fout<<"observation   ";
  for(int ch=0; ch<nch; ch++){
    fout<<yields[ch]["data"].val<<"  ";
  }
  fout<<endl;
  fout<<"------------"<<endl;

  //bin          eee   eee   eee   eee   eee   eee   eem   eem   eem   eem   eem   eem   emm   emm   emm   emm   emm   emm   mmm   mmm   mmm   mmm   mmm  mmm
  fout<<"bin          "; 
  for(int ch=0; ch<nch; ch++){
    for(int isample=0; isample<(int)bkgSamples.size()+1; isample++){
      fout<<bin(ch)<<"  ";
    }
  }
  fout<<endl;

  //process      sig   wz    zjet  ttbar zz    zg    sig   wz    zjet  ttbar zz    zg    sig   wz    zjet  ttbar zz    zg    sig   wz    zjet  ttbar zz   zg
  fout<<"process      ";
  for(int ch=0; ch<nch; ch++){
    fout<<"sig   ";
    for(int isample=0; isample<(int)bkgSamples.size(); isample++){
      fout<<bkgSamples[isample].substr(0,3)<<"  ";
    }
  }
  fout<<endl;
    
  //process       0     1     2     3     4     5     0     1     2     3     4     5     0     1     2     3     4     5     0     1     2     3     4     5
  fout<<"process      ";
  for(int ch=0; ch<nch; ch++){
    for(int isample=0; isample<(int)bkgSamples.size()+1; isample++){
      fout<<isample<<"    ";
    }
  }
  fout<<endl;

  //rate         20.   90    1     1     3     1     20.   90    1     1     3     1     20.   90    1     1     3     1     20.   90    1     1     3     1
  fout<<"rate         ";
  for(int ch=0; ch<nch; ch++){
    fout<<yields[ch]["sig"].val<<"  ";
    for(int isample=0; isample<(int)bkgSamples.size(); isample++){
      fout<<yields[ch][bkgSamples[isample]].val<<"  ";
    }
  }
  fout<<endl;
  fout<<"------------"<<endl;

  //Print Stat Errors
  fout<<"MCStat   lnN ";
  for(int ch=0; ch<nch; ch++){
    fout<<yields[ch]["sig"].err<<"  ";
    for(int isample=0; isample<(int)bkgSamples.size(); isample++){
      fout<<yields[ch][bkgSamples[isample]].err<<"  ";
    }
  }
  fout<<endl;

  //Print Sys Errors
  PrintSysErrors("METRes"   , "lnN", gMETRes   , mass, bkgSamples.size(),fout);
  PrintSysErrors("METScale" , "lnN", gMETScale , mass, bkgSamples.size(),fout);
  PrintSysErrors("PU"       , "lnN", gPU       , mass, bkgSamples.size(),fout);
  PrintSysErrors("MuPtRes"  , "lnN", gMuPtRes  , mass, bkgSamples.size(),fout);
  PrintSysErrors("MuPtScale", "lnN", gMuPtScale, mass, bkgSamples.size(),fout);
  PrintSysErrors("ElEnScale", "lnN", gElEnScale, mass, bkgSamples.size(),fout);
  PrintSysErrors("PDF"      , "lnN", gPDF      , mass, bkgSamples.size(),fout);

  fout.close(); 
}

void
PrintSysErrors(const string & name, const string & type, gPair* gSys, const float mass, const int nbkg, ofstream & fout){
  fout<<name<<"   "<<type<<" ";
  for(int ch=0; ch<nch; ch++){
    fout<<1. + gSys[ch].first->Eval(mass)<<"  ";
    for(int isample=0; isample<nbkg; isample++){
      fout<<1. + gSys[ch].second->Eval(mass)<<"  ";
    }
  }
  fout<<endl;
}

TGraphErrors*
MakeSigEffGraph(TFile* fIn, const string & moreCuts, float LtOffset, float WindOffset){
  TGraphErrors* g = new TGraphErrors(19);//??
  for(int mass=200, i=0; mass<=2000; mass+=100, ++i){
    string analysisCuts = AnalysisCuts(mass, LtOffset, WindOffset);
    string cuts = "(" + analysisCuts + " && " + moreCuts + ") * weight ";
    //get tree 
    TTree* tree = getTree(fIn, Form("WprimeToWZTo3LNu_M-%i", mass), "tEvts_MET");
    float nSigEvtsUnweighted = tree->Draw("WZMass", cuts.c_str(), "goff");

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
MakeSigEffErrGraph(TFile* fIn, const string & moreCuts, float LtOffset, float WindOffset){
  TGraphErrors* g = new TGraphErrors(19);//??
  for(int mass=200, i=0; mass<=2000; mass+=100, ++i){
    string analysisCuts = AnalysisCuts(mass, LtOffset, WindOffset);
    string cuts = "(" + analysisCuts + " && " + moreCuts + ") * weight ";
    //get tree 
    TTree* tree = getTree(fIn, Form("WprimeToWZTo3LNu_M-%i", mass), "tEvts_MET");
    float nSigEvtsUnweighted = tree->Draw("WZMass", cuts.c_str(), "goff");
            
    double     Eff = nSigEvtsUnweighted / nGen;
    double statEff = TMath::Sqrt(Eff * (1-Eff)/nGen); 

    g->SetPoint(  i, mass, statEff);
  }
  return g;
}

gPair
getGraphPair(string sigFile, string bkgFile, int ch){
  string format = "%lg";//mass
  for(int i=0; i<ch; i++) format += " %*lg %*lg";//value and err
  format += " %lg %lg";
  //cout<<"for ch: "<<ch<<" format is:"<<format<<endl;
  return make_pair(new TGraphErrors(sigFile.c_str(), format.c_str()), new TGraphErrors(bkgFile.c_str(), format.c_str()));
}

map<string, Value>
getYields(TFile* fIn, const int & mass, const string & cuts, const vector<string> & bkgSamples,
          const TGraph* gxsec, const TGraphErrors* geff, const TGraphErrors* gefferr){
  map<string, Value> yields;

  //Data
  yields["data"] = GetNEvtsAndError(fIn, "data", "tEvts_MET", cuts);
  
  //Signal
  float lumi = GetLumiUsed(fIn);
  float sigEff    = geff->Eval(mass);
  float sigEffErr = gefferr->Eval(mass);
  float sigxsec = gxsec->Eval(mass);
  float nSigEvt    = lumi * sigxsec * sigEff; 
  float nSigEvtErr = 1. + (sigEff > 0. ? sigEffErr / sigEff : 0.); 
  
  yields["sig"] = Value(nSigEvt, nSigEvtErr);
  yields["sigeff"] = Value(sigEff, sigEffErr);
  
  //Background
  Value vBkg(0,0);
  for(int isample=0; isample<(int)bkgSamples.size(); isample++){
    Value v = GetNEvtsAndError(fIn, bkgSamples[isample], "tEvts_MET", cuts);
    yields[bkgSamples[isample]] = Value(v.val, 1. + (v.val > 0. ? v.err / v.val : 0.));
    vBkg += v;
  }
  yields["bkg"] = vBkg;

  return yields;
}
