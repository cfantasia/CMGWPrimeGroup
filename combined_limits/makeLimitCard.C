#include "../root_macros/common.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include <fstream>

typedef pair<TGraphErrors*, TGraphErrors*> gPair;

TGraphErrors* MakeSigEffGraph   (TFile* fIn, const string & moreCuts, float LtOffset, float WindOffset);
TGraphErrors* MakeSigEffErrGraph(TFile* fIn, const string & moreCuts, float LtOffset, float WindOffset);
void PrintSysErrors(const string & name, const string & type, gPair* gSys, const float mass, const int nbkg, ofstream & fout);
gPair getGraphPair(string sigFile, string bkgFile, int ch);
const float nGen = 6688;//60192 * 4 / 9 / nch (4/9 to remove tau events, 4=number of channels
const int nch=4;

void
makeLimitCard(string inName, string outName, int mass, float LtOffset=0., float WindOffset=0.){
  TFile *fIn = TFile::Open(inName.c_str(), "read"); assert(fIn);

  float lumi = GetLumiUsed(fIn);

  TGraph* gxsec = new TGraph("../Limits/xSec_WZ.dat", "%lg %lg");
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

  //Make Sig Eff Graph
  TGraphErrors* geff   [nch];
  TGraphErrors* gefferr[nch];
  for(int ch=0; ch<nch; ch++){
    geff[ch]    = MakeSigEffGraph   (fIn, Form("EvtType == %i", ch), LtOffset, WindOffset);
    gefferr[ch] = MakeSigEffErrGraph(fIn, Form("EvtType == %i", ch), LtOffset, WindOffset);
  }

  map<string, Value> yields[nch];

  vector<string> samples;
  samples.push_back("WZJetsTo3LNu");
  samples.push_back("DYJetsToLL");
  samples.push_back("TTJets");
  samples.push_back("ZZ");
  samples.push_back("GVJets");

  string outfile(outName.c_str());
  ofstream fout(outfile.c_str());
  if(!fout) { 
    cout << "Cannot open file " << outfile << endl; 
    abort();
  } 
    
  string analysisCuts = AnalysisCuts(mass, LtOffset, WindOffset);//Form("WZMass > %.0f && WZMass < %.0f && Lt > %.0f", minMass, maxMass, minLt);
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
  PrintSysErrors("METRes"   , "lnN", gMETRes   , mass, samples.size(),fout);
  PrintSysErrors("METScale" , "lnN", gMETScale , mass, samples.size(),fout);
  PrintSysErrors("PU"       , "lnN", gPU       , mass, samples.size(),fout);
  PrintSysErrors("MuPtRes"  , "lnN", gMuPtRes  , mass, samples.size(),fout);
  PrintSysErrors("MuPtScale", "lnN", gMuPtScale, mass, samples.size(),fout);
  PrintSysErrors("ElEnScale", "lnN", gElEnScale, mass, samples.size(),fout);
  PrintSysErrors("PDF"      , "lnN", gPDF      , mass, samples.size(),fout);

  fout.close(); 
}

void
PrintSysErrors(const string & name, const string & type, gPair* gSys, const float mass, const int nbkg, ofstream & fout){
  //fout<<"#";//Cory: take this out when its ready
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
MakeSigEffErrGraph(TFile* fIn, const string & moreCuts, float LtOffset, float WindOffset){
  TGraphErrors* g = new TGraphErrors(19);//??
  for(int mass=200, i=0; mass<=2000; mass+=100, ++i){
    string analysisCuts = AnalysisCuts(mass, LtOffset, WindOffset);
    string cuts = "(" + analysisCuts + " && " + moreCuts + ") * weight ";
    //get tree 
    TTree* tree = getTree(fIn, Form("WprimeToWZTo3LNu_M-%i", mass), "tEvts_MET");
    float nSigEvtsUnweighted = tree->Draw("WZMass", cuts.c_str());
            
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
