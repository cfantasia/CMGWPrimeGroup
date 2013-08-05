#include "../root_macros/common.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include <fstream>
#include <map>
#include <string>

using namespace std;

typedef pair<TGraphErrors*, TGraphErrors*> gPair;
const float nGen = 6688;//60192 * 4 / 9 / nch (4/9 to remove tau events, 4=number of channels
const int nch=4;

struct Object{
  double    nEvt;
  double statErr;
  double totSysErr;
  double lumiErr;
  double  totErr;
  map<string, double> sysErr;
  double relStatErr(){ return statErr/nEvt;}
  double relSysErr(const string & src){ return sysErr[src]/nEvt;}
};

vector<string> sysSrc;

TGraphErrors* MakeSigEffGraph   (TFile* fIn, const string & moreCuts, float LtOffset=0., float WindOffset=0., int nChUsed=1);
TGraphErrors* MakeSigEffErrGraph(TFile* fIn, const string & moreCuts, float LtOffset=0., float WindOffset=0., int nChUsed=1);
void PrintCardHead(ofstream & fout, int nChUsed=nch);
void PrintCardTail(const vector<string> & bkgSamples, const int mass, ofstream & fout, int startCh=0, int endCh=nch);
void PrintSysErrors(const string & name, const string & type, gPair* gSys, const float mass, const int nbkg, ofstream & fout, int startCh=0, int endCh=nch);
gPair getGraphPair(string sigFile, string bkgFile, int ch);
map<string, Value>
getYields(TFile* fIn, const int & mass, const string & cuts, const vector<string> & bkgSamples,
          const TGraph* gXSecCh, const TGraphErrors* gEffCh, const TGraphErrors* gEffChErr);

map<string, Value>
getSysErrs(const int & mass, const vector<string> & bkgSamples, map<string, Value> & yields, int ch);

map<string, Value>
getSysErrSum(map<string, Value> sysErrs[nch], const vector<string> & bkgSamples);
void
StoreSysErrors(map<string, Value> & sysErrs, gPair gSys, const float mass, const vector<string> & bkgSamples);
void
StoreSysErrorsTail(map<string, Value> & sysErrs, const vector<string> & bkgSamples, int ch);


void PrintCard(vector<map<string, Object> > & yields, const string & outName, const vector<string> & bkgSamples );


void
makeLimitCard(string inName, string outName, int mass, float LtOffset=0., float WindOffset=0.){
  TFile *fIn = TFile::Open(inName.c_str(), "read"); assert(fIn);
  
  float lumi = GetLumiUsed(fIn);

  TGraph* gXSec = new TGraph("../Limits/xSec_WZ.dat", "%lg %lg");
  TGraph* gXSecCh = new TGraph(*gXSec);
  Scale(gXSecCh, 1./nch);

  //Make Sig Eff Graph
  TGraphErrors* gEffCh   [nch];
  TGraphErrors* gEffChErr[nch];
  for(int ch=0; ch<nch; ch++){
    gEffCh[ch]    = MakeSigEffGraph   (fIn, Form("EvtType == %i", ch), LtOffset, WindOffset);
    gEffChErr[ch] = MakeSigEffErrGraph(fIn, Form("EvtType == %i", ch), LtOffset, WindOffset);
  }
  
  map<string, Value> yields[nch];
  vector< map<string, Object> > objs(4);

  vector<string> bkgSamples = BkgSamples();

  //Sys Sources
  sysSrc.push_back("METRes");
  //sysSrc.push_back("METScale");
  sysSrc.push_back("PU");
  sysSrc.push_back("MuPtRes");
  sysSrc.push_back("MuPtScale");
  sysSrc.push_back("ElEnScale");
  sysSrc.push_back("PDF");
  sysSrc.push_back("ElTrig");
  sysSrc.push_back("ElReco");
  sysSrc.push_back("ElIDIso");
  sysSrc.push_back("MuTrig");
  sysSrc.push_back("MuReco");
  sysSrc.push_back("MuIDIso");
  sysSrc.push_back("MatrixM");
  sysSrc.push_back("ZGxsec");
  sysSrc.push_back("ZZxsec");
  sysSrc.push_back("WZxsec");
  sysSrc.push_back("WZNLO");
  sysSrc.push_back("lumi");

  
  ofstream fout(outName.c_str());
  if(!fout) { 
    cout << "Cannot open file " << outName << endl; 
    abort();
  } 
  PrintCardHead(fout);
  
  string analysisCuts = AnalysisCuts(mass, LtOffset, WindOffset);
  fout<<"# Cuts are "<<analysisCuts<<" for mass="<<mass<<endl;
  
  for(int ch=0; ch<nch; ch++){
    string chCuts = Form("EvtType == %i", ch);
    string cuts = "(" + analysisCuts + " && " + chCuts + ") * weight"; 
    yields[ch] = getYields(fIn, mass, cuts, bkgSamples, gXSecCh, gEffCh[ch], gEffChErr[ch]);
  }//Ch loop
  
  gPair gMETRes[nch], gMETScale[nch], gPU[nch], gMuPtRes[nch], gMuPtScale[nch], gElEnScale[nch], gPDF[nch];
  for(int ch=0; ch<nch; ch++){
    gMETRes[ch]   = getGraphPair("../Systematics/SysSigMETRes.dat"   , "../Systematics/SysBkgMETRes.dat", ch);
    //gMETScale[ch] = getGraphPair("../Systematics/SysSigMETScale.dat" , "../Systematics/SysBkgMETScale.dat", ch);
    gPU[ch]       = getGraphPair("../Systematics/SysSigPU.dat"       , "../Systematics/SysBkgPU.dat", ch);
    gMuPtRes[ch]  = getGraphPair("../Systematics/SysSigMuPtRes.dat"  , "../Systematics/SysBkgMuPtRes.dat", ch);
    gMuPtScale[ch]= getGraphPair("../Systematics/SysSigMuPtScale.dat", "../Systematics/SysBkgMuPtScale.dat", ch);
    gElEnScale[ch]= getGraphPair("../Systematics/SysSigElEnScale.dat", "../Systematics/SysBkgElEnScale.dat", ch);
    gPDF[ch]      = getGraphPair("../Systematics/SysSigPDF.dat"      , "../Systematics/SysBkgPDF.dat", ch);
  }

  //Everything is now set up and can be printed out however we want.


  //Test new procedure
/*
  //Store Sys Errors
  StoreSysErrors(objs, "METRes   ", "lnN", gMETRes   , mass, bkgSamples.size(),fout);
  StoreSysErrors(objs, "METScale ", "lnN", gMETScale , mass, bkgSamples.size(),fout);
  StoreSysErrors(objs, "PU       ", "lnN", gPU       , mass, bkgSamples.size(),fout);
  StoreSysErrors(objs, "MuPtRes  ", "lnN", gMuPtRes  , mass, bkgSamples.size(),fout);
  StoreSysErrors(objs, "MuPtScale", "lnN", gMuPtScale, mass, bkgSamples.size(),fout);
  StoreSysErrors(objs, "ElEnScale", "lnN", gElEnScale, mass, bkgSamples.size(),fout);
  StoreSysErrors(objs, "PDF      ", "lnN", gPDF      , mass, bkgSamples.size(),fout);
  StoreSysErrorsTail(objs, sysErrs, bkgSamples, ch);

  PrintCard(objs, "aaa.out", bkgSamples);
*/

  ///////////////////////
  /////////Print Out Card
  ///////////////////////

  //First the indivual channel card
  for(int ch=0; ch<nch; ch++){
    string outChName = outName;
    outChName.replace(outChName.find("WZ_"), 3, Form("WZCh%i_", ch));
    ofstream fCh(outChName.c_str());
    if(!fCh) { 
      cout << "Cannot open file " << outChName << endl; 
      abort();
    } 
    PrintCardHead(fCh, 1);

    fCh<<"# ch="<<ch
        <<"  lumi="<<lumi
       <<"  sigEff="<<max(0., gEffCh[ch]->Eval(mass))
       <<"  sigEffErr="<<max(0., gEffChErr[ch]->Eval(mass))
       <<"  sigxsec="<<max(0.,gXSecCh->Eval(mass))
       <<"  nSigEvt="<<lumi * max(0., gXSecCh->Eval(mass)) * max(0., gEffCh[ch]->Eval(mass))
       <<"  nSigEvtErr="<<1. + (gEffCh[ch]->Eval(mass) > 0. ? max(0., gEffChErr[ch]->Eval(mass) / gEffCh[ch]->Eval(mass)) : 0.)
        <<endl;
    //Print Data yields
    fCh<<"bin           "<<bin(ch)<<"   "<<endl;
    fCh<<"observation   "<<yields[ch]["data"].val<<"  "<<endl;
    fCh<<"------------"<<endl;

    //Print channel (bin) labels
    //bin          eee   eee   eee   eee   eee   eee   eem   eem   eem   eem   eem   eem   emm   emm   emm   emm   emm   emm   mmm   mmm   mmm   mmm   mmm  mmm
    fCh<<"bin          "; 
    for(int isample=0; isample<(int)bkgSamples.size()+1; isample++){
      fCh<<bin(ch)<<"  ";
    }
    fCh<<endl;

    //Print Sample name
    //process      sig   wz    zjet  ttbar zz    zg    sig   wz    zjet  ttbar zz    zg    sig   wz    zjet  ttbar zz    zg    sig   wz    zjet  ttbar zz   zg
    fCh<<"process      ";
    fCh<<"sig   ";
    for(int isample=0; isample<(int)bkgSamples.size(); isample++){
      fCh<<bkgSamples[isample].substr(0,3)<<"  ";
    }
    fCh<<endl;
    
    //Print Bin number
    //process       0     1     2     3     4     5     0     1     2     3     4     5     0     1     2     3     4     5     0     1     2     3     4     5
    fCh<<"process      ";
    for(int isample=0; isample<(int)bkgSamples.size()+1; isample++){
      fCh<<isample<<"    ";
    }
    fCh<<endl;

    //Print Sig/Bkg yields
    //rate         20.   90    1     1     3     1     20.   90    1     1     3     1     20.   90    1     1     3     1     20.   90    1     1     3     1
    fCh<<"rate         ";
    fCh<<yields[ch]["sig"].val<<"  ";
    for(int isample=0; isample<(int)bkgSamples.size(); isample++){
      fCh<<yields[ch][bkgSamples[isample]].val<<"  ";
    }
    fCh<<endl;
    fCh<<"------------"<<endl;

    //Print Stat Errors
    fCh<<"MCStat    lnN ";
    const Value & vSig = yields[ch]["sig"];
    fCh<<Value(1. + (vSig.val > 0. ? vSig.err / vSig.val : 0.), -4)<<"  ";
    for(int isample=0; isample<(int)bkgSamples.size(); isample++){
      const Value & vBkg = yields[ch][bkgSamples[isample]];
      fCh<<Value(1. + (vBkg.val > 0. ? vBkg.err / vBkg.val : 0.), -4)<<"  ";
    }
    fCh<<endl;

    //Print Sys Errors
    PrintSysErrors("METRes   ", "lnN", gMETRes   , mass, bkgSamples.size(),fCh,ch,ch+1);
    //PrintSysErrors("METScale ", "lnN", gMETScale , mass, bkgSamples.size(),fCh,ch,ch+1);
    PrintSysErrors("PU       ", "lnN", gPU       , mass, bkgSamples.size(),fCh,ch,ch+1);
    PrintSysErrors("MuPtRes  ", "lnN", gMuPtRes  , mass, bkgSamples.size(),fCh,ch,ch+1);
    PrintSysErrors("MuPtScale", "lnN", gMuPtScale, mass, bkgSamples.size(),fCh,ch,ch+1);
    PrintSysErrors("ElEnScale", "lnN", gElEnScale, mass, bkgSamples.size(),fCh,ch,ch+1);
    PrintSysErrors("PDF      ", "lnN", gPDF      , mass, bkgSamples.size(),fCh,ch,ch+1);
    PrintCardTail(bkgSamples, mass, fCh, ch, ch+1);
    fCh.close(); 
  
  }//Ch card loop
  
  //////////////////////////
  // Next the merged card //
  //////////////////////////

  //Some debug info
  for(int ch=0; ch<nch; ch++){
    fout<<"# ch="<<ch
        <<"  lumi="<<lumi
        <<"  sigEff="<<max(0., gEffCh[ch]->Eval(mass))
        <<"  sigEffErr="<<max(0., gEffChErr[ch]->Eval(mass))
        <<"  sigxsec="<<max(0., gXSecCh->Eval(mass))
        <<"  nSigEvt="<<lumi * max(0., gXSecCh->Eval(mass)) * max(0., gEffCh[ch]->Eval(mass))
        <<"  nSigEvtErr="<<lumi * max(0., gXSecCh->Eval(mass)) * max(0., gEffChErr[ch]->Eval(mass))
        <<endl;
  }

  //Print Data yields
  fout<<"bin           ";
  for(int ch=0; ch<nch; ch++) fout<<bin(ch)<<"   ";
  fout<<endl;
  fout<<"observation   ";
  for(int ch=0; ch<nch; ch++){
    fout<<yields[ch]["data"].val<<"  ";
  }
  fout<<endl;
  fout<<"------------"<<endl;

  //Print channel (bin) labels
  //bin          eee   eee   eee   eee   eee   eee   eem   eem   eem   eem   eem   eem   emm   emm   emm   emm   emm   emm   mmm   mmm   mmm   mmm   mmm  mmm
  fout<<"bin          "; 
  for(int ch=0; ch<nch; ch++){
    for(int isample=0; isample<(int)bkgSamples.size()+1; isample++){
      fout<<bin(ch)<<"  ";
    }
  }
  fout<<endl;

  //Print Sample name
  //process      sig   wz    zjet  ttbar zz    zg    sig   wz    zjet  ttbar zz    zg    sig   wz    zjet  ttbar zz    zg    sig   wz    zjet  ttbar zz   zg
  fout<<"process      ";
  for(int ch=0; ch<nch; ch++){
    fout<<"sig   ";
    for(int isample=0; isample<(int)bkgSamples.size(); isample++){
      fout<<bkgSamples[isample].substr(0,3)<<"  ";
    }
  }
  fout<<endl;
    
  //Print Bin number
  //process       0     1     2     3     4     5     0     1     2     3     4     5     0     1     2     3     4     5     0     1     2     3     4     5
  fout<<"process      ";
  for(int ch=0; ch<nch; ch++){
    for(int isample=0; isample<(int)bkgSamples.size()+1; isample++){
      fout<<isample<<"    ";
    }
  }
  fout<<endl;

  //Print Sig/Bkg yields
  //rate         20.   90    1     1     3     1     20.   90    1     1     3     1     20.   90    1     1     3     1     20.   90    1     1     3     1
  fout<<"rate         ";
  for(int ch=0; ch<nch; ch++){
    fout<<yields[ch]["sig"].val<<"  ";
    for(int isample=0; isample<(int)bkgSamples.size(); isample++){
      fout<<yields[ch][bkgSamples[isample]].val<<"  ";//Good line
      //if(0 == bkgSamples[isample].compare("WZJetsTo3LNu")) fout<<"0.000001 ";
      //else                                    fout<<" 0. "<<"  ";//Temp check
    }
  }
  fout<<endl;
  fout<<"------------"<<endl;

  //Print Stat Errors
  fout<<"MCStat    lnN ";
  for(int ch=0; ch<nch; ch++){
    const Value & vSig = yields[ch]["sig"];
    fout<<Value(1. + (vSig.val > 0. ? vSig.err / vSig.val : 0.), -4)<<"  ";
    for(int isample=0; isample<(int)bkgSamples.size(); isample++){
      const Value & vBkg = yields[ch][bkgSamples[isample]];
      fout<<Value(1. + (vBkg.val > 0. ? vBkg.err / vBkg.val : 0.), -4)<<"  ";
    }
  }
  fout<<endl;

  //Print Sys Errors
  PrintSysErrors("METRes   ", "lnN", gMETRes   , mass, bkgSamples.size(),fout);
  //PrintSysErrors("METScale ", "lnN", gMETScale , mass, bkgSamples.size(),fout);
  PrintSysErrors("PU       ", "lnN", gPU       , mass, bkgSamples.size(),fout);
  PrintSysErrors("MuPtRes  ", "lnN", gMuPtRes  , mass, bkgSamples.size(),fout);
  PrintSysErrors("MuPtScale", "lnN", gMuPtScale, mass, bkgSamples.size(),fout);
  PrintSysErrors("ElEnScale", "lnN", gElEnScale, mass, bkgSamples.size(),fout);
  PrintSysErrors("PDF      ", "lnN", gPDF      , mass, bkgSamples.size(),fout);
  PrintCardTail(bkgSamples, mass, fout);
  fout.close(); 

  ////////////////////
  //// Summed Card ///
  ////////////////////

  string chCuts = Form("EvtType > -1");
  string cuts = "(" + analysisCuts + " && " + chCuts + ") * weight";
  
  TGraphErrors* geffSum      = MakeSigEffGraph   (fIn, "1", 0, 0, 4.);
  TGraphErrors* gefferrSum   = MakeSigEffErrGraph(fIn, "1", 0, 0, 4.);

  map<string, Value> yieldSum = getYields(fIn, mass, cuts, bkgSamples, gXSec, geffSum, gefferrSum);

  string outSumName = outName;
  outSumName.replace(outSumName.find("WZ_"), 3, "WZSum_");
  ofstream fSum(outSumName.c_str());
  if(!fSum) { 
    cout << "Cannot open file " << outSumName << endl; 
    abort();
  } 
  PrintCardHead(fSum, 1);

  //Print Data yields
  fSum<<"bin           ";
  fSum<<"sum"<<"   ";
  fSum<<endl;
  fSum<<"observation   ";
  fSum<<yieldSum["data"].val<<"  ";
  fSum<<endl;
  fSum<<"------------"<<endl;

  //Print channel (bin) labels
  fSum<<"bin          "; 
  for(int isample=0; isample<(int)bkgSamples.size()+1; isample++){
    fSum<<"sum"<<"  ";
  }
  fSum<<endl;

  //Print Sample name
  fSum<<"process      ";
  fSum<<"sig   ";
  for(int isample=0; isample<(int)bkgSamples.size(); isample++){
    fSum<<bkgSamples[isample].substr(0,3)<<"  ";
  }
  fSum<<endl;
    
  //Print Bin number
  fSum<<"process      ";
  for(int isample=0; isample<(int)bkgSamples.size()+1; isample++){
    fSum<<isample<<"    ";
  }
  fSum<<endl;

  //Print Sig/Bkg yields
  fSum<<"rate         ";
  fSum<<yieldSum["sig"].val<<"  ";
  for(int isample=0; isample<(int)bkgSamples.size(); isample++){
    fSum<<yieldSum[bkgSamples[isample]].val<<"  ";
  }
  fSum<<endl;
  fSum<<"------------"<<endl;

  //Print Stat Errors
  fSum<<"MCStat    lnN ";
  const Value & vSig = yieldSum["sig"];
  fSum<<Value(vSig.relErr(), -4)<<"  ";
  for(int isample=0; isample<(int)bkgSamples.size(); isample++){
    const Value & vBkg = yieldSum[bkgSamples[isample]];
    fSum<<Value(vBkg.relErr(), -4)<<"  ";
  }
  fSum<<endl;

  //Print Sys Errors
  //PrintSysErrors("METRes   ", "lnN", gMETRes   , mass, bkgSamples.size(),fSum);
  //PrintSysErrors("METScale ", "lnN", gMETScale , mass, bkgSamples.size(),fSum);
  //PrintSysErrors("PU       ", "lnN", gPU       , mass, bkgSamples.size(),fSum);
  //PrintSysErrors("MuPtRes  ", "lnN", gMuPtRes  , mass, bkgSamples.size(),fSum);
  //PrintSysErrors("MuPtScale", "lnN", gMuPtScale, mass, bkgSamples.size(),fSum);
  //PrintSysErrors("ElEnScale", "lnN", gElEnScale, mass, bkgSamples.size(),fSum);
  //PrintSysErrors("PDF      ", "lnN", gPDF      , mass, bkgSamples.size(),fSum);
  //PrintCardTail(bkgSamples, mass, fSum, );
  
  fSum.close(); 


}

void
PrintSysErrors(const string & name, const string & type, gPair* gSys, const float mass, const int nbkg, ofstream & fout, int startCh, int endCh){
  fout<<name<<" "<<type<<" ";
  for(int ch=startCh; ch<endCh; ch++){
    if(gSys[ch].first->Eval(mass) > 1e-5) fout<<Value(1. + gSys[ch].first->Eval(mass), -4)<<"  ";
    else                                  fout<<" -     ";
    for(int isample=0; isample<nbkg; isample++){
      if(gSys[ch].second->Eval(mass) > 1e-5) fout<<Value(1. + gSys[ch].second->Eval(mass), -4)<<"  ";
      else                                   fout<<" -     ";
    }
  }
  fout<<endl;
}

void
PrintCardHead(ofstream & fout, int nChUsed){
  fout<<"# Counting experiment with multiple channels"<<endl;
  fout<<"# Wprime with 4 channels"<<endl;
  fout<<"imax "<<nChUsed<<"  number of channels"<<endl;
  fout<<"jmax 5  number of backgrounds ('*' = automatic)"<<endl;
  fout<<"kmax 19  number of nuisance parameters (sources of systematical uncertainties)"<<endl;
  fout<<"------------"<<endl;
}

void
PrintCardTail(const vector<string> & bkgSamples, const int mass, ofstream & fout, int startCh, int endCh){
  const int nbkg = bkgSamples.size();

  fout<<"#"<<endl;
  fout<<"ElTrig   lnN ";
  for(int ch=startCh; ch<endCh; ch++){
    for(int isample=0; isample<nbkg+1; isample++){
      if(ch<=1) fout<<"1.02  ";
      else      fout<<" -    ";
    }
  }fout<<endl;
  
  fout<<"ElReco   lnN ";
  for(int ch=startCh; ch<endCh; ch++){
    for(int isample=0; isample<nbkg+1; isample++){
      if     (ch==0) fout<<"1.06  ";
      else if(ch==1) fout<<"1.04  ";
      else if(ch==2) fout<<"1.02  ";
      else      fout<<" -    ";
    }
  }fout<<endl;

  fout<<"ElIDIso  lnN ";
  for(int ch=startCh; ch<endCh; ch++){
    for(int isample=0; isample<nbkg+1; isample++){
      if     (ch==0) fout<<"1.03  ";
      else if(ch==1) fout<<"1.02  ";
      else if(ch==2) fout<<"1.01  ";
      else      fout<<" -    ";
    }
  }fout<<endl;

  fout<<"MuTrig   lnN ";
  for(int ch=startCh; ch<endCh; ch++){
    for(int isample=0; isample<nbkg+1; isample++){
      if(ch>=2) fout<<"1.05  ";
      else      fout<<" -    ";
    }
  }fout<<endl;

  fout<<"MuReco   lnN ";
  for(int ch=startCh; ch<endCh; ch++){
    for(int isample=0; isample<nbkg+1; isample++){
      if     (ch==1) fout<<"1.02  ";
      else if(ch==2) fout<<"1.04  ";
      else if(ch==3) fout<<"1.06  ";
      else      fout<<" -    ";
    }
  }fout<<endl;

  fout<<"MuIDIso  lnN ";
  for(int ch=startCh; ch<endCh; ch++){
    for(int isample=0; isample<nbkg+1; isample++){
      if     (ch==1) fout<<"1.03  ";
      else if(ch==2) fout<<"1.06  ";
      else if(ch==3) fout<<"1.09  ";
      else      fout<<" -    ";
    }
  }fout<<endl;

  fout<<"#"<<endl;

  fout<<"MatrixM  lnN ";
  for(int ch=startCh; ch<endCh; ch++){
    fout<<" -    ";//for signal
    for(int isample=0; isample<nbkg; isample++){
      if     (0 == bkgSamples[isample].compare("TTJets"))     fout<<"1.15  ";
      else if(0 == bkgSamples[isample].compare("DYJetsToLL")) fout<<"1.30  ";
      else                                               fout<<" -    ";
    }
  }fout<<endl;

  fout<<"ZGxsec   lnN ";
  for(int ch=startCh; ch<endCh; ch++){
    fout<<" -    ";//for signal
    for(int isample=0; isample<nbkg; isample++){
      //if(0 == bkgSamples[isample].compare("GVJets")) fout<<"1.054 ";//Good value
      if(0 == bkgSamples[isample].compare("GVJets")) fout<<"1.500 ";//Temp test
      else               fout<<" -    ";
    }
  }fout<<endl;

    fout<<"ZZxsec   lnN ";
    for(int ch=startCh; ch<endCh; ch++){
      fout<<" -    ";//for signal
      for(int isample=0; isample<nbkg; isample++){
        if(0 == bkgSamples[isample].compare("ZZ")) fout<<"1.30  ";
        //if(0 == bkgSamples[isample].compare("ZZ")) fout<<"1.15  ";
        else               fout<<" -    ";
      }
    }fout<<endl;

    fout<<"WZRatio  lnN ";
    for(int ch=startCh; ch<endCh; ch++){
      fout<<" -    ";//for signal
      for(int isample=0; isample<nbkg; isample++){
        if     (                mass <  600 && 0 == bkgSamples[isample].compare("WZJetsTo3LNu")) fout<<"1.050 ";
        else if(mass >=  600 && mass < 1000 && 0 == bkgSamples[isample].compare("WZJetsTo3LNu")) fout<<"1.050 ";
        else if(mass >= 1000                && 0 == bkgSamples[isample].compare("WZJetsTo3LNu")) fout<<"1.300 ";
        else               fout<<" -    ";
      }
    }fout<<endl;

    fout<<"WZQCDSc  lnN ";
    for(int ch=startCh; ch<endCh; ch++){
      fout<<" -    ";//for signal
      for(int isample=0; isample<nbkg; isample++){
        if     (                mass <  600 && 0 == bkgSamples[isample].compare("WZJetsTo3LNu")) fout<<"1.100 ";
        else if(mass >=  600 && mass < 1000 && 0 == bkgSamples[isample].compare("WZJetsTo3LNu")) fout<<"1.150 ";
        else if(mass >= 1000                && 0 == bkgSamples[isample].compare("WZJetsTo3LNu")) fout<<"1.150 ";
        else               fout<<" -    ";
      }
    }fout<<endl;

//     fout<<"WZxsec   lnN ";
//     for(int ch=startCh; ch<endCh; ch++){
//       fout<<" -    ";//for signal
//       for(int isample=0; isample<nbkg; isample++){
//         if(0 == bkgSamples[isample].compare("WZJetsTo3LNu")) fout<<"1.075 ";
//         else               fout<<" -    ";
//       }
//     }fout<<endl;

//   fout<<"WZNLO    lnN ";
//   for(int ch=startCh; ch<endCh; ch++){
//     fout<<" -    ";//for signal
//     for(int isample=0; isample<nbkg; isample++){
//       if(0 == bkgSamples[isample].compare("WZJetsTo3LNu")) fout<<"1.70  ";
//       else               fout<<" -    ";
//     }
//   }fout<<endl;

  fout<<"lumi     lnN ";
  for(int ch=startCh; ch<endCh; ch++){
    fout<<"1.044 ";//for signal
    for(int isample=0; isample<nbkg; isample++){
      fout<<"1.044 ";
    }
  }fout<<endl;
}

TGraphErrors*
MakeSigEffGraph(TFile* fIn, const string & moreCuts, float LtOffset, float WindOffset, int nChUsed){
  TGraphErrors* g = new TGraphErrors(19);//??
  for(int mass=200, i=0; mass<=2000; mass+=100, ++i){
    string analysisCuts = AnalysisCuts(mass, LtOffset, WindOffset);
    string cuts = "(" + analysisCuts + " && " + moreCuts + ") * weight ";
    //get tree 
    TTree* tree = getTree(fIn, Form("WprimeToWZTo3LNu_M-%i", mass), "tEvts_MET");
    float nSigEvtsUnweighted = tree->Draw("WZMass", cuts.c_str(), "goff");

    double     Eff = nSigEvtsUnweighted / (nChUsed*nGen);
    double statEff = TMath::Sqrt(Eff * (1-Eff)/(nChUsed*nGen)); 

    g->SetPoint     (  i, mass, Eff);
    g->SetPointError(  i, 0., statEff);
    delete tree;
  }
  return g;
}

TGraphErrors*
MakeSigEffErrGraph(TFile* fIn, const string & moreCuts, float LtOffset, float WindOffset, int nChUsed){
  TGraphErrors* g = new TGraphErrors(19);//??
  for(int mass=200, i=0; mass<=2000; mass+=100, ++i){
    string analysisCuts = AnalysisCuts(mass, LtOffset, WindOffset);
    string cuts = "(" + analysisCuts + " && " + moreCuts + ") * weight ";
    //get tree 
    TTree* tree = getTree(fIn, Form("WprimeToWZTo3LNu_M-%i", mass), "tEvts_MET");
    float nSigEvtsUnweighted = tree->Draw("WZMass", cuts.c_str(), "goff");
            
    double     Eff = nSigEvtsUnweighted / (nChUsed*nGen);
    double statEff = TMath::Sqrt(Eff * (1-Eff)/(nChUsed*nGen)); 

    g->SetPoint(  i, mass, statEff);
    delete tree;
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
          const TGraph* gXSecCh, const TGraphErrors* gEffCh, const TGraphErrors* gEffChErr){
  map<string, Value> yields;

  //Data
  yields["data"] = GetNEvtsAndError(fIn, "data", "tEvts_MET", cuts);
  
  //Signal
  float lumi = GetLumiUsed(fIn);
  float sigEff    = max(0., gEffCh->Eval(mass));
  float sigEffErr = max(0., gEffChErr->Eval(mass));
  float sigxsec = max(0., gXSecCh->Eval(mass));
  float nSigEvt    = lumi * sigxsec * sigEff; 
  float nSigEvtErr = lumi * sigxsec * sigEffErr;
  
  yields["sig"] = Value(nSigEvt, nSigEvtErr);
  yields["sigeff"] = Value(sigEff, sigEffErr);
  
  //Background
  Value vBkg(0,0);
  for(int isample=0; isample<(int)bkgSamples.size(); isample++){
    //Value v = GetNEvtsAndError(fIn, bkgSamples[isample], "tEvts_MET", cuts, true);//Apply Scale Factors
    Value v = GetNEvtsAndError(fIn, bkgSamples[isample], "tEvts_MET", cuts, false);//Don't Apply Scale Factors
    //yields[bkgSamples[isample]] = Value(v.val, 1. + (v.val > 0. ? v.err / v.val : 0.));
    yields[bkgSamples[isample]] = v;
    vBkg += v;
  }
  yields["bkg"] = vBkg;

  return yields;
}

map<string, Value>
getSysErrs(const int & mass, const vector<string> & bkgSamples, map<string, Value> & yields, int ch){
  map<string, Value> sysErrs;

  gPair gMETRes, gMETScale, gPU, gMuPtRes, gMuPtScale, gElEnScale, gPDF;
  
  gMETRes   = getGraphPair("../Systematics/SysSigMETRes.dat"   , "../Systematics/SysBkgMETRes.dat", ch);
  //gMETScale = getGraphPair("../Systematics/SysSigMETScale.dat" , "../Systematics/SysBkgMETScale.dat", ch);
  gPU       = getGraphPair("../Systematics/SysSigPU.dat"       , "../Systematics/SysBkgPU.dat", ch);
  gMuPtRes  = getGraphPair("../Systematics/SysSigMuPtRes.dat"  , "../Systematics/SysBkgMuPtRes.dat", ch);
  gMuPtScale= getGraphPair("../Systematics/SysSigMuPtScale.dat", "../Systematics/SysBkgMuPtScale.dat", ch);
  gElEnScale= getGraphPair("../Systematics/SysSigElEnScale.dat", "../Systematics/SysBkgElEnScale.dat", ch);
  gPDF      = getGraphPair("../Systematics/SysSigPDF.dat"      , "../Systematics/SysBkgPDF.dat", ch);
    
  //cout<<"sig eff before graphs is "<<sysErrs["sigeff"].err<<endl;
  


  StoreSysErrors(sysErrs, gMETRes   , mass, bkgSamples);
  //StoreSysErrors(sysErrs, gMETScale , mass, bkgSamples);
  StoreSysErrors(sysErrs, gPU       , mass, bkgSamples);
  StoreSysErrors(sysErrs, gMuPtRes  , mass, bkgSamples);
  StoreSysErrors(sysErrs, gMuPtScale, mass, bkgSamples);
  StoreSysErrors(sysErrs, gElEnScale, mass, bkgSamples);
  StoreSysErrors(sysErrs, gPDF      , mass, bkgSamples);
  //cout<<"sig eff after graphs is "<<sysErrs["sigeff"].err<<endl;
  StoreSysErrorsTail(sysErrs, bkgSamples, ch);
  //cout<<"sig eff after tails is "<<sysErrs["sigeff"].err<<endl;

  sysErrs["sigeff"].val = yields["sigeff"].val;
  sysErrs["sigeff"].err *= sysErrs["sigeff"].val; //convert to abs err
  //cout<<"sig eff scale is "<<sysErrs["sigeff"].err<<endl;

  Value vBkg(0,0);
  for(int isample=0; isample<(int)bkgSamples.size(); isample++){
    sysErrs[bkgSamples[isample]].val = yields[bkgSamples[isample]].val;
    sysErrs[bkgSamples[isample]].err *= sysErrs[bkgSamples[isample]].val; //convert to abs err

    vBkg += sysErrs[bkgSamples[isample]];
  }
  sysErrs["bkg"] = vBkg;

  return sysErrs;
}

map<string, Value>
getSysErrSum(map<string, Value> sysErrs[nch], const vector<string> & bkgSamples){
  map<string, Value> sysErrsSum;
  for(int ch=0; ch<nch; ++ch){
    sysErrsSum["sigeff"] += sysErrs[ch]["sigeff"];
    for(int isample=0; isample<(int)bkgSamples.size(); isample++){
      sysErrsSum[bkgSamples[isample]] += sysErrs[ch][bkgSamples[isample]];
    }  
    sysErrsSum["bkg"] += sysErrs[ch]["bkg"];
  }
  return sysErrsSum;
}

void
StoreSysErrors(map<string, Value> & sysErrs, gPair gSys, const float mass, const vector<string> & bkgSamples){
  const int nbkg = bkgSamples.size();

  sysErrs["sigeff"].err = AddInQuad(sysErrs["sigeff"].err, gSys.first->Eval(mass));
  for(int isample=0; isample<nbkg; isample++){
    sysErrs[bkgSamples[isample]].err = AddInQuad(sysErrs[bkgSamples[isample]].err, gSys.second->Eval(mass));
  }
}

void
StoreSysErrorsTail(map<string, Value> & sysErrs, const vector<string> & bkgSamples, int ch){
  const int nbkg = bkgSamples.size();

  //fout<<"ElTrig   lnN ";
  if(ch<=1) sysErrs["sigeff"].err = AddInQuad(sysErrs["sigeff"].err, 0.02);
  for(int isample=0; isample<nbkg; isample++){
    if(ch<=1) sysErrs[bkgSamples[isample]].err = AddInQuad(sysErrs[bkgSamples[isample]].err, 0.02);
  }
  
  //fout<<"ElReco   lnN ";
  if(ch<=2) sysErrs["sigeff"].err = AddInQuad(sysErrs["sigeff"].err, 0.02);
  for(int isample=0; isample<nbkg; isample++){
    if(ch<=2) sysErrs[bkgSamples[isample]].err = AddInQuad(sysErrs[bkgSamples[isample]].err, 0.02);
  }

  //fout<<"ElIDIso  lnN ";
  if(ch<=2) sysErrs["sigeff"].err = AddInQuad(sysErrs["sigeff"].err, 0.01);
  for(int isample=0; isample<nbkg; isample++){
    if(ch<=2) sysErrs[bkgSamples[isample]].err = AddInQuad(sysErrs[bkgSamples[isample]].err, 0.01);
  }

  //fout<<"MuTrig   lnN ";
  if(ch>=2) sysErrs["sigeff"].err = AddInQuad(sysErrs["sigeff"].err, 0.02);
  for(int isample=0; isample<nbkg; isample++){
    if(ch>=2) sysErrs[bkgSamples[isample]].err = AddInQuad(sysErrs[bkgSamples[isample]].err, 0.02);
  }

  //fout<<"MuReco   lnN ";
  if(ch>=1) sysErrs["sigeff"].err = AddInQuad(sysErrs["sigeff"].err, 0.02);
  for(int isample=0; isample<nbkg; isample++){
    if(ch>=1) sysErrs[bkgSamples[isample]].err = AddInQuad(sysErrs[bkgSamples[isample]].err, 0.02);
  }

  //fout<<"MuIDIso  lnN ";
  if(ch>=1) sysErrs["sigeff"].err = AddInQuad(sysErrs["sigeff"].err, 0.03);
  for(int isample=0; isample<nbkg; isample++){
    if(ch>=1) sysErrs[bkgSamples[isample]].err = AddInQuad(sysErrs[bkgSamples[isample]].err, 0.03);
  }

  //fout<<"#"<<endl;

  //fout<<"MatrixM  lnN ";
  //Nothing for signal
  for(int isample=0; isample<nbkg; isample++){
    if(0 == bkgSamples[isample].compare("TTJets"))     sysErrs[bkgSamples[isample]].err = AddInQuad(sysErrs[bkgSamples[isample]].err, 0.15);
    if(0 == bkgSamples[isample].compare("DYJetsToLL")) sysErrs[bkgSamples[isample]].err = AddInQuad(sysErrs[bkgSamples[isample]].err, 0.30);
  }

  //fout<<"ZGxsec   lnN ";
  //Nothing for signal
  for(int isample=0; isample<nbkg; isample++){
    if(0 == bkgSamples[isample].compare("GVJets")) sysErrs[bkgSamples[isample]].err = AddInQuad(sysErrs[bkgSamples[isample]].err, 0.054);
  }
  
  //fout<<"ZZxsec   lnN ";
  //Nothing for signal
  for(int isample=0; isample<nbkg; isample++){
    if(0 == bkgSamples[isample].compare("ZZ")) sysErrs[bkgSamples[isample]].err = AddInQuad(sysErrs[bkgSamples[isample]].err, 0.30);
  }

  //fout<<"WZxsec   lnN ";
  //Nothing for signal
  for(int isample=0; isample<nbkg; isample++){
    if(0 == bkgSamples[isample].compare("WZJetsTo3LNu")) sysErrs[bkgSamples[isample]].err = AddInQuad(sysErrs[bkgSamples[isample]].err, 0.075);
  }

  //fout<<"WZNLO    lnN ";
  //Nothing for signal
  for(int isample=0; isample<nbkg; isample++){
    if(0 == bkgSamples[isample].compare("WZJetsTo3LNu")) sysErrs[bkgSamples[isample]].err = AddInQuad(sysErrs[bkgSamples[isample]].err, 0.10);
  }

  /*//Don't include lumi in sys error
  //fout<<"lumi     lnN ";
  for(int ch=startCh; ch<endCh; ch++){
    sysErrs[bkgSamples[isample]] = AddInQuad(sysErrs[bkgSamples[isample]], 0.044);//for signal
    for(int isample=0; isample<nbkg; isample++){
      sysErrs[bkgSamples[isample]] = AddInQuad(sysErrs[bkgSamples[isample]], 0.044);
    }
  }
  */
}


void
PrintCard(vector<map<string, Object> > & yields, const string & outName, const vector<string> & bkgSamples ){
  ofstream fOut(outName.c_str());
  if(!fOut) { 
    cout << "Cannot open file " << outName << endl; 
    abort();
  } 

  int nChUsed = yields.size();

  PrintCardHead(fOut, nChUsed);
  /*
  for(int ch=0; ch<nChUsed; ch++){
    fOut<<"# ch="<<ch
        <<"  lumi="<<lumi
        <<"  sigEff="<<max(0., gEffOut[ch]->Eval(mass))
        <<"  sigEffErr="<<max(0., gEffOutErr[ch]->Eval(mass))
        <<"  sigxsec="<<gXSecCh->Eval(mass)
        <<"  nSigEvt="<<lumi * gXSecCh->Eval(mass) * max(0., gEffOut[ch]->Eval(mass))
        <<"  nSigEvtErr="<<1. + (gEffOut[ch]->Eval(mass) > 0. ? max(0., gEffOutErr[ch]->Eval(mass) / gEffOut[ch]->Eval(mass)) : 0.)
        <<endl;
  }
  */

  //Print Data yields
  fOut<<"bin           ";
  for(int ch=0; ch<nChUsed; ch++) 
    fOut<<bin(ch)<<"   ";
  fOut<<endl;

  fOut<<"observation   ";
  for(int ch=0; ch<nChUsed; ch++) 
    fOut<<yields[ch]["data"].nEvt<<"  ";
  fOut<<endl;

  fOut<<"------------"<<endl;

  //Print channel (bin) labels
  fOut<<"bin          "; 
  for(int ch=0; ch<nChUsed; ch++) 
    for(int isample=0; isample<(int)bkgSamples.size()+1; isample++) 
      fOut<<bin(ch)<<"  ";
  fOut<<endl;

    //Print Sample name eg  sig   wz    zjet  ttbar zz    zg 
  fOut<<"process      ";
  for(int ch=0; ch<nChUsed; ch++){
    fOut<<"sig   ";
    for(int isample=0; isample<(int)bkgSamples.size(); isample++)
      fOut<<bkgSamples[isample].substr(0,3)<<"  ";
  }
  fOut<<endl;
    
  //Print Bin number eg  0     1     2     3     4     5
  fOut<<"process      ";
  for(int ch=0; ch<nChUsed; ch++) 
    for(int isample=0; isample<(int)bkgSamples.size()+1; isample++)
      fOut<<isample<<"    ";
  fOut<<endl;
  
  //Print Sig/Bkg yields eg     20.   90    1     1     3     1 
  fOut<<"rate         ";
  for(int ch=0; ch<nChUsed; ch++){
    fOut<<yields[ch]["sig"].nEvt<<"  ";
    for(int isample=0; isample<(int)bkgSamples.size(); isample++)
      fOut<<yields[ch][bkgSamples[isample]].nEvt<<"  ";
  }
  fOut<<endl;

  fOut<<"------------"<<endl;

  //Print Stat Errors
  fOut<<"MCStat    lnN ";
  for(int ch=0; ch<nChUsed; ch++){
    fOut<<Value(yields[ch]["sig"].relStatErr(), -4)<<"  ";
    for(int isample=0; isample<(int)bkgSamples.size(); isample++)
      fOut<<Value(yields[ch][bkgSamples[isample]].relStatErr(), -4)<<"  ";
  }
  fOut<<endl;

  //Print Sys Errors
  for (int iSys=0; iSys<(int)sysSrc.size(); iSys++){//Loop over sys errors
    fOut<<sysSrc[iSys]<<" lnN ";
    for(int ch=0; ch<nChUsed; ch++){
      fOut<<Value(yields[ch]["sig"].relSysErr(sysSrc[iSys]), -4)<<"  ";
      for(int isample=0; isample<(int)bkgSamples.size(); isample++)
        fOut<<Value(yields[ch][bkgSamples[isample]].relSysErr(sysSrc[iSys]), -4)<<"  ";
    }
    fOut<<endl;
  }
  
  fOut.close(); 
}
