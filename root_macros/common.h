#ifndef _common_h_
#define _common_h_

#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TTree.h"
#include <iomanip>
#include <iostream>
#include "TLegend.h"
#include "TMath.h"

struct Value{
  float val;
  float err;
  bool printErr;
  int forcePrecision;

  Value():val(0),err(0),printErr(false), forcePrecision(0){}//Let's not mess up log10
  Value(float v):val(v),err(v),printErr(false), forcePrecision(0){}
  Value(float v, float e, bool b=false):val(v),err(e),printErr(b), forcePrecision(0){
    //if(e < 0.) err = sqrt(val);
  }
  Value(Value v, float f){
    val = v.val;
    err = f < 0. ? sqrt(val) : f;
    printErr = v.printErr;
    forcePrecision = v.forcePrecision;
  }
  void setPrintErr(bool b){
    printErr = b;
  }
  int findPrecision() const{
    if(forcePrecision) return forcePrecision;
    if(err < 0.){
      int nExtra = abs((int)err);
      Value temp(val);
      return temp.findPrecision()+nExtra-1;//-1 is since we always print 1 sig fig
    }else if(err == 0) return 0;
    return std::max(-1.*floor(log10(err)), 0.); 
  }

  double relErr() const{
    return 1. + (val > 0. ? err / val : 0.);
  }

  friend std::ostream& operator << (std::ostream &o, const Value & v){
    ios_base::fmtflags flags = o.flags();//Get Current config
    o<<setiosflags(ios::fixed) << setprecision(v.findPrecision())<<v.val;
    if(v.printErr) o<<" +/- "<<v.err;
    o.flags(flags);//Reset config
    return o;
  }

  Value operator+( const Value &rhs ) const{
    Value result = *this;     // Make a copy of myself.  Same as Value result(*this);
    result += rhs;            // Use += to add other to the copy.
    return result;
  }
  Value & operator+=(const Value &rhs) {
    val += rhs.val;
    err = sqrt(pow(err,2) + pow(rhs.err,2));
    return *this;
  }
  Value operator-( const Value &rhs ) const{
    Value result = *this;     // Make a copy of myself.  Same as Value result(*this);
    result -= rhs;            // Use -= to add other to the copy.
    return result;
  }
  Value & operator-=(const Value &rhs) {
    val -= rhs.val;
    err = sqrt(pow(err,2) - pow(rhs.err,2));//Cory: Wrong?  What should it be?
    return *this;
  }
  Value operator*( const double &rhs ) const{
    Value result = *this;     
    result *= rhs;            
    return result;
  }
  Value & operator*=(const double &rhs) {
    val *= rhs;
    err *= rhs;
    return *this;
  }
  
};

float
AddInQuad(float a, float b){
  return sqrt(a*a + b*b);
}

float WZScaleFactor(int EvtType, float ZLep1Pt, float ZLep1Eta, float ZLep2Pt, float ZLep2Eta, float WLepPt, float WLepEta);

/*
TH1F* get_sum_of_hists(TFile* f, const std::string & sample,
                       const std::string& objName, int rebinme, float weight){
  const std::vector<std::string> samples(1, sample);
  return get_sum_of_hists(f, samples, objName, rebinme, weight);
}
*/
TH1F* get_sum_of_hists(TFile* f, const std::vector<std::string> & samples,
                       const std::string& objName, int rebinme, float weight){
  std::vector<float> weights(samples.size(), weight);
  return get_sum_of_hists(f, samples, objName, rebinme, weight);
}

TH1F* get_sum_of_hists(TFile* f, const std::vector<std::string> & samples,
                       const std::string& objName, int rebinme=0, const std::vector<float> & weights=std::vector<float>()){
  TH1F* hall=NULL;
  const int dim = samples.size();
  if(dim == 0) return NULL;

  for (unsigned j=0; j != samples.size(); ++j) {
    std::string histname = samples[j] + "/" + objName;
    TH1F* h = (TH1F*)f->Get(histname.c_str());
    if(h == NULL){
      std::cout<<"Failed Getting "<<histname<<std::endl;
      abort();
    }
    TH1F* hist = (TH1F*)h->Clone("hist");


    if(!hist->GetSumw2N()) hist->Sumw2();
    if (rebinme){
      float binwidth = hist->GetBinWidth(1);
      std::string oldbin = Form("Events / %.0f", binwidth);
      std::string newbin = Form("Events / %.0f", binwidth*rebinme);

      hist->Rebin(rebinme);

      std::string title = hist->GetYaxis()->GetTitle();
      std::string::size_type pos = title.find(oldbin);
      title.replace( pos, oldbin.size(), newbin );
      hist->SetYTitle(title.c_str());
      
    }
    //if(weights.size() == samples.size()) printf("scaling %s by %.2f!!!\n", samples[j].c_str(), weights[j]);
    if(weights.size() == samples.size()) hist->Scale(weights[j]);

    if(j==0) hall = (TH1F*)hist->Clone("hall");
    else     hall->Add(hist);

  }
  return hall;
}

TH2F* get_sum_of_hists2D(TFile* f, const std::vector<std::string> & samples,
                         const std::string& objName, int rebinme=0, float weight=1.){
  TH2F* hall=NULL;
  const int dim = samples.size();
  if(dim == 0) return NULL;

  for (unsigned j=0; j != samples.size(); ++j) {
    std::string histname = samples[j] + "/" + objName;
    TH2F* h = (TH2F*)f->Get(histname.c_str());
    if(h == NULL){
      std::cout<<"Failed Getting "<<histname<<std::endl;
      abort();
    }
    TH2F* hist = (TH2F*)h->Clone("hist");


    if(!hist->GetSumw2N()) hist->Sumw2();
    if (rebinme){
      float binwidth = hist->GetBinWidth(1);
      std::string oldbin = Form("Events / %.0f", binwidth);
      std::string newbin = Form("Events / %.0f", binwidth*rebinme);

      hist->Rebin(rebinme);

      std::string title = hist->GetYaxis()->GetTitle();
      std::string::size_type pos = title.find(oldbin);
      title.replace( pos, oldbin.size(), newbin );
      hist->SetYTitle(title.c_str());
      
    }
    if(j==0) hall = (TH2F*)hist->Clone("hall");
    else     hall->Add(hist);

  }
  hall->Scale(weight);
  return hall;
}

void get_sum_of_hists(TFile* f, const std::vector<std::string> & samples,
                      const std::string& objName, const std::string& variable,
                      const std::string& cuts, TH1F & hist, 
                      const std::vector<float> & weights=std::vector<float>()){
  hist.Sumw2();
  bool doWeight = weights.size() == samples.size();
  for(unsigned j=0; j<samples.size(); ++j){
    std::string fullName = samples[j] + "/" + objName; 
    TTree* tree = (TTree*) f->Get(fullName.c_str()); 
    if(!tree){
      std::cout<<"Failed getting tree named "<<fullName<<std::endl;
      assert(tree != NULL);
    }
    tree->Draw(variable.c_str(), cuts.c_str(), "goff");
    int n = tree->GetSelectedRows();
    for(int ientry=0; ientry<n; ++ientry){
      float weight = tree->GetW()[ientry];
      if(doWeight) weight *= weights[j];
      const float var = tree->GetVal(0)[ientry];
      hist.Fill(var, weight);
    }
    delete tree;
  }
}

void get_sum_of_hists(TFile* f, const std::vector<std::string> & samples,
                      const std::string& objName, const std::string& variable,
                      const std::string& cuts, TH2 & hist){
  hist.Sumw2();
  for(unsigned j=0; j<samples.size(); ++j){
    std::string fullName = samples[j] + "/" + objName; 
    //printf("Looking for %s\n", fullName.c_str());
    TTree* tree = (TTree*) f->Get(fullName.c_str()); assert(tree != NULL);
    tree->Draw(Form("weight:%s",variable.c_str()), cuts.c_str(), "goff");
    int n = tree->GetSelectedRows();
    for(int ientry=0; ientry<n; ++ientry){
      const float weight = tree->GetVal(0)[ientry];
      const float varX = tree->GetVal(1)[ientry];
      const float varY = tree->GetVal(2)[ientry];
      hist.Fill(varX, varY, weight);
    }
    delete tree;
  }
}


void get_sum_of_hists(TFile* f, const std::string & sample,
                      const std::string& objName, const std::string& variable,
                      const std::string& cuts, TH1F & hist){
  const std::vector<string> samples(1, sample);
  get_sum_of_hists(f, samples, objName, variable, cuts, hist);
}

float
GetLumiUsed(TFile* f){
  TH1F* hLumi = (TH1F*) f->Get("lumi_ipb");
  bool valid = hLumi;
  if(!hLumi) std::cout<<" Failed getting hLumi"<<std::endl;
  return valid ? hLumi->GetBinContent(1) / hLumi->GetBinContent(2) : -1;
}

float
GetSampleInfo(TH1F* h, std::string binName){
  TAxis* axis = h->GetXaxis();
  int bin = axis->FindBin(binName.c_str()); 
  if(bin == 0) std::cout<<"Didn't find bin with name "<<binName<<std::endl;
  //printf("bin named %s is bin number %i\n", binName.c_str(), bin);
  float value = h->GetBinContent(bin);
  //printf("value: %f \n", value);
  //Most parameters are not to be added (so average)
  if(binName.find("Number of Events in root Files") == std::string::npos &&
     binName.find("Number of Files Merged")         == std::string::npos ){
    int nMerged = h->GetBinContent(axis->FindBin("Number of Files Merged"));
    //printf("value: %f and nMerged: %i\n", value, nMerged);
    value /= nMerged;
  }
  return value;
}


////////////
TTree* getTree(TFile* f, const vector<string> & samples, const string & tName){
  gROOT->cd();
  TList list;
  for(unsigned i=0; i<samples.size(); ++i){
    TTree* t = (TTree*) f->Get(Form("%s/%s",samples[i].c_str(), tName.c_str()));
    if(!t)  cout<<f->GetName()<<" "<<samples[i]<<" "<<tName<<endl;
    assert(t);
    list.Add(t);
  }
  return TTree::MergeTrees(&list);
}

TTree* getTree(TFile* f, const string & sample, const string & tName){
  const vector<string> samples(1, sample);
  return getTree(f, samples, tName);
}

//////////////////
Value
GetNEvtsAndError(TTree* tree, const string & cuts, const bool useScaleFactors=false){
  if(useScaleFactors) tree->Draw("EvtType:ZLep1Pt:ZLep1Eta:ZLep2Pt:ZLep2Eta:WLepPt:WLepEta", cuts.c_str(), "para goff");
  else                tree->Draw("weight", cuts.c_str(), "para goff");
  int n = tree->GetSelectedRows();
  if(n>0 && tree->GetW() == NULL) abort();
  TH1F hist("hist", "Dummy Hist", 1, 0, 10);
  hist.Sumw2();
  for(int ientry=0; ientry<n; ++ientry){
    float weight = tree->GetW()[ientry];
    if(useScaleFactors){
      int treeIdx = 0;
      const int    EvtType  = round(tree->GetVal(treeIdx++)[ientry]);
      const double ZLep1Pt  = tree->GetVal(treeIdx++)[ientry];
      const double ZLep1Eta = tree->GetVal(treeIdx++)[ientry];
      const double ZLep2Pt  = tree->GetVal(treeIdx++)[ientry];
      const double ZLep2Eta = tree->GetVal(treeIdx++)[ientry];
      const double WLepPt   = tree->GetVal(treeIdx++)[ientry];
      const double WLepEta  = tree->GetVal(treeIdx++)[ientry];

      weight *= WZScaleFactor(EvtType, ZLep1Pt,ZLep1Eta,ZLep2Pt,ZLep2Eta,WLepPt,WLepEta);
    }
    hist.Fill(1, weight);
    //printf("weight is %.3f\n", weight);
  }
  //printf("new: %.2f old: %.2f\n", hist.GetBinError(1), sqrt(hist.GetBinContent(1)));
  //printf("val: %.4f err: %.4f\n", hist.GetBinContent(1), hist.GetBinError(1));
  return Value(hist.GetBinContent(1), hist.GetBinError(1));
}

Value
GetNEvtsAndError(TFile* f, const vector<string> & samples, const string & tName, const string & cuts, const bool useScaleFactors=false){
  TTree* t = getTree(f, samples, tName);
  return GetNEvtsAndError(t, cuts, useScaleFactors);
}

Value
GetNEvtsAndError(TFile* f, const string & sample, const string & tName, const string & cuts, const bool useScaleFactors=false){
  const vector<string> samples(1, sample);
  return GetNEvtsAndError(f, samples, tName, cuts, useScaleFactors);
}

float
GetNEvts(TTree* tree, const string & cuts, const bool useScaleFactors=false){
  Value val = GetNEvtsAndError(tree, cuts, useScaleFactors); 
  return val.val;
}

float
GetNEvts(TFile* f, const vector<string> & samples, const string & tName, const string & cuts, const bool useScaleFactors=false){
  TTree* t = getTree(f, samples, tName);
  assert(t);
  return GetNEvts(t, cuts, useScaleFactors);
}

float
GetNEvts(TFile* f, const string & sample, const string & tName, const string & cuts, const bool useScaleFactors=false){
  const vector<string> samples(1, sample);
  return GetNEvts(f, samples, tName, cuts, useScaleFactors);
}

Value
ShiftErr(const Value & orig, const Value & mod, const double & corrCoef=0.){
  Value ret;
  if(orig.val){
    ret.val = fabs(1. - (mod.val / orig.val));
    ret.err = (mod.val / orig.val) * sqrt(pow(mod.err/mod.val,2) + pow(orig.err/orig.val,2) - 2 * (mod.err/mod.val) * (orig.err/orig.val) * corrCoef);
  }else{
    ret.val = 1.;
    ret.err = 0.;
  }
  return ret;
}

float
CalcCorrCoef(TTree* t1, TTree* t2, const string & cuts){
  TH1F hTT("hTT", "", 1, 0, 10), hTF("hTF", "", 1, 0, 10), hFT("hFT", "", 1, 0, 10);
  //cout<<"cuts are "<<cuts<<endl;
  int n = t1->Draw("Run:Lumi:Event", cuts.c_str(), "goff");
  for(int ientry=0; ientry<n; ++ientry){
    int idx = 0;
    unsigned long int run    = round(t1->GetVal(idx++)[ientry]);
    unsigned long int lumi   = round(t1->GetVal(idx++)[ientry]);
    unsigned long int event  = round(t1->GetVal(idx++)[ientry]);
    float weight = t1->GetW()[ientry];
    
    string newCuts = cuts + Form("*(Run==%lu && Lumi==%lu && Event==%lu)", run, lumi, event);
    if(t2->Draw("Run:Lumi:Event", newCuts.c_str(), "goff")){
      hTT.Fill(1, weight);
    }else{
      hTF.Fill(1, weight);
    }
  }

  n = t2->Draw("Run:Lumi:Event", cuts.c_str(), "goff");
  for(int ientry=0; ientry<n; ++ientry){
    int idx = 0;
    unsigned long int run    = round(t2->GetVal(idx++)[ientry]);
    unsigned long int lumi   = round(t2->GetVal(idx++)[ientry]);
    unsigned long int event  = round(t2->GetVal(idx++)[ientry]);
    float weight = t2->GetW()[ientry];
    
    string newCuts = cuts + Form("*(Run==%lu && Lumi==%lu && Event==%lu)", run, lumi, event);
    if(t1->Draw("Run:Lumi:Event", newCuts.c_str(), "goff")){
      //Remember, we've already counted TT events
    }else{
      hFT.Fill(1, weight);
    }
  }

  //printf("TT: %.2f\n",hTT.GetBinContent(1));
  //printf("TF: %.2f\n",hTF.GetBinContent(1));
  //printf("FT: %.2f\n",hFT.GetBinContent(1));

  float sigma_TT = hTT.GetBinError(1);
  float sigma_TF = hTF.GetBinError(1);
  float sigma_FT = hFT.GetBinError(1);

  float sigma_n = AddInQuad(sigma_TT,sigma_TF);
  float sigma_m = AddInQuad(sigma_TT,sigma_FT);

  float retVal = (sigma_TT*sigma_TT) / (sigma_m*sigma_n);
  //printf("corrCoef = %.2f\n", retVal);
  return retVal;
}


float 
roundToNearest(float value, int mark){ //Cool trick
  return round(value/mark)*mark;
}

void
Scale(TGraph* g, double scaleFactor){
  for (int i=0;i<g->GetN();i++) g->GetY()[i] *= scaleFactor;
}

void
prepLegend(TLegend* leg){
  leg->SetTextSize(0.05);
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  //leg->SetNColumns(2);
  //leg->SetColumnSeparation(0.05);
}

///WprimeWZ Specific Stuff
float WindowWidth(float mass){
  return max(0.2*mass-20., 1.45*mass - 1500);
  //return max(50., 0.5*mass - 100);
}

float LtCut(float mass){
  return min(500., 0.5*mass+25.);
  //return min(400., 0.6*mass - 50.);
}

string 
AnalysisCuts(float mass, float LtOffset=0., float WindOffset=0.){
    float winWidth = WindowWidth(mass) + WindOffset;
    float minMass = mass - winWidth/2.;
    float maxMass = mass + winWidth/2.;
    float minLt   = LtCut(mass) + LtOffset;
    return Form("WZMass > %.0f && WZMass < %.0f && Lt > %.0f", minMass, maxMass, minLt);
}

string bin(int ch){
  if(ch == 0) return "eee";
  if(ch == 1) return "eem";
  if(ch == 2) return "mme";
  if(ch == 3) return "mmm";
  else abort();
  return "";
}

vector<string> BkgSamples(){
  vector<string> bkgSamples;
  bkgSamples.push_back("WZJetsTo3LNu");
  bkgSamples.push_back("DYJetsToLL");
  bkgSamples.push_back("TTJets");
  bkgSamples.push_back("ZZ");
  bkgSamples.push_back("GVJets");
  return bkgSamples;
}


///////////////////
// Scale Factors //
///////////////////
float laserTurnOnClusEt(float et);
float noLaserTurnOnClusEt(float et);
float barrelTurnOnClusEt(float et);

//float totalClusEtTurnOn(float et,float detEta,float totLumi=19300)
float totalClusEtTurnOn(float et,float detEta,float totLumi=18258)
{
  if(fabs(detEta)<1.5) return barrelTurnOnClusEt(et);
  else{
    float fracPostLaser=(totLumi-325.)/totLumi;
    return fracPostLaser*laserTurnOnClusEt(et) + (1-fracPostLaser)*noLaserTurnOnClusEt(et);
  }

}


float barrelTurnOnClusEt(float et)
{
  float maxEff = 0.996;
  maxEff=1.; //its 1 within errors
  float midPoint = 34.76;
  float turnOn=0.85;
  return 0.5*maxEff*(1+TMath::Erf((et-midPoint)/(sqrt(2)*turnOn)));
}

float laserTurnOnClusEt(float et)
{
  float maxEff = 0.9948;
  maxEff=1.; //its 1 within errors
  float midPoint = 32.74;
  float turnOn=2.36;
  return 0.5*maxEff*(1+TMath::Erf((et-midPoint)/(sqrt(2)*turnOn)));

}


float noLaserTurnOnClusEt(float et)
{
  float maxEff = 0.979;
  maxEff=1.; //its 1 within errors more or less
  float midPoint = 37.2;
  float turnOn=2.3;
  return 0.5*maxEff*(1+TMath::Erf((et-midPoint)/(sqrt(2)*turnOn)));

}

float
ElectronTriggerScaleFactor(float pt1, float eta1, float pt2, float eta2){
  return 0.995*
    totalClusEtTurnOn(pt1, eta1)*
    totalClusEtTurnOn(pt2, eta2);
}

float
ZElectronIDIsoScaleFactor(float pt, float eta){
  if      (eta < 0.8){
    if(pt>50) return 1.007;
    if(pt>40) return 1.007;
    if(pt>30) return 1.003;
    if(pt>20) return 1.004;
  }else if(eta < 1.442){
    if(pt>50) return 0.995;
    if(pt>40) return 0.992;
    if(pt>30) return 0.984;
    if(pt>20) return 0.975;
  }else if(eta < 1.556){
    if(pt>50) return 0.993;
    if(pt>40) return 0.991;
    if(pt>30) return 1.006;
    if(pt>20) return 1.034;
  }else if(eta < 2.0){
    if(pt>50) return 1.007;
    if(pt>40) return 1.006;
    if(pt>30) return 0.990;
    if(pt>20) return 0.983;
  }else if(eta < 2.5){
    if(pt>50) return 1.009;
    if(pt>40) return 1.013;
    if(pt>30) return 1.022;
    if(pt>20) return 1.025;
  }
  abort();
  return 1.;
}

float
ZElectronScaleFactor(float pt, float eta){
  return ZElectronIDIsoScaleFactor(pt, eta);
}

float
WElectronScaleFactor(float pt, float eta){
  if      (eta < 0.8){
    if(pt>50) return 1.008;
    if(pt>40) return 1.008;
    if(pt>30) return 1.004;
    if(pt>20) return 1.005;
  }else if(eta < 1.442){
    if(pt>50) return 0.999;
    if(pt>40) return 0.994;
    if(pt>30) return 0.991;
    if(pt>20) return 0.981;
  }else if(eta < 1.556){
    if(pt>50) return 0.994;
    if(pt>40) return 0.989;
    if(pt>30) return 0.998;
    if(pt>20) return 1.044;
  }else if(eta < 2.0){
    if(pt>50) return 1.006;
    if(pt>40) return 1.004;
    if(pt>30) return 0.992;
    if(pt>20) return 0.980;
  }else if(eta < 2.5){
    if(pt>50) return 1.009;
    if(pt>40) return 1.005;
    if(pt>30) return 1.019;
    if(pt>20) return 1.017;
  }
  abort();
  return 1.;
}

float
MuonTriggerScaleFactor(float pt1, float eta1, float pt2, float eta2){
  
  pt1 += 0.;//supress unused variable warning
  pt2 += 0.;//supress unused variable warning
  float leg1(0.), leg2(0.);
  if      (eta1 < 0.9){
    leg1 = 0.9884;
  }else if(eta1 < 1.2){
    leg1 = 0.9836;
  }else if(eta1 < 2.4){
    leg1 = 0.9879;
  }else{
    abort();
  }
    
  if      (eta2 < 0.9){
    leg2 = 0.9820;
  }else if(eta2 < 1.2){
    leg2 = 0.9836;
  }else if(eta2 < 2.4){
    leg2 = 0.9836;
  }else{
    abort();
  }

  return leg1 * leg2;
  /*
  //Updated for HLT_Mu17_TrkMu8
  if(pt2 > 20){
    if(0.0 < eta1 && eta1 < 0.9 && 0.0 < eta2 && eta2 < 0.9) return 0.98;
    if(0.0 < eta1 && eta1 < 0.9 && 0.9 < eta2 && eta2 < 1.2) return 0.97;
    if(0.0 < eta1 && eta1 < 0.9 && 1.2 < eta2 && eta2 < 2.1) return 0.99;
    if(0.0 < eta1 && eta1 < 0.9 && 2.1 < eta2 && eta2 < 2.4) return 1.01;

    if(0.9 < eta1 && eta1 < 1.2 && 0.0 < eta2 && eta2 < 0.9) return 0.98;
    if(0.9 < eta1 && eta1 < 1.2 && 0.9 < eta2 && eta2 < 1.2) return 1.00;
    if(0.9 < eta1 && eta1 < 1.2 && 1.2 < eta2 && eta2 < 2.1) return 0.99;
    if(0.9 < eta1 && eta1 < 1.2 && 2.1 < eta2 && eta2 < 2.4) return 1.00;

    if(1.2 < eta1 && eta1 < 2.1 && 0.0 < eta2 && eta2 < 0.9) return 1.00;
    if(1.2 < eta1 && eta1 < 2.1 && 0.9 < eta2 && eta2 < 1.2) return 1.01;
    if(1.2 < eta1 && eta1 < 2.1 && 1.2 < eta2 && eta2 < 2.1) return 1.00;
    if(1.2 < eta1 && eta1 < 2.1 && 2.1 < eta2 && eta2 < 2.4) return 0.98;

    if(2.1 < eta1 && eta1 < 2.4 && 0.0 < eta2 && eta2 < 0.9) return 1.00;
    if(2.1 < eta1 && eta1 < 2.4 && 0.9 < eta2 && eta2 < 1.2) return 1.00;
    if(2.1 < eta1 && eta1 < 2.4 && 1.2 < eta2 && eta2 < 2.1) return 1.01;
    if(2.1 < eta1 && eta1 < 2.4 && 2.1 < eta2 && eta2 < 2.4) return 1.07;
  }else{
    if(0.0 < eta1 && eta1 < 0.9 && 0.0 < eta2 && eta2 < 0.9) return 0.99;
    if(0.0 < eta1 && eta1 < 0.9 && 0.9 < eta2 && eta2 < 1.2) return 1.01;
    if(0.0 < eta1 && eta1 < 0.9 && 1.2 < eta2 && eta2 < 2.1) return 1.02;
    if(0.0 < eta1 && eta1 < 0.9 && 2.1 < eta2 && eta2 < 2.4) return 0.96;

    if(0.9 < eta1 && eta1 < 1.2 && 0.0 < eta2 && eta2 < 0.9) return 1.01;
    if(0.9 < eta1 && eta1 < 1.2 && 0.9 < eta2 && eta2 < 1.2) return 1.00;
    if(0.9 < eta1 && eta1 < 1.2 && 1.2 < eta2 && eta2 < 2.1) return 0.95;
    if(0.9 < eta1 && eta1 < 1.2 && 2.1 < eta2 && eta2 < 2.4) return 0.94;

    if(1.2 < eta1 && eta1 < 2.1 && 0.0 < eta2 && eta2 < 0.9) return 1.01;
    if(1.2 < eta1 && eta1 < 2.1 && 0.9 < eta2 && eta2 < 1.2) return 0.98;
    if(1.2 < eta1 && eta1 < 2.1 && 1.2 < eta2 && eta2 < 2.1) return 1.02;
    if(1.2 < eta1 && eta1 < 2.1 && 2.1 < eta2 && eta2 < 2.4) return 1.01;

    if(2.1 < eta1 && eta1 < 2.4 && 0.0 < eta2 && eta2 < 0.9) return 0.96;
    if(2.1 < eta1 && eta1 < 2.4 && 0.9 < eta2 && eta2 < 1.2) return 1.04;
    if(2.1 < eta1 && eta1 < 2.4 && 1.2 < eta2 && eta2 < 2.1) return 0.95;
    if(2.1 < eta1 && eta1 < 2.4 && 2.1 < eta2 && eta2 < 2.4) return 1.04;
  }
  */

  printf("Aborting b/c mu trig sf is wrong. Eta1 = %.2f Eta2 = %.2f\n", eta1, eta2);
  abort();
  return 0.;


}

float
ZMuonIDScaleFactor(float pt, float eta){
  if      (eta < 2.1){
    if(pt>90) return 0.9926;
    if(pt>80) return 0.9883;
    if(pt>70) return 0.9854;
    if(pt>60) return 0.9827;
    if(pt>50) return 0.9871;
    if(pt>40) return 0.9909;
    if(pt>30) return 0.9888;
    if(pt>20) return 0.9876;
    if(pt>10) return 0.9525;
  }else if(eta < 2.4){
    if(pt>90) return 0.9783;
    if(pt>80) return 0.9638;
    if(pt>70) return 0.9841;
    if(pt>60) return 0.9857;
    if(pt>50) return 0.9827;
    if(pt>40) return 0.9922;
    if(pt>30) return 0.9914;
    if(pt>20) return 0.9836;
    if(pt>10) return 0.9049;
  }
  printf("Aborting b/c pt=%.0f and eta=%.2f\n", pt, eta);
  abort();
  return 1.;
}

float
ZMuonIsoScaleFactor(float pt, float eta){
  if      (eta < 2.1){
    if(pt>90) return 0.9977;
    if(pt>80) return 1.0010;
    if(pt>70) return 0.9973;
    if(pt>60) return 0.9951;
    if(pt>50) return 0.9937;
    if(pt>40) return 0.9906;
    if(pt>30) return 0.9863;
    if(pt>20) return 0.9738;
    if(pt>10) return 0.9273;
  }else if(eta < 2.4){
    if(pt>90) return 1.0137;
    if(pt>80) return 1.0163;
    if(pt>70) return 1.0089;
    if(pt>60) return 1.0309;
    if(pt>50) return 1.0444;
    if(pt>40) return 1.0678;
    if(pt>30) return 1.1139;
    if(pt>20) return 1.1617;
    if(pt>10) return 1.1931;
  }
  printf("Aborting b/c pt=%.0f and eta=%.2f\n", pt, eta);
  abort();
  return 1.;
}

float
ZMuonScaleFactor(float pt, float eta){
  return ZMuonIDScaleFactor(pt, eta) * ZMuonIsoScaleFactor(pt, eta);
}

float
WMuonIDScaleFactor(float pt, float eta){
  if      (eta < 0.9){
    //if(pt>45) return 0.9899;
    if(pt>10) return 0.9932;
  }else if(eta < 1.2){
    //if(pt>45) return 0.9884;
    if(pt>10) return 0.9911;
  }else if(eta < 2.1){
    //if(pt>45) return 0.9958;
    if(pt>10) return 0.9975;
  }else if(eta < 2.4){
    //if(pt>45) return 0.9910;
    if(pt>10) return 0.9946;
  }
  abort();
  return 1.;
}

float
WMuonIsoScaleFactor(float pt, float eta){
  if      (eta < 0.9){
    if(pt>10) return 1.0004;
  }else if(eta < 1.2){
    if(pt>10) return 1.0031;
  }else if(eta < 2.4){
    if(pt>10) return 1.0050;//Cory: No factors below 20 GeV so assuming the same
  }
  abort();
  return 1.;
}

float
WMuonScaleFactor(float pt, float eta){
  return WMuonIDScaleFactor(pt, eta) * WMuonIsoScaleFactor(pt, eta);
}

float
WZScaleFactor(int EvtType, float ZLep1Pt, float ZLep1Eta, float ZLep2Pt, float ZLep2Eta, float WLepPt, float WLepEta){
  ZLep1Eta = fabs(ZLep1Eta);
  ZLep2Eta = fabs(ZLep2Eta);
  WLepEta  = fabs(WLepEta);
  
  if(EvtType == 0) return ZElectronScaleFactor(ZLep1Pt, ZLep1Eta)*ElectronTriggerScaleFactor(ZLep1Pt, ZLep1Eta, ZLep2Pt, ZLep2Eta)*
                          ZElectronScaleFactor(ZLep2Pt, ZLep2Eta)*WElectronScaleFactor(WLepPt, WLepEta);

  if(EvtType == 1) return ZElectronScaleFactor(ZLep1Pt, ZLep1Eta)*ElectronTriggerScaleFactor(ZLep1Pt, ZLep1Eta, ZLep2Pt, ZLep2Eta)*
                          ZElectronScaleFactor(ZLep2Pt, ZLep2Eta)*    WMuonScaleFactor(WLepPt, WLepEta);

  if(EvtType == 2) return ZMuonScaleFactor    (ZLep1Pt, ZLep1Eta)*    MuonTriggerScaleFactor(ZLep1Pt, ZLep1Eta, ZLep2Pt, ZLep2Eta)*
                          ZMuonScaleFactor(ZLep2Pt, ZLep2Eta)*WElectronScaleFactor(WLepPt, WLepEta);

  if(EvtType == 3) return ZMuonScaleFactor    (ZLep1Pt, ZLep1Eta)*    MuonTriggerScaleFactor(ZLep1Pt, ZLep1Eta, ZLep2Pt, ZLep2Eta)*
                          ZMuonScaleFactor(ZLep2Pt, ZLep2Eta)*    WMuonScaleFactor(WLepPt, WLepEta);

  std::cout<<"Aborting b/c Event Type is "<<EvtType<<std::endl;
  abort();
  return 0.;
}

#endif//#define _common_h_
