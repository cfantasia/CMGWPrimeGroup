#ifndef _common_h_
#define _common_h_

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include <iomanip>
#include <iostream>

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

};

TH1F* get_sum_of_hists(TFile* f, const std::vector<std::string> & samples,
                       const std::string& objName, int rebinme, float weight){
  vector<float> weights(samples.size(), weight);
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

Value
GetNEvtsAndError(TTree* tree, const string & cuts){
  tree->Draw("weight", cuts.c_str(), "goff");
  int n = tree->GetSelectedRows();
  TH1F hist("hist", "Dummy Hist", 1, 0, 10);
  hist.Sumw2();
  for(int ientry=0; ientry<n; ++ientry){
    float weight = tree->GetW()[ientry];
    hist.Fill(1, weight);
    //printf("weight is %.3f\n", weight);
  }
  //printf("new: %.2f old: %.2f\n", hist.GetBinError(1), sqrt(hist.GetBinContent(1)));
  return Value(hist.GetBinContent(1), hist.GetBinError(1));
}

float
GetNEvts(TTree* tree, const string & cuts){
  Value val = GetNEvtsAndError(tree, cuts);
  return val.val;
}

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

#endif//#define _common_h_
