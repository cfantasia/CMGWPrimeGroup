#ifndef _common_h_
#define _common_h_

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

TH1F* get_sum_of_hists(TFile* f, const std::vector<std::string> & samples,
                       const std::string& objName, int rebinme=0, float weight=1.){
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
    if(j==0) hall = (TH1F*)hist->Clone("hall");
    else     hall->Add(hist);

  }
  hall->Scale(weight);
  return hall;
}

void get_sum_of_hists(TFile* f, const std::vector<std::string> & samples,
                      const std::string& objName, const std::string& variable,
                      const std::string& cuts, TH1F & hist){
  for(unsigned j=0; j<samples.size(); ++j){
    std::string fullName = samples[j] + "/" + objName; //std::cout<<"tree name is "<<fullName<<std::endl;
    TTree* tree = (TTree*) f->Get(fullName.c_str()); assert(tree != NULL);
    tree->Draw(Form("weight:%s",variable.c_str()), cuts.c_str(), "goff");
    int n = tree->GetSelectedRows();
    for(int ientry=0; ientry<n; ++ientry){
      const float weight = tree->GetVal(0)[ientry];
      const float var = tree->GetVal(1)[ientry];
      hist.Fill(var, weight);
    }
    delete tree;
  }
}

void get_sum_of_hists(TFile* f, const std::vector<std::string> & samples,
                      const std::string& objName, const std::string& variable,
                      const std::string& cuts, TH2 & hist){
  for(unsigned j=0; j<samples.size(); ++j){
    std::string fullName = samples[j] + "/" + objName; 
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
  if(!hLumi) cout<<" Failed getting hLumi"<<endl;
  return valid ? hLumi->GetBinContent(1) / hLumi->GetBinContent(2) : -1;
}

#endif//#define _common_h_
