#ifndef _common_h_
#define _common_h_

#include "TH1F.h"
#include "TTree.h"

TH1F* get_sum_of_hists(TFile* f, const std::vector<std::string> & samples,
                       const std::string& objName, int rebinme=0, float weight=1.){
  TH1F* hall=NULL;
  const int dim = samples.size();
  if(dim == 0) return NULL;

  for (unsigned j=0; j != samples.size(); ++j) {
    std::string histname = samples[j] + "/" + objName;
    TH1F* h = (TH1F*)f->Get(histname.c_str());
    TH1F* hist = (TH1F*)h->Clone("hist");

    if(hist == NULL){
      cout<<"Failed Getting "<<histname<<endl;
      abort();
    }

    if(!hist->GetSumw2N()) hist->Sumw2();
    if (rebinme){
      float binwidth = hist->GetBinWidth(1);
      string oldbin = Form("Events / %.0f", binwidth);
      string newbin = Form("Events / %.0f", binwidth*rebinme);

      hist->Rebin(rebinme);

      string title = hist->GetYaxis()->GetTitle();
      string::size_type pos = title.find(oldbin);
      title.replace( pos, oldbin.size(), newbin );
      hist->SetYTitle(title.c_str());
      
    }
    if(j==0) hall = (TH1F*)hist->Clone("hall");
    else     hall->Add(hist);

  }
  hall->Scale(weight);
  return hall;
}

////
TH1F* get_sum_of_hists(TFile* f, const std::vector<std::string> & samples,
                       const std::string& objName, const std::string& variable,
                       const std::string& cuts, const std::string& histParam){
  TH1F* hall = NULL;
  for(unsigned j=0; j<samples.size(); ++j){
    std::string fullName = samples[j] + "/" + objName; 
    TTree* tree = (TTree*) f->Get(fullName.c_str()); assert(tree != NULL);
    tree->Draw(Form("%s>>hist%s", variable.c_str(), histParam.c_str()),cuts.c_str());
    TH1F* hist = (TH1F*)gDirectory->Get("hist");
    if(j==0) hall = (TH1F*)hist->Clone("hall");
    else     hall->Add(hist);
    delete tree;
  }
  return hall;
}
//

float
GetLumiUsed(TFile* f){
  TH1F* hist = (TH1F*) f->Get("lumi_ipb");
  return hist ? hist->GetBinContent(1) : 0;
}

#endif//#define _common_h_
