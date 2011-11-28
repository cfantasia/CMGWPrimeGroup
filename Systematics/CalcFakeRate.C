//Usage: root -b -q 'CalcFakeRate.C+(file, useData)'

#include <vector>
#include "THStack.h"
#include "TStyle.h"
#include "TLine.h"
#include "../Limits/consts.h"

void
CalcFakeRate(string infile, bool useData=true){  
  TFile *f = TFile::Open(infile.c_str(), "read"); assert(f);

  vector<string> allsamples;
  if(!useData){
    allsamples.push_back("WJetsToLNu");//corrupted file
    allsamples.push_back("TTJets");
    allsamples.push_back("ZZ");
    allsamples.push_back("GVJets");
    allsamples.push_back("WWTo2L2Nu");
    allsamples.push_back("WZJetsTo3LNu");
    allsamples.push_back("DYJetsToLL");
  }else{
    allsamples.push_back("data");
  }
  
  float in_tot(0), out_tot(0);
  for(unsigned i=0; i<allsamples.size(); ++i){
    string hist_name = allsamples[i] + "/hNumEvts";
    TH1F* hist = (TH1F*) f->Get(hist_name.c_str()); assert(hist);
    int lastbin = hist->GetNbinsX();
    

    float in  = hist->GetBinContent(lastbin);
    float out = hist->GetBinContent(lastbin-1);

    in_tot += in;
    out_tot += out;

    cout<<"Sample: "<<allsamples[i]<<" = "<<in/out<<" = pass/total = "<<in<<"/"<<out<<endl;
  }

  cout<<"Total: "<<in_tot/out_tot*100<<"% = in_tot/out_tot*100% = "<<in_tot<<"/"<<out_tot<<"*100%\n";
}
