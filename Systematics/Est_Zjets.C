//Usage: root -b -l -q 'Est_Zjets.C+(file, useElec, useData)'
//eg:  root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 0, 1)'

#include <vector>
#include "THStack.h"
#include "TStyle.h"
#include "TLine.h"
#include "../Limits/consts.h"

void
Est_Zjets(string infile, bool useElec, bool useData=true, int etaMode=-1, float minpt=-1, float maxpt=9e9){  
  TFile *f = TFile::Open(infile.c_str(), "read"); assert(f);

  vector<string> allsamples;
  if(!useData){
    allsamples.push_back("WJetsToLNu");
    allsamples.push_back("TTJets");
    allsamples.push_back("ZZ");
    allsamples.push_back("GVJets");
    allsamples.push_back("WWTo2L2Nu");
    allsamples.push_back("WZJetsTo3LNu");
    allsamples.push_back("DYJetsToLL");
  }else{
    allsamples.push_back("data");

  }
  TH1F* allTThist=NULL;
  TH1F* allTFhist=NULL;

  string histTTName, histTFName;

  if(etaMode == -1){
    if(useElec){
      allTThist = get_sum_of_hists(f, allsamples,"hZeeMassTT_AllCuts");
      allTFhist = get_sum_of_hists(f, allsamples,"hZeeMassTF_AllCuts");
    }else{
      allTThist = get_sum_of_hists(f, allsamples,"hZmmMassTT_AllCuts");
      allTFhist = get_sum_of_hists(f, allsamples,"hZmmMassTF_AllCuts");
    }
  }else{
    if(useElec){
      histTTName = etaMode ? "hPtVsZeeBarrelMassTT" : "hPtVsZeeEndCapMassTT";
      histTFName = etaMode ? "hPtVsZeeBarrelMassTF" : "hPtVsZeeEndCapMassTF";
    }else{
      histTTName = "hPtVsZmmMassTT";
      histTFName = "hPtVsZmmMassTF";
    }
    histTTName += "_AllCuts";
    histTFName += "_AllCuts";

    TH2F* histTT = get_sum_of_hists2D(f, allsamples, histTTName, 0);
    TH2F* histTF = get_sum_of_hists2D(f, allsamples, histTFName, 0);

    int yminbin,ymaxbin;
    yminbin = histTT->GetYaxis()->FindBin(minpt);
    ymaxbin = histTT->GetYaxis()->FindBin(maxpt)-1;
        
    allTThist = (TH1F*) histTT->ProjectionX("histTT", yminbin, ymaxbin, "e");
    allTFhist = (TH1F*) histTF->ProjectionX("histTF", yminbin, ymaxbin, "e");

    allTThist->Scale(0.5);//Hack bc I forgot to weight TT events by half
  }


  TF1* linearTT = new TF1("linearTT", fline, FitWind_low, FitWind_high, 2);
  TF1* linearTF = new TF1("linearTF", fline, FitWind_low, FitWind_high, 2);

  if(0) cout<<"Fitting Background TT"<<endl;
  TCanvas* c1 = new TCanvas("c1", "Z Mass TT");
  reject = true;
  allTThist->Fit(linearTT, "MRELLQ");
  reject = false;
  double       B_TT = linearTT->Integral     (ZWind_low,ZWind_high);
  double sigma_B_TT = linearTT->IntegralError(ZWind_low,ZWind_high);
  c1->Print(Form("Eff_TT_useElec%i_useData%i.pdf", useElec, useData));

  if(0) cout<<"Fitting Background TF"<<endl;
  TCanvas* c2 = new TCanvas("c2", "Z Mass TF");
  reject = true;
  allTFhist->Fit(linearTF, "MRELLQ");
  reject = false;
  double       B_TF = linearTF->Integral     (ZWind_low,ZWind_high);
  double sigma_B_TF = linearTF->IntegralError(ZWind_low,ZWind_high);
  c2->Print(Form("Eff_TF_useElec%i_useData%i.pdf", useElec, useData));
              
  int  lowBin = allTThist->FindBin(ZWind_low );  
  int highBin = allTThist->FindBin(ZWind_high);

  double sigma_N_TT = 0;
  double sigma_N_TF = 0;
  double N_TT = allTThist->IntegralAndError(lowBin,highBin,sigma_N_TT);
  double N_TF = allTFhist->IntegralAndError(lowBin,highBin,sigma_N_TF);

  ///Calculate Eff////////////////////

  double S_TT = N_TT - B_TT;
  double S_TF = N_TF - B_TF;

  double num = 2*(S_TT);
  double denom = (S_TF) + 2*(S_TT);
  double e_tight = num / denom;

  double sigma2_TT = sigma_N_TT*sigma_N_TT + sigma_B_TT*sigma_B_TT;
  double sigma2_TF = sigma_N_TF*sigma_N_TF + sigma_B_TF*sigma_B_TF;

  double term1 = 4 * pow(S_TF,2) * sigma2_TT;
  double term2 = 4 * pow(S_TT,2) * sigma2_TF;
  double sigma_e = sqrt(term1 + term2) / (denom*denom);

  cout<<"\n\nUsing ";
  if(useElec) cout<<"Electrons";
  else        cout<<"Muons";
  if(useData) cout<<" from Data"<<endl;
  else        cout<<" from MC"<<endl;
  if(etaMode != -1) cout<<"etaMode="<<etaMode<<" from pt="<<minpt<<"->"<<maxpt<<endl;
  cout<<"  Bkg TT:"<<B_TT<<" +/- " <<sigma_B_TT<<endl
      <<"  Bkg TF:"<<B_TF<<" +/- " <<sigma_B_TF<<endl
      <<"  Tot TT:"<<N_TT<<" +/- " <<sigma_N_TT<<endl
      <<"  Tot TF:"<<N_TF<<" +/- " <<sigma_N_TF<<endl
      <<"  Sig TT:"<<S_TT<<endl
      <<"  Sig TF:"<<S_TF<<endl;
  printf("  e_tight: %.4f  +/- %.4f\n", e_tight,sigma_e);

  //N_loose and N_tight are only observables
  //N_loose = N_lep + N_jet
  //N_tight = e_tight * N_lep + P_fake * N_jet
  //get e_tight and p_fake from data
  //e_tight = 2*(S_TT) / [ (S_TF) + 2*(S_TT) ]
  //doesn't work: P_fake = (N_tight - e_tight*N_lep) / (N_loose - N_lep)
  //i think it does work since there is an assumption which is pretty true


}

