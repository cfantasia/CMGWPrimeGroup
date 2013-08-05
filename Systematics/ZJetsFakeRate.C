//Author: Cory Fantasia 2012
//Purpose: Calculate fake rates on z+jets sample
//Usage: root -b -l -q 'ZJetsFakeRate.C+(file, useData, doSystematics)'
//e.g. root -b -l -q 'ZJetsFakeRate.C+("../../../EWKWZ-FF.root", 1)'

#include <vector>
#include "THStack.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "../Limits/consts.h"

const bool tex = false;

void
ZJetsFakeRate(string infile, bool useData=true, bool doSystematics=false){  
  bool doWLep = (infile.find("-Wiso") != string::npos); 
  cout<<"Do W Lep = "<<doWLep<<endl;

  //open file
  TFile *f = TFile::Open(infile.c_str(), "read"); assert(f);
  
  //setup samples
  vector<string> allsamples;
  if(!useData){
    allsamples.push_back("DYJetsToLL");
    if(!doSystematics){
      allsamples.push_back("WZJetsTo3LNu");//This one
      allsamples.push_back("GVJets");//Def not
      allsamples.push_back("WWTo2L2Nu");//Def not
      allsamples.push_back("ZZ");//Def not
      //allsamples.push_back("WJetsToLNu");
      allsamples.push_back("TTJets");
    }
  }else{
    allsamples.push_back("data");
  }
  


  TTree* t = NULL;
  //Do MC individually for pt agnostic numbers
  const int nchannel = 4;
  const int startchannel = 0;
  for(int ichannel=startchannel; ichannel<nchannel+startchannel; ++ichannel){
    bool useElectrons = ichannel%2 == false;
    cout<<"  ======== "<<(useElectrons ? "Electron" : "Muon")<<" ==========\n";
    string cuts = Form("weight*(ZTightCode==3)*(EvtType==%i)*(MET<20)*(WTransMass<20.)", ichannel);
    //cout<<"Cuts are : "<<cuts<<endl;

    if(!useData){
      for(int i=0; i<(int)allsamples.size(); ++i){
        //get tree from file
        t = getTree(f, allsamples[i], "tEvts_ValidW"); assert(t);
      
        float tot(0), pass(0);
        t->Draw("WTightCode", cuts.c_str(), "goff");
        int n = t->GetSelectedRows();
        for(int ientry=0; ientry<n; ++ientry){
          float weight = t->GetW()[ientry];
          int idx = 0;
          int WTightCode = round(t->GetVal(idx++)[ientry]);
          //assert(channel == ichannel || channel == (ichannel+2));
          tot += weight;
          if(WTightCode == true) pass += weight;
        }
        //calc fake rate
        float mean    = tot>0 ? pass / tot : 0.;//hFakeRateEff->GetY()       [bin-1];
        float errUp   = TEfficiency::ClopperPearson(tot, pass, 0.68, true) - mean;//hFakeRateEff->GetErrorYhigh(bin-1);
        float errDown = mean - TEfficiency::ClopperPearson(tot, pass, 0.68, false);//hFakeRateEff->GetErrorYlow(bin-1);
        errDown = max(errDown, (float) 0.);
        float err = (errUp + errDown)/2;
        printf(" %s &  %4.2f & %4.2f & %.2f \\pm %.2f\\%% \\\\\n", allsamples[i].c_str(), tot, pass, mean*100, err*100);
      }//sample loop
    }
    //get tree from file
    t = getTree(f, allsamples, "tEvts_ValidW"); assert(t);
    //Do sum for pt agnostic numbers
    float tot(0), pass(0);
    t->Draw("WTightCode", cuts.c_str(), "goff");
    int n = t->GetSelectedRows();
    for(int ientry=0; ientry<n; ++ientry){
      float weight = t->GetW()[ientry];
      int idx = 0;
      int WTightCode = round(t->GetVal(idx++)[ientry]);
      //assert(channel == ichannel || channel == (ichannel+2));
      tot += weight;
      if(WTightCode == true) pass += weight;
    }
    //calc fake rate
    float mean    = tot>0 ? pass / tot : 0.;//hFakeRateEff->GetY()       [bin-1];
    float errUp   = TEfficiency::ClopperPearson(tot, pass, 0.68, true) - mean;//hFakeRateEff->GetErrorYhigh(bin-1);
    float errDown = mean - TEfficiency::ClopperPearson(tot, pass, 0.68, false);//hFakeRateEff->GetErrorYlow(bin-1);
    errDown = max(errDown, (float) 0.);
    float err = (errUp + errDown)/2;
    printf(" %s & %4.2f & %4.2f & %.2f \\pm %.2f\\%% \\\\\n", useData ? "Data" : "Combined", tot, pass, mean*100, err*100);
  }//channel loop
  

  //Now do pt dep fake rates

  //loop over flavors
  for(int ichannel=startchannel; ichannel<nchannel+startchannel; ++ichannel){

    //loop over eta bins
    bool useElectrons = ichannel%2 == false;
    float* eta = NULL;
    int neta = 0;
    float etaElec[] = {0., 1.5, 2.5}; int netaElec = 3;
    //float etaMuon[] = {0., 1.0, 1.479, 2.0, 2.5}; int netaMuon = 5;
    float etaMuon[] = {0., 2.1, 2.5}; int netaMuon = 3;

    if(useElectrons){
      eta = etaElec; neta = netaElec;
    }else{
      eta = etaMuon; neta = netaMuon;
    }

    for(int ieta=0; ieta<neta-1; ++ieta){
      string sChannel = useElectrons ? "Electron" : "Muon";
      string sEta = ""; 
      if     (useElectrons && ieta == 0) sEta = "Barrel";
      else if(useElectrons && ieta == 1) sEta = "EndCap";
      else                               sEta = Form("%.2f < eta < %.2f", eta[ieta], eta[ieta+1]);
      printf("-------%s %s------\n",sChannel.c_str(), sEta.c_str());
      TCanvas c1;
      c1.Clear();
      c1.SetLogx(1);
      //c1.SetGrid();
      string title = useElectrons ? "Electron Fake Rate" : "Muon Fake Rate";
      title += useData ? " from Data" : " from MC";
      title += !useData && doSystematics ? " Systematics" : "";
      title += " using Z+Jets method";
      title += doWLep ? " on W lepton" : " on Z lepton";

      //loop over pt bins
      //float pt [] = {10., 15., 20., 25., 30., 100}; const int npt  = 6;//to match Alicia
      float pt [] = {10., 15., 40., 100}; const int npt  = 4;
      //float pt [] = {10., 15., 25., 35., 50., 1000}; const int npt  = 6;
      TGraphAsymmErrors* hFakeRatePt = new TGraphAsymmErrors(npt);
      hFakeRatePt->SetMaximum(1.1);
      hFakeRatePt->SetMinimum(0.);
      hFakeRatePt->SetTitle((title + ";p_{T} (GeV);Rate").c_str());
      for(int ipt=npt-2; ipt>=0; --ipt){

        float tot(0), pass(0);
        string cuts = Form("weight*(ZTightCode==3)*(EvtType==%i)*(MET<30)*(WTransMass<30.)*(NLeps<4)*(WLepPt >= %.0f && WLepPt < %.0f)*(abs(WLepEta) >= %.3f && abs(WLepEta) < %.3f)", 
                           ichannel, pt[ipt], pt[ipt+1], eta[ieta], eta[ieta+1]);
        //cout<<"Cuts are : "<<cuts<<endl;
        t->Draw("WTightCode", cuts.c_str(), "goff");
        int n = t->GetSelectedRows();
        for(int ientry=0; ientry<n; ++ientry){
          float weight = t->GetW()[ientry];
          int idx = 0;
          int WTightCode = round(t->GetVal(idx++)[ientry]);
          //assert(channel == ichannel || channel == (ichannel+2));
          tot += weight;
          if(WTightCode == true) pass += weight;

        }
        //calc fake rate
        float mean    = tot>0 ? pass / tot : 0.;//hFakeRateEff->GetY()       [bin-1];
        float errUp   = TEfficiency::ClopperPearson(tot, pass, 0.68, true) - mean;//hFakeRateEff->GetErrorYhigh(bin-1);
        float errDown = mean - TEfficiency::ClopperPearson(tot, pass, 0.68, false);//hFakeRateEff->GetErrorYlow(bin-1);
        errDown = max(errDown, (float) 0.);
        float err = (errUp + errDown)/2;
        
        if(!tex) printf("if(pt > %.0f) return Value(%.4f, %.4f); //channel %i, pt %3.0f eta %.1f, (mean, err)= (%.4f, %.4f) (%.2f / %.2f) +%.4f -%.4f\n", 
                        pt[ipt], mean, err, ichannel, pt[ipt], eta[ieta],  mean, err, pass, tot, errUp, errDown);
        else printf("& $%.3f_{%.3f}^{%.3f} ", mean, errDown, errUp);

        float x = (pt[ipt] + pt[ipt+1])/2.;
        hFakeRatePt->SetPoint      (ipt, x, mean);
        hFakeRatePt->SetPointEYhigh(ipt, errUp);
        hFakeRatePt->SetPointEYlow (ipt, errDown);
        hFakeRatePt->SetPointEXhigh(ipt, pt[ipt+1]- x);
        hFakeRatePt->SetPointEXlow (ipt, x - pt[ipt]);
  
      }//pt loop
      hFakeRatePt->Draw("ap*");
      string outName = Form("ZJetsFakeRatePt-Eta%.1fto%.1f-", eta[ieta], eta[ieta+1]);
      replace(outName.begin(), outName.end(), '.', 'p');
      outName += useElectrons ? "Elec" : "Muon";
      outName += useData ? "Data" : "MC";
      outName += !useData && doSystematics ? "Sys" : "";
      outName += doWLep ? "WLep" : "ZLep";
      outName += ".pdf";
      c1.Print(outName.c_str());
    }//eta loop
  }//channel loop
}

