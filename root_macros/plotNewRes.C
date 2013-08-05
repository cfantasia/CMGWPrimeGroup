#include "UserCode/CMGWPrimeGroup/root_macros/common.h"
#include "TPaveStats.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include <fstream>

void
plotNewRes(){
  TH1::StatOverflows(kTRUE);

  TH1F* hWZRes = new TH1F("hWZRes", "Res", 200, -500, 500);


  TTree* tEvts = new TTree("tLimit", "Limits");
  tEvts->ReadFile("fixRes.txt");
  
  string cuts = Form("(EvtType==2)");
  int n = tEvts->Draw("EvtType:ZLep1PtGen:ZLep1EtaGen:ZLep1PhiGen:ZLep2PtGen:ZLep2EtaGen:ZLep2PhiGen:bilZLep1Pt:ZLep1PtErr:ZLep1Eta:ZLep1Phi:bilZLep2Pt:ZLep2PtErr:ZLep2Eta:ZLep2Phi:WLepPtGen:WLepEtaGen:WLepPhiGen:WNeuPtGen:WNeuEtaGen:WNeuPhiGen:WLepPt:WLepPtErr:WLepEta:WLepPhi:WNeuPt:WNeuPtErr:WNeuEta:WNeuPhi", 
                      cuts.c_str(), "para goff");
  printf("There are %i events total\n", n);
  for(int ievt=0; ievt<n; ++ievt){  //loop over events
    printf("--------------------\n");
    int treeIdx = 0;
    const int EvtType  = round(tEvts->GetVal(treeIdx++)[ievt]);

    const double ZLep1PtGen  = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep1EtaGen = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep1PhiGen = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2PtGen  = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2EtaGen = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2PhiGen = tEvts->GetVal(treeIdx++)[ievt];
    
    const double ZLep1Pt  = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep1PtErr= tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep1Eta = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep1Phi = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2Pt  = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2PtErr= tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2Eta = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2Phi = tEvts->GetVal(treeIdx++)[ievt];

    const double WLepPtGen  = tEvts->GetVal(treeIdx++)[ievt];
    const double WLepEtaGen = tEvts->GetVal(treeIdx++)[ievt];
    const double WLepPhiGen = tEvts->GetVal(treeIdx++)[ievt];
    const double WNeuPtGen  = tEvts->GetVal(treeIdx++)[ievt];
    const double WNeuEtaGen = tEvts->GetVal(treeIdx++)[ievt];
    const double WNeuPhiGen = tEvts->GetVal(treeIdx++)[ievt];
    
    const double WLepPt  = tEvts->GetVal(treeIdx++)[ievt];
    const double WLepPtErr= tEvts->GetVal(treeIdx++)[ievt];
    const double WLepEta = tEvts->GetVal(treeIdx++)[ievt];
    const double WLepPhi = tEvts->GetVal(treeIdx++)[ievt];
    const double WNeuPt  = tEvts->GetVal(treeIdx++)[ievt];
    const double WNeuPtErr= tEvts->GetVal(treeIdx++)[ievt];
    const double WNeuEta = tEvts->GetVal(treeIdx++)[ievt];
    const double WNeuPhi = tEvts->GetVal(treeIdx++)[ievt];

    const float weight       = tEvts->GetW()[ievt];
    
    //Gen
    TLorentzVector l1_gen, l2_gen, zsum_gen;
    if(ZLep1PtGen > ZLep2PtGen){
      l1_gen.SetPtEtaPhiM(   ZLep1PtGen,ZLep1EtaGen,ZLep1PhiGen,0.); 
      l2_gen.SetPtEtaPhiM(   ZLep2PtGen,ZLep2EtaGen,ZLep2PhiGen,0.); 
    }else{
      l2_gen.SetPtEtaPhiM(   ZLep1PtGen,ZLep1EtaGen,ZLep1PhiGen,0.); 
      l1_gen.SetPtEtaPhiM(   ZLep2PtGen,ZLep2EtaGen,ZLep2PhiGen,0.); 
    }
    TLorentzVector l3_gen, l4_gen, wsum_gen;
    l3_gen.SetPtEtaPhiM(   WLepPtGen,WLepEtaGen,WLepPhiGen,0.); 
    l4_gen.SetPtEtaPhiM(   WNeuPtGen,WNeuEtaGen,WNeuPhiGen,0.); 

    TLorentzVector sum_gen = l1_gen + l2_gen + l3_gen + l4_gen;

    //Old
    TLorentzVector l1_old, l2_old, zsum_old;
    l1_old.SetPtEtaPhiM(   ZLep1Pt,ZLep1Eta,ZLep1Phi,0.); 
    l2_old.SetPtEtaPhiM(   ZLep2Pt,ZLep2Eta,ZLep2Phi,0.); 

    TLorentzVector l3_old, l4_old, wsum_old;
    l3_old.SetPtEtaPhiM(   WLepPt,WLepEta,WLepPhi,0.); 
    l4_old.SetPtEtaPhiM(   WNeuPt,WNeuEta,WNeuPhi,0.); 

    TLorentzVector sum_old = l1_old + l2_old + l3_old + l4_old;

    printf("Gen: %.1f New: %.1f Diff: %.1f\n", sum_gen.M(), sum_old.M(), sum_gen.M() - sum_old.M());

    hWZRes->Fill(sum_gen.M() - sum_old.M());
  }
  
  hWZRes->Draw();
  printf("Resolution is %.2f\n", hWZRes->GetRMS());
}
