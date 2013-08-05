//Author: Cory Fantasia 2013
//Purpose: Plot resolution of lepton pt, met, wz mass
//Usage: root -b -l -q 'PlotResolution.C+("../../../WprimeWZ.root", 2000)'

#include "UserCode/CMGWPrimeGroup/root_macros/common.h"
#include "TPaveStats.h"
#include "TCanvas.h"
#include "TH1F.h"

void shiftStats(TPaveStats* stat1, TPaveStats* stat2);
float findMax(TH1F* a, TH1F* b);
float findMax(TH1F* a, TH1F* b, TH1F* c, TH1F* d);

void
PlotResolution(string inName="", int mass=2000){
  TH1::StatOverflows(kTRUE);
  string sample = Form("WprimeToWZTo3LNu_M-%i", mass);
  string outName = Form("WprimeWZ%i", mass);

  if(inName.empty()) inName = "../../../WprimeWZ.root";
  TFile *_file0 = TFile::Open(inName.c_str()); assert(_file0);
  TPaveStats *stat1,*stat2,*stat3,*stat4; 

  TCanvas* cLepRes = new TCanvas("Res");

  TH1F* WLepPt_muon = new TH1F("WLepPt_muon", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WFlavorGen && ZFlavorGen && EvtType%2 == 1", *WLepPt_muon);
  WLepPt_muon->SetMarkerColor(kBlue); WLepPt_muon->SetLineColor(kBlue); WLepPt_muon->Draw(""); cLepRes->Update(); 
  stat3 = (TPaveStats*) (WLepPt_muon->GetListOfFunctions()->FindObject("stats")); stat3->SetTextColor(kBlue);
  //shiftStats(stat2, stat3);
  TH1F* WLepPt_elec = new TH1F("WLepPt_elec", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WFlavorGen && ZFlavorGen && EvtType%2 == 0", *WLepPt_elec);
  WLepPt_elec->SetMarkerColor(kRed); WLepPt_elec->SetLineColor(kRed); WLepPt_elec->Draw("sames"); cLepRes->Update(); 
  stat4 = (TPaveStats*) (WLepPt_elec->GetListOfFunctions()->FindObject("stats")); stat4->SetTextColor(kRed);
  shiftStats(stat3, stat4);

  WLepPt_muon->SetMaximum(findMax(WLepPt_muon, WLepPt_elec));
  cLepRes->SaveAs((outName + "_WLeptonPtResolution.pdf").c_str());
  ////////////////////////////////////

  TH1F* MET_mmm = new TH1F("MET_mmm", "MET; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "MET-WNeuPtGen", "WFlavorGen && ZFlavorGen && EvtType == 3", *MET_mmm);
  MET_mmm->SetMarkerColor(kOrange+2  ); MET_mmm->SetLineColor(kOrange+2  ); MET_mmm->Draw(); cLepRes->Update(); 
  stat1 = (TPaveStats*) (MET_mmm->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kOrange+2);
  TH1F* MET_emm = new TH1F("MET_emm", "MET; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "MET-WNeuPtGen", "WFlavorGen && ZFlavorGen && EvtType == 2", *MET_emm);
  MET_emm->SetMarkerColor(kGreen); MET_emm->SetLineColor(kGreen); MET_emm->Draw("sames"); cLepRes->Update(); 
  stat2 = (TPaveStats*) (MET_emm->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kGreen);
  shiftStats(stat1, stat2);
  TH1F* MET_eem = new TH1F("MET_eem", "MET; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "MET-WNeuPtGen", "WFlavorGen && ZFlavorGen && EvtType == 1", *MET_eem);
  MET_eem->SetMarkerColor(kBlue); MET_eem->SetLineColor(kBlue); MET_eem->Draw("sames"); cLepRes->Update(); 
  stat3 = (TPaveStats*) (MET_eem->GetListOfFunctions()->FindObject("stats")); stat3->SetTextColor(kBlue);
  shiftStats(stat2, stat3);
  TH1F* MET_eee = new TH1F("MET_eee", "MET; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "MET-WNeuPtGen", "WFlavorGen && ZFlavorGen && EvtType == 0", *MET_eee);
  MET_eee->SetMarkerColor(kRed); MET_eee->SetLineColor(kRed); MET_eee->Draw("sames"); cLepRes->Update(); 
  stat4 = (TPaveStats*) (MET_eee->GetListOfFunctions()->FindObject("stats")); stat4->SetTextColor(kRed);
  shiftStats(stat3, stat4);

  MET_mmm->SetMaximum(findMax(MET_mmm, MET_emm, MET_eem, MET_eee));
  cLepRes->SaveAs((outName + "_METResolution.pdf").c_str());
  ////////////////////////////////////

  TH1F* ZLep1Pt_muon = new TH1F("ZLep1Pt_muon", "Z Leading Lepton; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "ZLep1Pt-max(ZLep2PtGen,ZLep1PtGen)", "WFlavorGen && ZFlavorGen && EvtType >= 2", *ZLep1Pt_muon);
  ZLep1Pt_muon->SetMarkerColor(kBlue  ); ZLep1Pt_muon->SetLineColor(kBlue  ); ZLep1Pt_muon->Draw(); cLepRes->Update(); 
  stat1 = (TPaveStats*) (ZLep1Pt_muon->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kBlue);
  TH1F* ZLep1Pt_elec = new TH1F("ZLep1Pt_elec", "Z Leading Lepton; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "ZLep1Pt-max(ZLep2PtGen,ZLep1PtGen)", "WFlavorGen && ZFlavorGen && EvtType <= 1", *ZLep1Pt_elec);
  ZLep1Pt_elec->SetMarkerColor(kRed); ZLep1Pt_elec->SetLineColor(kRed); ZLep1Pt_elec->Draw("sames"); cLepRes->Update(); 
  stat2 = (TPaveStats*) (ZLep1Pt_elec->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kRed);
  shiftStats(stat1, stat2);

  ZLep1Pt_muon->SetMaximum(findMax(ZLep1Pt_muon, ZLep1Pt_elec));
  cLepRes->SaveAs((outName + "_ZLeadingLeptonPtResolution.pdf").c_str());
  ////////////////////////////////////

  TH1F* ZLep2Pt_muon = new TH1F("ZLep2Pt_muon", "Z Trailing Lepton; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "ZLep2Pt-min(ZLep2PtGen,ZLep1PtGen)", "WFlavorGen && ZFlavorGen && EvtType >= 2", *ZLep2Pt_muon);
  ZLep2Pt_muon->SetMarkerColor(kBlue  ); ZLep2Pt_muon->SetLineColor(kBlue  ); ZLep2Pt_muon->Draw(); cLepRes->Update(); 
  stat1 = (TPaveStats*) (ZLep2Pt_muon->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kBlue);
  TH1F* ZLep2Pt_elec = new TH1F("ZLep2Pt_elec", "Z Trailing Lepton; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "ZLep2Pt-min(ZLep2PtGen,ZLep1PtGen)", "WFlavorGen && ZFlavorGen && EvtType <= 1", *ZLep2Pt_elec);
  ZLep2Pt_elec->SetMarkerColor(kRed); ZLep2Pt_elec->SetLineColor(kRed); ZLep2Pt_elec->Draw("sames"); cLepRes->Update(); 
  stat2 = (TPaveStats*) (ZLep2Pt_elec->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kRed);
  shiftStats(stat1, stat2);

  ZLep2Pt_muon->SetMaximum(findMax(ZLep2Pt_muon, ZLep2Pt_elec));
  cLepRes->SaveAs((outName + "_ZTrailingLeptonPtResolution.pdf").c_str());
  ////////////////////////////////////

  TCanvas* cWZRes = new TCanvas("WZ Resolution");

  TH1F* WZRes_mmm = new TH1F("WZRes_mmm", "WZ Mass; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "WZMass-WprimeGenMass", "WFlavorGen && ZFlavorGen && EvtType == 3", *WZRes_mmm);
  WZRes_mmm->SetMarkerColor(kOrange+2  ); WZRes_mmm->SetLineColor(kOrange+2  ); WZRes_mmm->Draw(); cWZRes->Update(); 
  stat1 = (TPaveStats*) (WZRes_mmm->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kOrange+2);
  TH1F* WZRes_emm = new TH1F("WZRes_emm", "WZ Mass; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "WZMass-WprimeGenMass", "WFlavorGen && ZFlavorGen && EvtType == 2", *WZRes_emm);
  WZRes_emm->SetMarkerColor(kGreen); WZRes_emm->SetLineColor(kGreen); WZRes_emm->Draw("sames"); cWZRes->Update(); 
  stat2 = (TPaveStats*) (WZRes_emm->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kGreen);
  shiftStats(stat1, stat2);
  TH1F* WZRes_eem = new TH1F("WZRes_eem", "WZ Mass; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "WZMass-WprimeGenMass", "WFlavorGen && ZFlavorGen && EvtType == 1", *WZRes_eem);
  WZRes_eem->SetMarkerColor(kBlue); WZRes_eem->SetLineColor(kBlue); WZRes_eem->Draw("sames"); cWZRes->Update(); 
  stat3 = (TPaveStats*) (WZRes_eem->GetListOfFunctions()->FindObject("stats")); stat3->SetTextColor(kBlue);
  shiftStats(stat2, stat3);
  TH1F* WZRes_eee = new TH1F("WZRes_eee", "WZ Mass; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "WZMass-WprimeGenMass", "WFlavorGen && ZFlavorGen && EvtType == 0", *WZRes_eee);
  WZRes_eee->SetMarkerColor(kRed); WZRes_eee->SetLineColor(kRed); WZRes_eee->Draw("sames"); cWZRes->Update(); 
  stat4 = (TPaveStats*) (WZRes_eee->GetListOfFunctions()->FindObject("stats")); stat4->SetTextColor(kRed);
  shiftStats(stat3, stat4);

  WZRes_mmm->SetMaximum(findMax(WZRes_mmm, WZRes_emm, WZRes_eem, WZRes_eee));
  cWZRes->SaveAs((outName + "_WZResolution.pdf").c_str());
  printf(" mass %i maxRes %.0f\n", mass, max(
           max(WZRes_mmm->GetRMS(), WZRes_emm->GetRMS()), 
           max(WZRes_eem->GetRMS(), WZRes_eee->GetRMS())
           ));
  /////////////////////////////////////
  
  TCanvas* cResPt = new TCanvas("Resolution as Function of Pt");

  TH1F* WLepPt0_100_muon = new TH1F("WLepPt0_100_muon", "0 < W Lepton Pt < 100; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 0 && WLepPtGen < 100 && WFlavorGen && ZFlavorGen && EvtType%2 == 1", *WLepPt0_100_muon);
  WLepPt0_100_muon->SetMarkerColor(kBlue  ); WLepPt0_100_muon->SetLineColor(kBlue  ); WLepPt0_100_muon->Draw(); cResPt->Update(); 
  stat1 = (TPaveStats*) (WLepPt0_100_muon->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kBlue);
  TH1F* WLepPt0_100_elec = new TH1F("WLepPt0_100_elec", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 0 && WLepPtGen < 100 && WFlavorGen && ZFlavorGen && EvtType%2 == 0", *WLepPt0_100_elec);
  WLepPt0_100_elec->SetMarkerColor(kRed); WLepPt0_100_elec->SetLineColor(kRed); WLepPt0_100_elec->Draw("sames"); cResPt->Update(); 
  stat2 = (TPaveStats*) (WLepPt0_100_elec->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kRed);
  shiftStats(stat1, stat2); 

  WLepPt0_100_muon->SetMaximum(findMax(WLepPt0_100_muon, WLepPt0_100_elec));
  cResPt->SaveAs((outName + "_WLeptonPtResolution_0-100.pdf").c_str());
  /////////////////////////////////////

  TH1F* WLepPt100_200_muon = new TH1F("WLepPt100_200_muon", "100 < W Lepton Pt < 200; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 100 && WLepPtGen < 200 && WFlavorGen && ZFlavorGen && EvtType%2 == 1", *WLepPt100_200_muon);
  WLepPt100_200_muon->SetMarkerColor(kBlue  ); WLepPt100_200_muon->SetLineColor(kBlue  ); WLepPt100_200_muon->Draw(); cResPt->Update(); 
  stat1 = (TPaveStats*) (WLepPt100_200_muon->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kBlue);
  TH1F* WLepPt100_200_elec = new TH1F("WLepPt100_200_elec", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 100 && WLepPtGen < 200 && WFlavorGen && ZFlavorGen && EvtType%2 == 0", *WLepPt100_200_elec);
  WLepPt100_200_elec->SetMarkerColor(kRed); WLepPt100_200_elec->SetLineColor(kRed); WLepPt100_200_elec->Draw("sames"); cResPt->Update(); 
  stat2 = (TPaveStats*) (WLepPt100_200_elec->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kRed);
  shiftStats(stat1, stat2); 

  WLepPt100_200_muon->SetMaximum(findMax(WLepPt100_200_muon, WLepPt100_200_elec));
  cResPt->SaveAs((outName + "_WLeptonPtResolution_100-200.pdf").c_str());
  /////////////////////////////////////

  TH1F* WLepPt200_300_muon = new TH1F("WLepPt200_300_muon", "200 < W Lepton Pt < 300; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 200 && WLepPtGen < 300 && WFlavorGen && ZFlavorGen && EvtType%2 == 1", *WLepPt200_300_muon);
  WLepPt200_300_muon->SetMarkerColor(kBlue  ); WLepPt200_300_muon->SetLineColor(kBlue  ); WLepPt200_300_muon->Draw(); cResPt->Update(); 
  stat1 = (TPaveStats*) (WLepPt200_300_muon->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kBlue);
  TH1F* WLepPt200_300_elec = new TH1F("WLepPt200_300_elec", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 200 && WLepPtGen < 300 && WFlavorGen && ZFlavorGen && EvtType%2 == 0", *WLepPt200_300_elec);
  WLepPt200_300_elec->SetMarkerColor(kRed); WLepPt200_300_elec->SetLineColor(kRed); WLepPt200_300_elec->Draw("sames"); cResPt->Update(); 
  stat2 = (TPaveStats*) (WLepPt200_300_elec->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kRed);
  shiftStats(stat1, stat2); 

  WLepPt200_300_muon->SetMaximum(findMax(WLepPt200_300_muon, WLepPt200_300_elec));
  cResPt->SaveAs((outName + "_WLeptonPtResolution_200-300.pdf").c_str());
  /////////////////////////////////////

  TH1F* WLepPt300_400_muon = new TH1F("WLepPt300_400_muon", "300 < W Lepton Pt < 400; Resolution (Reco - Gen) (GeV); a.u.", 300, -300, 300); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 300 && WLepPtGen < 400 && WFlavorGen && ZFlavorGen && EvtType%2 == 1", *WLepPt300_400_muon);
  WLepPt300_400_muon->SetMarkerColor(kBlue  ); WLepPt300_400_muon->SetLineColor(kBlue  ); WLepPt300_400_muon->Draw(); cResPt->Update(); 
  stat1 = (TPaveStats*) (WLepPt300_400_muon->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kBlue);
  TH1F* WLepPt300_400_elec = new TH1F("WLepPt300_400_elec", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 300, -300, 300); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 300 && WLepPtGen < 400 && WFlavorGen && ZFlavorGen && EvtType%2 == 0", *WLepPt300_400_elec);
  WLepPt300_400_elec->SetMarkerColor(kRed); WLepPt300_400_elec->SetLineColor(kRed); WLepPt300_400_elec->Draw("sames"); cResPt->Update(); 
  stat2 = (TPaveStats*) (WLepPt300_400_elec->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kRed);
  shiftStats(stat1, stat2); 

  WLepPt300_400_muon->SetMaximum(findMax(WLepPt300_400_muon, WLepPt300_400_elec));
  cResPt->SaveAs((outName + "_WLeptonPtResolution_300-400.pdf").c_str());
  /////////////////////////////////////

  TH1F* WLepPt400_500_muon = new TH1F("WLepPt400_500_muon", "400 < W Lepton < 500; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 400 && WLepPtGen < 500 && WFlavorGen && ZFlavorGen && EvtType%2 == 1", *WLepPt400_500_muon);
  WLepPt400_500_muon->SetMarkerColor(kBlue  ); WLepPt400_500_muon->SetLineColor(kBlue  ); WLepPt400_500_muon->Draw(); cResPt->Update(); 
  stat1 = (TPaveStats*) (WLepPt400_500_muon->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kBlue);
  TH1F* WLepPt400_500_elec = new TH1F("WLepPt400_500_elec", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 400 && WLepPtGen < 500 && WFlavorGen && ZFlavorGen && EvtType%2 == 0", *WLepPt400_500_elec);
  WLepPt400_500_elec->SetMarkerColor(kRed); WLepPt400_500_elec->SetLineColor(kRed); WLepPt400_500_elec->Draw("sames"); cResPt->Update(); 
  stat2 = (TPaveStats*) (WLepPt400_500_elec->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kRed);
  shiftStats(stat1, stat2); 

  WLepPt400_500_muon->SetMaximum(findMax(WLepPt400_500_muon, WLepPt400_500_elec));
  cResPt->SaveAs((outName + "_WLeptonPtResolution_400-500.pdf").c_str());
  /////////////////////////////////////

  TH1F* WLepPt500_600_muon = new TH1F("WLepPt500_600_muon", "500 < W Lepton < 600; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 500 && WLepPtGen < 600 && WFlavorGen && ZFlavorGen && EvtType%2 == 1", *WLepPt500_600_muon);
  WLepPt500_600_muon->SetMarkerColor(kBlue  ); WLepPt500_600_muon->SetLineColor(kBlue  ); WLepPt500_600_muon->Draw(); cResPt->Update(); 
  stat1 = (TPaveStats*) (WLepPt500_600_muon->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kBlue);
  TH1F* WLepPt500_600_elec = new TH1F("WLepPt500_600_elec", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 500 && WLepPtGen < 600 && WFlavorGen && ZFlavorGen && EvtType%2 == 0", *WLepPt500_600_elec);
  WLepPt500_600_elec->SetMarkerColor(kRed); WLepPt500_600_elec->SetLineColor(kRed); WLepPt500_600_elec->Draw("sames"); cResPt->Update(); 
  stat2 = (TPaveStats*) (WLepPt500_600_elec->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kRed);
  shiftStats(stat1, stat2); 

  WLepPt500_600_muon->SetMaximum(findMax(WLepPt500_600_muon, WLepPt500_600_elec));
  cResPt->SaveAs((outName + "_WLeptonPtResolution_500-600.pdf").c_str());
  /////////////////////////////////////

  TH1F* WLepPt600_800_muon = new TH1F("WLepPt600_800_muon", "600 < W Lepton < 800; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 600 && WLepPtGen < 800 && WFlavorGen && ZFlavorGen && EvtType%2 == 1", *WLepPt600_800_muon);
  WLepPt600_800_muon->SetMarkerColor(kBlue  ); WLepPt600_800_muon->SetLineColor(kBlue  ); WLepPt600_800_muon->Draw(); cResPt->Update(); 
  stat1 = (TPaveStats*) (WLepPt600_800_muon->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kBlue);
  TH1F* WLepPt600_800_elec = new TH1F("WLepPt600_800_elec", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 600 && WLepPtGen < 800 && WFlavorGen && ZFlavorGen && EvtType%2 == 0", *WLepPt600_800_elec);
  WLepPt600_800_elec->SetMarkerColor(kRed); WLepPt600_800_elec->SetLineColor(kRed); WLepPt600_800_elec->Draw("sames"); cResPt->Update(); 
  stat2 = (TPaveStats*) (WLepPt600_800_elec->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kRed);
  shiftStats(stat1, stat2); 

  WLepPt600_800_muon->SetMaximum(findMax(WLepPt600_800_muon, WLepPt600_800_elec));
  cResPt->SaveAs((outName + "_WLeptonPtResolution_600-800.pdf").c_str());
  /////////////////////////////////////

  TH1F* WLepPt800_5000_muon = new TH1F("WLepPt800_5000_muon", "800 < W Lepton < 5000; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 800 && WLepPtGen < 5000 && WFlavorGen && ZFlavorGen && EvtType%2 == 1", *WLepPt800_5000_muon);
  WLepPt800_5000_muon->SetMarkerColor(kBlue  ); WLepPt800_5000_muon->SetLineColor(kBlue  ); WLepPt800_5000_muon->Draw(); cResPt->Update(); 
  stat1 = (TPaveStats*) (WLepPt800_5000_muon->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kBlue);
  TH1F* WLepPt800_5000_elec = new TH1F("WLepPt800_5000_elec", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 800 && WLepPtGen < 5000 && WFlavorGen && ZFlavorGen && EvtType%2 == 0", *WLepPt800_5000_elec);
  WLepPt800_5000_elec->SetMarkerColor(kRed); WLepPt800_5000_elec->SetLineColor(kRed); WLepPt800_5000_elec->Draw("sames"); cResPt->Update(); 
  stat2 = (TPaveStats*) (WLepPt800_5000_elec->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kRed);
  shiftStats(stat1, stat2); 

  WLepPt800_5000_muon->SetMaximum(findMax(WLepPt800_5000_muon, WLepPt800_5000_elec));
  cResPt->SaveAs((outName + "_WLeptonPtResolution_800-5000.pdf").c_str());
  /////////////////////////////////////

  TCanvas* cMETScaledRes = new TCanvas("Scaled MET Resolution");

  ///Scaled MET Resolution
  TH1F* METScaled_mmm = new TH1F("METScaled_mmm", "METScaled; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "min(MET, MET*(80.398/WTransMass)^2) -WNeuPtGen", "WFlavorGen && ZFlavorGen && EvtType == 3", *METScaled_mmm);
  METScaled_mmm->SetMarkerColor(kOrange+2  ); METScaled_mmm->SetLineColor(kOrange+2  ); METScaled_mmm->Draw(); cMETScaledRes->Update(); 
  stat1 = (TPaveStats*) (METScaled_mmm->GetListOfFunctions()->FindObject("stats")); assert( stat1); stat1->SetTextColor(kOrange+2);
  TH1F* METScaled_emm = new TH1F("METScaled_emm", "METScaled; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "min(MET, MET*(80.398/WTransMass)^2) -WNeuPtGen", "WFlavorGen && ZFlavorGen && EvtType == 2", *METScaled_emm);
  METScaled_emm->SetMarkerColor(kGreen); METScaled_emm->SetLineColor(kGreen); METScaled_emm->Draw("sames"); cMETScaledRes->Update(); 
  stat2 = (TPaveStats*) (METScaled_emm->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kGreen);
  shiftStats(stat1, stat2);
  TH1F* METScaled_eem = new TH1F("METScaled_eem", "METScaled; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "min(MET, MET*(80.398/WTransMass)^2) -WNeuPtGen", "WFlavorGen && ZFlavorGen && EvtType == 1", *METScaled_eem);
  METScaled_eem->SetMarkerColor(kBlue); METScaled_eem->SetLineColor(kBlue); METScaled_eem->Draw("sames"); cMETScaledRes->Update(); 
  stat3 = (TPaveStats*) (METScaled_eem->GetListOfFunctions()->FindObject("stats")); stat3->SetTextColor(kBlue);
  shiftStats(stat2, stat3);
  TH1F* METScaled_eee = new TH1F("METScaled_eee", "METScaled; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "min(MET, MET*(80.398/WTransMass)^2) -WNeuPtGen", "WFlavorGen && ZFlavorGen && EvtType == 0", *METScaled_eee);
  METScaled_eee->SetMarkerColor(kRed); METScaled_eee->SetLineColor(kRed); METScaled_eee->Draw("sames"); cMETScaledRes->Update(); 
  stat4 = (TPaveStats*) (METScaled_eee->GetListOfFunctions()->FindObject("stats")); stat4->SetTextColor(kRed);
  shiftStats(stat3, stat4);
  METScaled_mmm->SetMaximum(findMax(METScaled_mmm, METScaled_emm, METScaled_eem, METScaled_eee));
  cMETScaledRes->SaveAs((outName + "_METScaledResolution.pdf").c_str());
  ////////////////////////////////////
}


void 
shiftStats(TPaveStats* stat1, TPaveStats* stat2){
  float height = stat1->GetY2NDC() - stat1->GetY1NDC();
  stat2->SetY2NDC(stat1->GetY1NDC() );
  stat2->SetY1NDC(stat2->GetY2NDC() - height);
  stat2->Draw();
}

float
findMax(TH1F* a, TH1F* b){
  Double_t max = 0;
  max = TMath::Max(max, a->GetMaximum());
  max = TMath::Max(max, b->GetMaximum());
  return max*1.1;
}

float
findMax(TH1F* a, TH1F* b, TH1F* c, TH1F* d){
  Double_t max = 0;
  max = TMath::Max(max, a->GetMaximum());
  max = TMath::Max(max, b->GetMaximum());
  max = TMath::Max(max, c->GetMaximum());
  max = TMath::Max(max, d->GetMaximum());
  return max*1.1;
}
