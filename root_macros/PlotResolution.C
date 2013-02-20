#include "UserCode/CMGWPrimeGroup/root_macros/common.h"
#include "TPaveStats.h"
#include "TCanvas.h"
#include "TH1F.h"

void shiftStats(TPaveStats* stat1, TPaveStats* stat2);
float findMax(TH1F* a, TH1F* b, TH1F* c, TH1F* d);

void
PlotRes(int mass){
  TH1::StatOverflows(kTRUE);
  string sample = Form("WprimeToWZTo3LNu_M-%i", mass);
  string outName = Form("WprimeWZ%i", mass);

  TFile *_file0 = TFile::Open("WprimeWZ.root");
  TPaveStats *stat1,*stat2,*stat3,*stat4; 

  TCanvas* cLepRes = new TCanvas("Res");

  TH1F* WLepPt_mmm = new TH1F("WLepPt_mmm", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "EvtType == 3", *WLepPt_mmm);
  WLepPt_mmm->SetMarkerColor(kYellow  ); WLepPt_mmm->SetLineColor(kYellow  ); WLepPt_mmm->Draw(); cLepRes->Update(); 
  stat1 = (TPaveStats*) (WLepPt_mmm->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kYellow);
  TH1F* WLepPt_emm = new TH1F("WLepPt_emm", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "EvtType == 2", *WLepPt_emm);
  WLepPt_emm->SetMarkerColor(kGreen); WLepPt_emm->SetLineColor(kGreen); WLepPt_emm->Draw("sames"); cLepRes->Update(); 
  stat2 = (TPaveStats*) (WLepPt_emm->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kGreen);
  shiftStats(stat1, stat2); 
  TH1F* WLepPt_eem = new TH1F("WLepPt_eem", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "EvtType == 1", *WLepPt_eem);
  WLepPt_eem->SetMarkerColor(kBlue); WLepPt_eem->SetLineColor(kBlue); WLepPt_eem->Draw("sames"); cLepRes->Update(); 
  stat3 = (TPaveStats*) (WLepPt_eem->GetListOfFunctions()->FindObject("stats")); stat3->SetTextColor(kBlue);
  shiftStats(stat2, stat3);
  TH1F* WLepPt_eee = new TH1F("WLepPt_eee", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "EvtType == 0", *WLepPt_eee);
  WLepPt_eee->SetMarkerColor(kRed); WLepPt_eee->SetLineColor(kRed); WLepPt_eee->Draw("sames"); cLepRes->Update(); 
  stat4 = (TPaveStats*) (WLepPt_eee->GetListOfFunctions()->FindObject("stats")); stat4->SetTextColor(kRed);
  shiftStats(stat3, stat4);
  WLepPt_mmm->SetMaximum(findMax(WLepPt_mmm, WLepPt_emm, WLepPt_eem, WLepPt_eee));

  cLepRes->SaveAs((outName + "_WLeptonPtResolution.png").c_str());
  ////////////////////////////////////

  cLepRes->cd(2);
  TH1F* MET_mmm = new TH1F("MET_mmm", "MET; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "MET-WNeuPtGen", "EvtType == 3", *MET_mmm);
  MET_mmm->SetMarkerColor(kYellow  ); MET_mmm->SetLineColor(kYellow  ); MET_mmm->Draw(); cLepRes->Update(); 
  stat1 = (TPaveStats*) (MET_mmm->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kYellow);
  TH1F* MET_emm = new TH1F("MET_emm", "MET; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "MET-WNeuPtGen", "EvtType == 2", *MET_emm);
  MET_emm->SetMarkerColor(kGreen); MET_emm->SetLineColor(kGreen); MET_emm->Draw("sames"); cLepRes->Update(); 
  stat2 = (TPaveStats*) (MET_emm->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kGreen);
  shiftStats(stat1, stat2);
  TH1F* MET_eem = new TH1F("MET_eem", "MET; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "MET-WNeuPtGen", "EvtType == 1", *MET_eem);
  MET_eem->SetMarkerColor(kBlue); MET_eem->SetLineColor(kBlue); MET_eem->Draw("sames"); cLepRes->Update(); 
  stat3 = (TPaveStats*) (MET_eem->GetListOfFunctions()->FindObject("stats")); stat3->SetTextColor(kBlue);
  shiftStats(stat2, stat3);
  TH1F* MET_eee = new TH1F("MET_eee", "MET; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "MET-WNeuPtGen", "EvtType == 0", *MET_eee);
  MET_eee->SetMarkerColor(kRed); MET_eee->SetLineColor(kRed); MET_eee->Draw("sames"); cLepRes->Update(); 
  stat4 = (TPaveStats*) (MET_eee->GetListOfFunctions()->FindObject("stats")); stat4->SetTextColor(kRed);
  shiftStats(stat3, stat4);
  MET_mmm->SetMaximum(findMax(MET_mmm, MET_emm, MET_eem, MET_eee));

  cLepRes->SaveAs((outName + "_METResolution.png").c_str());
  ////////////////////////////////////

  cLepRes->cd(3);
  TH1F* ZLep1Pt_mmm = new TH1F("ZLep1Pt_mmm", "Z Leading Lepton; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "ZLep1Pt-max(ZLep2PtGen,ZLep1PtGen)", "EvtType == 3", *ZLep1Pt_mmm);
  ZLep1Pt_mmm->SetMarkerColor(kYellow  ); ZLep1Pt_mmm->SetLineColor(kYellow  ); ZLep1Pt_mmm->Draw(); cLepRes->Update(); 
  stat1 = (TPaveStats*) (ZLep1Pt_mmm->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kYellow);
  TH1F* ZLep1Pt_emm = new TH1F("ZLep1Pt_emm", "Z Leading Lepton; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "ZLep1Pt-max(ZLep2PtGen,ZLep1PtGen)", "EvtType == 2", *ZLep1Pt_emm);
  ZLep1Pt_emm->SetMarkerColor(kGreen); ZLep1Pt_emm->SetLineColor(kGreen); ZLep1Pt_emm->Draw("sames"); cLepRes->Update(); 
  stat2 = (TPaveStats*) (ZLep1Pt_emm->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kGreen);
  shiftStats(stat1, stat2);
  TH1F* ZLep1Pt_eem = new TH1F("ZLep1Pt_eem", "Z Leading Lepton; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "ZLep1Pt-max(ZLep2PtGen,ZLep1PtGen)", "EvtType == 1", *ZLep1Pt_eem);
  ZLep1Pt_eem->SetMarkerColor(kBlue); ZLep1Pt_eem->SetLineColor(kBlue); ZLep1Pt_eem->Draw("sames"); cLepRes->Update(); 
  stat3 = (TPaveStats*) (ZLep1Pt_eem->GetListOfFunctions()->FindObject("stats")); stat3->SetTextColor(kBlue);
  shiftStats(stat2, stat3);
  TH1F* ZLep1Pt_eee = new TH1F("ZLep1Pt_eee", "Z Leading Lepton; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "ZLep1Pt-max(ZLep2PtGen,ZLep1PtGen)", "EvtType == 0", *ZLep1Pt_eee);
  ZLep1Pt_eee->SetMarkerColor(kRed); ZLep1Pt_eee->SetLineColor(kRed); ZLep1Pt_eee->Draw("sames"); cLepRes->Update(); 
  stat4 = (TPaveStats*) (ZLep1Pt_eee->GetListOfFunctions()->FindObject("stats")); stat4->SetTextColor(kRed);
  shiftStats(stat3, stat4);
  ZLep1Pt_mmm->SetMaximum(findMax(ZLep1Pt_mmm, ZLep1Pt_emm, ZLep1Pt_eem, ZLep1Pt_eee));

  cLepRes->SaveAs((outName + "_ZLeadingLeptonPtResolution.png").c_str());
  ////////////////////////////////////

  cLepRes->cd(4);
  TH1F* ZLep2Pt_mmm = new TH1F("ZLep2Pt_mmm", "Z Trailing Lepton; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "ZLep2Pt-min(ZLep2PtGen,ZLep1PtGen)", "EvtType == 3", *ZLep2Pt_mmm);
  ZLep2Pt_mmm->SetMarkerColor(kYellow  ); ZLep2Pt_mmm->SetLineColor(kYellow  ); ZLep2Pt_mmm->Draw(); cLepRes->Update(); 
  stat1 = (TPaveStats*) (ZLep2Pt_mmm->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kYellow);
  TH1F* ZLep2Pt_emm = new TH1F("ZLep2Pt_emm", "Z Trailing Lepton; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "ZLep2Pt-min(ZLep2PtGen,ZLep1PtGen)", "EvtType == 2", *ZLep2Pt_emm);
  ZLep2Pt_emm->SetMarkerColor(kGreen); ZLep2Pt_emm->SetLineColor(kGreen); ZLep2Pt_emm->Draw("sames"); cLepRes->Update(); 
  stat2 = (TPaveStats*) (ZLep2Pt_emm->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kGreen);
  shiftStats(stat1, stat2);
  TH1F* ZLep2Pt_eem = new TH1F("ZLep2Pt_eem", "Z Trailing Lepton; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "ZLep2Pt-min(ZLep2PtGen,ZLep1PtGen)", "EvtType == 1", *ZLep2Pt_eem);
  ZLep2Pt_eem->SetMarkerColor(kBlue); ZLep2Pt_eem->SetLineColor(kBlue); ZLep2Pt_eem->Draw("sames"); cLepRes->Update(); 
  stat3 = (TPaveStats*) (ZLep2Pt_eem->GetListOfFunctions()->FindObject("stats")); stat3->SetTextColor(kBlue);
  shiftStats(stat2, stat3);
  TH1F* ZLep2Pt_eee = new TH1F("ZLep2Pt_eee", "Z Trailing Lepton; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "ZLep2Pt-min(ZLep2PtGen,ZLep1PtGen)", "EvtType == 0", *ZLep2Pt_eee);
  ZLep2Pt_eee->SetMarkerColor(kRed); ZLep2Pt_eee->SetLineColor(kRed); ZLep2Pt_eee->Draw("sames"); cLepRes->Update(); 
  stat4 = (TPaveStats*) (ZLep2Pt_eee->GetListOfFunctions()->FindObject("stats")); stat4->SetTextColor(kRed);
  shiftStats(stat3, stat4);
  ZLep2Pt_mmm->SetMaximum(findMax(ZLep2Pt_mmm, ZLep2Pt_emm, ZLep2Pt_eem, ZLep2Pt_eee));

  cLepRes->SaveAs((outName + "_ZTrailingLeptonPtResolution.png").c_str());
  ////////////////////////////////////

  TCanvas* cWZRes = new TCanvas("WZ Resolution");

  TH1F* WZRes_mmm = new TH1F("WZRes_mmm", "WZ Mass; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "WZMass-WprimeGenMass", "EvtType == 3", *WZRes_mmm);
  WZRes_mmm->SetMarkerColor(kYellow  ); WZRes_mmm->SetLineColor(kYellow  ); WZRes_mmm->Draw(); cWZRes->Update(); 
  stat1 = (TPaveStats*) (WZRes_mmm->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kYellow);
  TH1F* WZRes_emm = new TH1F("WZRes_emm", "WZ Mass; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "WZMass-WprimeGenMass", "EvtType == 2", *WZRes_emm);
  WZRes_emm->SetMarkerColor(kGreen); WZRes_emm->SetLineColor(kGreen); WZRes_emm->Draw("sames"); cWZRes->Update(); 
  stat2 = (TPaveStats*) (WZRes_emm->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kGreen);
  shiftStats(stat1, stat2);
  TH1F* WZRes_eem = new TH1F("WZRes_eem", "WZ Mass; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "WZMass-WprimeGenMass", "EvtType == 1", *WZRes_eem);
  WZRes_eem->SetMarkerColor(kBlue); WZRes_eem->SetLineColor(kBlue); WZRes_eem->Draw("sames"); cWZRes->Update(); 
  stat3 = (TPaveStats*) (WZRes_eem->GetListOfFunctions()->FindObject("stats")); stat3->SetTextColor(kBlue);
  shiftStats(stat2, stat3);
  TH1F* WZRes_eee = new TH1F("WZRes_eee", "WZ Mass; Resolution (Reco - Gen) (GeV); a.u.", 100, -500, 500); get_sum_of_hists(_file0, sample, "tEvts_MET", "WZMass-WprimeGenMass", "EvtType == 0", *WZRes_eee);
  WZRes_eee->SetMarkerColor(kRed); WZRes_eee->SetLineColor(kRed); WZRes_eee->Draw("sames"); cWZRes->Update(); 
  stat4 = (TPaveStats*) (WZRes_eee->GetListOfFunctions()->FindObject("stats")); stat4->SetTextColor(kRed);
  shiftStats(stat3, stat4);
  WZRes_mmm->SetMaximum(findMax(WZRes_mmm, WZRes_emm, WZRes_eem, WZRes_eee));

  cWZRes->SaveAs((outName + "_WZResolution.png").c_str());

  /////////////////////////////////////
  
  TCanvas* cResPt = new TCanvas("Resolution as Function of Pt");

  TH1F* WLepPt0_100_mmm = new TH1F("WLepPt0_100_mmm", "0 < W Lepton Pt < 100; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 0 && WLepPtGen < 100 && EvtType == 3", *WLepPt0_100_mmm);
  WLepPt0_100_mmm->SetMarkerColor(kYellow  ); WLepPt0_100_mmm->SetLineColor(kYellow  ); WLepPt0_100_mmm->Draw(); cResPt->Update(); 
  stat1 = (TPaveStats*) (WLepPt0_100_mmm->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kYellow);
  TH1F* WLepPt0_100_emm = new TH1F("WLepPt0_100_emm", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 0 && WLepPtGen < 100 && EvtType == 2", *WLepPt0_100_emm);
  WLepPt0_100_emm->SetMarkerColor(kGreen); WLepPt0_100_emm->SetLineColor(kGreen); WLepPt0_100_emm->Draw("sames"); cResPt->Update(); 
  stat2 = (TPaveStats*) (WLepPt0_100_emm->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kGreen);
  shiftStats(stat1, stat2); 
  TH1F* WLepPt0_100_eem = new TH1F("WLepPt0_100_eem", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 0 && WLepPtGen < 100 && EvtType == 1", *WLepPt0_100_eem);
  WLepPt0_100_eem->SetMarkerColor(kBlue); WLepPt0_100_eem->SetLineColor(kBlue); WLepPt0_100_eem->Draw("sames"); cResPt->Update(); 
  stat3 = (TPaveStats*) (WLepPt0_100_eem->GetListOfFunctions()->FindObject("stats")); stat3->SetTextColor(kBlue);
  shiftStats(stat2, stat3);
  TH1F* WLepPt0_100_eee = new TH1F("WLepPt0_100_eee", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 0 && WLepPtGen < 100 && EvtType == 0", *WLepPt0_100_eee);
  WLepPt0_100_eee->SetMarkerColor(kRed); WLepPt0_100_eee->SetLineColor(kRed); WLepPt0_100_eee->Draw("sames"); cResPt->Update(); 
  stat4 = (TPaveStats*) (WLepPt0_100_eee->GetListOfFunctions()->FindObject("stats")); stat4->SetTextColor(kRed);
  shiftStats(stat3, stat4);
  WLepPt0_100_mmm->SetMaximum(findMax(WLepPt0_100_mmm, WLepPt0_100_emm, WLepPt0_100_eem, WLepPt0_100_eee));

  cResPt->SaveAs((outName + "_WLeptonPtResolution_0-100.png").c_str());
  /////////////////////////////////////

  TH1F* WLepPt100_200_mmm = new TH1F("WLepPt100_200_mmm", "100 < W Lepton Pt < 200; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 100 && WLepPtGen < 200 && EvtType == 3", *WLepPt100_200_mmm);
  WLepPt100_200_mmm->SetMarkerColor(kYellow  ); WLepPt100_200_mmm->SetLineColor(kYellow  ); WLepPt100_200_mmm->Draw(); cResPt->Update(); 
  stat1 = (TPaveStats*) (WLepPt100_200_mmm->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kYellow);
  TH1F* WLepPt100_200_emm = new TH1F("WLepPt100_200_emm", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 100 && WLepPtGen < 200 && EvtType == 2", *WLepPt100_200_emm);
  WLepPt100_200_emm->SetMarkerColor(kGreen); WLepPt100_200_emm->SetLineColor(kGreen); WLepPt100_200_emm->Draw("sames"); cResPt->Update(); 
  stat2 = (TPaveStats*) (WLepPt100_200_emm->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kGreen);
  shiftStats(stat1, stat2); 
  TH1F* WLepPt100_200_eem = new TH1F("WLepPt100_200_eem", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 100 && WLepPtGen < 200 && EvtType == 1", *WLepPt100_200_eem);
  WLepPt100_200_eem->SetMarkerColor(kBlue); WLepPt100_200_eem->SetLineColor(kBlue); WLepPt100_200_eem->Draw("sames"); cResPt->Update(); 
  stat3 = (TPaveStats*) (WLepPt100_200_eem->GetListOfFunctions()->FindObject("stats")); stat3->SetTextColor(kBlue);
  shiftStats(stat2, stat3);
  TH1F* WLepPt100_200_eee = new TH1F("WLepPt100_200_eee", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 100 && WLepPtGen < 200 && EvtType == 0", *WLepPt100_200_eee);
  WLepPt100_200_eee->SetMarkerColor(kRed); WLepPt100_200_eee->SetLineColor(kRed); WLepPt100_200_eee->Draw("sames"); cResPt->Update(); 
  stat4 = (TPaveStats*) (WLepPt100_200_eee->GetListOfFunctions()->FindObject("stats")); stat4->SetTextColor(kRed);
  shiftStats(stat3, stat4);
  WLepPt100_200_mmm->SetMaximum(findMax(WLepPt100_200_mmm, WLepPt100_200_emm, WLepPt100_200_eem, WLepPt100_200_eee));

  cResPt->SaveAs((outName + "_WLeptonPtResolution_100-200.png").c_str());
  /////////////////////////////////////

  TH1F* WLepPt200_400_mmm = new TH1F("WLepPt200_400_mmm", "200 < W Lepton Pt < 400; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 200 && WLepPtGen < 400 && EvtType == 3", *WLepPt200_400_mmm);
  WLepPt200_400_mmm->SetMarkerColor(kYellow  ); WLepPt200_400_mmm->SetLineColor(kYellow  ); WLepPt200_400_mmm->Draw(); cResPt->Update(); 
  stat1 = (TPaveStats*) (WLepPt200_400_mmm->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kYellow);
  TH1F* WLepPt200_400_emm = new TH1F("WLepPt200_400_emm", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 200 && WLepPtGen < 400 && EvtType == 2", *WLepPt200_400_emm);
  WLepPt200_400_emm->SetMarkerColor(kGreen); WLepPt200_400_emm->SetLineColor(kGreen); WLepPt200_400_emm->Draw("sames"); cResPt->Update(); 
  stat2 = (TPaveStats*) (WLepPt200_400_emm->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kGreen);
  shiftStats(stat1, stat2); 
  TH1F* WLepPt200_400_eem = new TH1F("WLepPt200_400_eem", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 200 && WLepPtGen < 400 && EvtType == 1", *WLepPt200_400_eem);
  WLepPt200_400_eem->SetMarkerColor(kBlue); WLepPt200_400_eem->SetLineColor(kBlue); WLepPt200_400_eem->Draw("sames"); cResPt->Update(); 
  stat3 = (TPaveStats*) (WLepPt200_400_eem->GetListOfFunctions()->FindObject("stats")); stat3->SetTextColor(kBlue);
  shiftStats(stat2, stat3);
  TH1F* WLepPt200_400_eee = new TH1F("WLepPt200_400_eee", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 200 && WLepPtGen < 400 && EvtType == 0", *WLepPt200_400_eee);
  WLepPt200_400_eee->SetMarkerColor(kRed); WLepPt200_400_eee->SetLineColor(kRed); WLepPt200_400_eee->Draw("sames"); cResPt->Update(); 
  stat4 = (TPaveStats*) (WLepPt200_400_eee->GetListOfFunctions()->FindObject("stats")); stat4->SetTextColor(kRed);
  shiftStats(stat3, stat4);
  WLepPt200_400_mmm->SetMaximum(findMax(WLepPt200_400_mmm, WLepPt200_400_emm, WLepPt200_400_eem, WLepPt200_400_eee));

  cResPt->SaveAs((outName + "_WLeptonPtResolution_200-400.png").c_str());
  /////////////////////////////////////

  TH1F* WLepPt400_5000_mmm = new TH1F("WLepPt400_5000_mmm", "400 < W Lepton < 5000; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 400 && WLepPtGen < 5000 && EvtType == 3", *WLepPt400_5000_mmm);
  WLepPt400_5000_mmm->SetMarkerColor(kYellow  ); WLepPt400_5000_mmm->SetLineColor(kYellow  ); WLepPt400_5000_mmm->Draw(); cResPt->Update(); 
  stat1 = (TPaveStats*) (WLepPt400_5000_mmm->GetListOfFunctions()->FindObject("stats")); stat1->SetTextColor(kYellow);
  TH1F* WLepPt400_5000_emm = new TH1F("WLepPt400_5000_emm", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 400 && WLepPtGen < 5000 && EvtType == 2", *WLepPt400_5000_emm);
  WLepPt400_5000_emm->SetMarkerColor(kGreen); WLepPt400_5000_emm->SetLineColor(kGreen); WLepPt400_5000_emm->Draw("sames"); cResPt->Update(); 
  stat2 = (TPaveStats*) (WLepPt400_5000_emm->GetListOfFunctions()->FindObject("stats")); stat2->SetTextColor(kGreen);
  shiftStats(stat1, stat2); 
  TH1F* WLepPt400_5000_eem = new TH1F("WLepPt400_5000_eem", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 400 && WLepPtGen < 5000 && EvtType == 1", *WLepPt400_5000_eem);
  WLepPt400_5000_eem->SetMarkerColor(kBlue); WLepPt400_5000_eem->SetLineColor(kBlue); WLepPt400_5000_eem->Draw("sames"); cResPt->Update(); 
  stat3 = (TPaveStats*) (WLepPt400_5000_eem->GetListOfFunctions()->FindObject("stats")); stat3->SetTextColor(kBlue);
  shiftStats(stat2, stat3);
  TH1F* WLepPt400_5000_eee = new TH1F("WLepPt400_5000_eee", "W Lepton; Resolution (Reco - Gen) (GeV); a.u.", 200, -200, 200); get_sum_of_hists(_file0, sample, "tEvts_MET", "WLepPt-WLepPtGen", "WLepPtGen > 400 && WLepPtGen < 5000 && EvtType == 0", *WLepPt400_5000_eee);
  WLepPt400_5000_eee->SetMarkerColor(kRed); WLepPt400_5000_eee->SetLineColor(kRed); WLepPt400_5000_eee->Draw("sames"); cResPt->Update(); 
  stat4 = (TPaveStats*) (WLepPt400_5000_eee->GetListOfFunctions()->FindObject("stats")); stat4->SetTextColor(kRed);
  shiftStats(stat3, stat4);
  WLepPt400_5000_mmm->SetMaximum(findMax(WLepPt400_5000_mmm, WLepPt400_5000_emm, WLepPt400_5000_eem, WLepPt400_5000_eee));

  cResPt->SaveAs((outName + "_WLeptonPtResolution_400-5000.png").c_str());
  /////////////////////////////////////

}


void 
shiftStats(TPaveStats* stat1, TPaveStats* stat2){
  float height = stat1->GetY2NDC() - stat1->GetY1NDC();
  stat2->SetY2NDC(stat1->GetY1NDC() );
  stat2->SetY1NDC(stat2->GetY2NDC() - height);
  stat2->Draw();
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
