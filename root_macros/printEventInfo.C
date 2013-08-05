//Author: Cory Fantasia 2013
//Purpose: Print important event info for select events
//Usage: root -b -l -q printEventInfo.C

#include "TFile.h"
#include "TTree.h"

void
printEventInfo(){
  TFile *_file0 = TFile::Open("../../../WprimeWZ.root");
  TTree* tEvts = (TTree*) _file0->Get("data/tEvts_MET");
  int n = tEvts->Draw("Run:Lumi:Event:WZMass:ZMass:WTransMass:Zpt:Wpt:Lt:ZLep1Pt:ZLep1Eta:ZLep1Phi:ZLep2Pt:ZLep2Eta:ZLep2Phi:WLepPt:WLepEta:WLepPhi:MET:METPhi:EvtType", 
                      "(Lt>500 || MET>300 || WZMass>1000)", "para goff");

  for(int ievt=0; ievt<n; ++ievt){
    int treeIdx = 0;
    const double Run   = tEvts->GetVal(treeIdx++)[ievt];
    const double Lumi  = tEvts->GetVal(treeIdx++)[ievt];
    const double Event = tEvts->GetVal(treeIdx++)[ievt];
    const double WZMass   = tEvts->GetVal(treeIdx++)[ievt];
    const double ZMass   = tEvts->GetVal(treeIdx++)[ievt];
    const double WTransMass   = tEvts->GetVal(treeIdx++)[ievt];
    const double Zpt   = tEvts->GetVal(treeIdx++)[ievt];
    const double Wpt   = tEvts->GetVal(treeIdx++)[ievt];
    const double Lt   = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep1Pt   = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep1Eta   = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep1Phi   = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2Pt   = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2Eta   = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2Phi   = tEvts->GetVal(treeIdx++)[ievt];
    const double WLepPt   = tEvts->GetVal(treeIdx++)[ievt];
    const double WLepEta   = tEvts->GetVal(treeIdx++)[ievt];
    const double WLepPhi   = tEvts->GetVal(treeIdx++)[ievt];
    const double MET   = tEvts->GetVal(treeIdx++)[ievt];
    const double METPhi   = tEvts->GetVal(treeIdx++)[ievt];
    const int EvtType   = round(tEvts->GetVal(treeIdx++)[ievt]);
    
    printf(" Run & Lumi & Event \\\\ \n");
    printf("%.0f & %.0f & %.0f \\\\ \n", Run, Lumi, Event);

    printf("$M_{WZ}$ (\\GeV)& $M_Z$ (\\GeV) & $M_W^T$ (\\GeV) \\\\ \n");
    printf("%.0f & %.0f & %.0f \\\\ \n", WZMass, ZMass, WTransMass);

    printf("$p_{T}^{Z}$ (\\GeV) & $p_{T}^{W}$ (\\GeV) & $L_{T}$ (\\GeV) \\\\ \n");
    printf("%.0f & %.0f & %.0f \\\\ \n", Zpt, Wpt, Lt);

    printf("$p_{T}^{%s}$ & $\\eta$ & $\\phi$ \\\\ \n", EvtType>=1.5 ? "\\mu" : "e");
    printf("%.0f & %.1f & %.1f \\\\ \n", ZLep1Pt, ZLep1Eta, ZLep1Phi);
    
    printf("$p_{T}^{%s}$ & $\\eta$ & $\\phi$ \\\\ \n", EvtType>=1.5 ? "\\mu" : "e");
    printf("%.0f & %.1f & %.1f \\\\ \n", ZLep2Pt, ZLep2Eta, ZLep2Phi);
    
    printf("$p_{T}^{%s}$ & $\\eta$ & $\\phi$ \\\\ \n", EvtType%2 ? "\\mu" : "e");
    printf("%.0f & %.1f & %.1f \\\\ \n", WLepPt, WLepEta, WLepPhi);
    
    printf("$\\MET$ & $\\eta$ & $\\phi$ \\\\ \n");
    printf("%.0f & - & %.1f \\\\ \n", MET, METPhi);

    printf("\\hline");
  }

}
