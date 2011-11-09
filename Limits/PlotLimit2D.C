#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "../root_macros/CMSStyle.C"
#include "consts.h"
#include "TGraph2D.h"

const float minPI = 150;
const float maxPI = 500;
const float maxRHO = 700;
const float minRHO = 300;
const int RHO_INC = 1;
const int PI_INC = 1;

void PlotLimit2D() {
  TCanvas* c1 = new TCanvas("c1","",600,500);
  gStyle->SetOptStat(0);

  TTree* tLimit = new TTree("tLimit", "Limits");
  tLimit->ReadFile("nLimit.txt");
  tLimit->Draw("SignalCode:Mass:Lumi:ObsLimit:ExpLimit", "", "para goff");
  int n = tLimit->GetSelectedRows(); assert(n); 
  const Double_t* mass = tLimit->GetVal(1);
  const Double_t* lumi = tLimit->GetVal(2);
  const Double_t* ObsLimit = tLimit->GetVal(3);
  const Double_t* ExpLimit = tLimit->GetVal(4);
  TGraph* gExpLim = new TGraph(n, mass, ExpLimit);
  TGraph* gObsLim = new TGraph(n, mass, ObsLimit);
  cout<<"Done with importing limits "<<endl;

  //Fill Xsec 2D graph
  TTree* tXsec = new TTree("tXsec", "X Sec");
  tXsec->ReadFile("xSec_TCWZ.dat");
  tXsec->Draw("Rho:Pi:Xsec", "", "para goff");
  int nXsec = tXsec->GetSelectedRows(); assert(nXsec); 
  Double_t* gRho = tXsec->GetVal(0);
  Double_t* gPi =  tXsec->GetVal(1);
  Double_t* gXsec = tXsec->GetVal(2);
  TGraph2D *gXsec2D = new TGraph2D(nXsec, gRho, gPi, gXsec);
  cout<<"Done with 2D xsec graph "<<endl;

  vector<float> masses, expPiLims, obsPiLims;
  for(int rho=minRHO; rho<=maxRHO; rho+=RHO_INC){
    float limExp = maxPI;
    float limObs = maxPI;
    masses.push_back(rho);
    for(int pi=minPI; pi<=maxPI; pi+=PI_INC){
      if(gXsec2D->Interpolate(rho,pi) > gExpLim->Eval(rho)){
        limExp = pi;
        break;
      }
    }
    expPiLims.push_back(limExp);

    for(int pi=minPI; pi<=maxPI; pi+=PI_INC){
      if(gXsec2D->Interpolate(rho,pi) > gObsLim->Eval(rho)){
        limObs = pi;
        break;
      }
    }
    obsPiLims.push_back(limObs);
  }
  cout<<"Done with looping "<<endl;
  TGraph* exp = new TGraph(expPiLims.size()+2);//n+2 to close the curve
  TGraph* obs = new TGraph(expPiLims.size()+2);


  for(int i=0; i<expPiLims.size(); ++i){///////
    exp->SetPoint(i, masses[i], expPiLims[i]);
    obs->SetPoint(i, masses[i], obsPiLims[i]);

    cout<<"mass:ExpPi:ObsPi = "<<masses[i]<<"\t"<<expPiLims[i]<<"\t\t"<<obsPiLims[i]<<endl;
  }
  exp->SetPoint(expPiLims.size(), maxRHO, maxPI);
  obs->SetPoint(expPiLims.size(), maxRHO, maxPI);
  exp->SetPoint(expPiLims.size()+1, minRHO, maxPI);
  obs->SetPoint(expPiLims.size()+1, minRHO, maxPI);

  obs->SetFillStyle(3003);
  exp->SetFillStyle(3004);
  obs->SetFillColor(kRed);
  exp->SetFillColor(4);
  obs->SetLineColor(kRed);
  exp->SetLineColor(4);
  
  TH2F* frame = new TH2F("frame", "", 100, minRHO, maxRHO, 100, 150, maxPI);
  TAxis* ax = frame->GetXaxis();
  TAxis* ay = frame->GetYaxis();
  ax->SetLabelFont(132);
  ax->SetTitleFont(132);
  ay->SetLabelFont(132);
  ax->SetTitleFont(132);
  ax->SetLabelSize(0.05);
  ay->SetLabelSize(0.05);
  ax->SetTitleSize(0.06);
  ay->SetTitleSize(0.06);
  ax->SetTitle("m(#rho_{T}) (GeV)");
  ay->SetTitle("m(#pi_{T}) (GeV)");
  frame->Draw();
  obs->Draw("F");
  exp->Draw("F");
  exp->Draw("L");
  obs->Draw("L");

  TLegend* leg = new TLegend(0.50, 0.20, 0.90, 0.36);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(exp, "95% C.L. limit (exp)", "fl");
  leg->AddEntry(obs, "95% C.L. limit (obs)", "fl");
  leg->Draw();


  TLatex* text = new TLatex(400, 180, "CMS Preliminary 2011 #sqrt{s} = 7 TeV");
  text->SetTextSize(0.05);
  text->Draw();

  TLatex* text2 = new TLatex(520, 280, Form("#int L dt = %.2f fb^{-1}",lumi[0]/1000));
  text2->SetTextSize(0.05);
  text2->Draw();

  TLine* line1 = new TLine(minRHO, minRHO - 80.4, maxPI+80.4, maxPI);
  line1->SetLineStyle(2);
  line1->SetLineWidth(2);
  line1->Draw();
/*
  TLine* line2 = new TLine(300, 200, 700, 500);
  line2->SetLineStyle(1);
  line2->SetLineWidth(2);
  line2->Draw();
*/
  TLatex* text3 = new TLatex(317, 249, "m(#pi_{T}) = m(#rho_{T}) - m(W)");
  text3->SetTextSize(0.06);
  text3->SetTextAngle(43);
  text3->Draw(); 

  c1->RedrawAxis();
  c1->Print("tcLimit.eps");
  c1->Print("tcLimit.gif");

}
/*
float 
  FindLimit(const TGraph2D* xsec, const TGraph* limit, const float def=-1){
  float lower = min(xsec->GetX()[0], limit->GetX()[0]);
  float upper = max(xsec->GetX()[xsec->GetN()-1], limit->GetX()[limit->GetN()-1]);
  const float INC = 1;
  for(float mass=lower; mass<upper; mass+=INC){
    if(xsec->Interpolate(rho,pi) > limit->Eval(rho)) return pi;
}
*/
