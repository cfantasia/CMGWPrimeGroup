#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "../root_macros/CMSStyle.C"
#include "consts.h"
#include "TGraph2D.h"

const float minPI = 60;
const float maxPI = 1000;
const float maxRHO = 1000;
const float minRHO = 150;
const int RHO_INC = 1;
const int PI_INC = 1;

void PlotLimit2D() {
  gErrorIgnoreLevel = kWarning;
  CMSstyle();

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
  cout<<"Done with importing limits.  There were "<<n<<endl;

  //Fill Xsec 2D graph
  TTree* tXsec = new TTree("tXsec", "X Sec");
  tXsec->ReadFile("xSec_TCWZ.dat");
  //tXsec->ReadFile("xSec_TCWZ-sinchi1d2.dat");
  tXsec->Draw("Rho:Pi:Xsec", "", "para goff");
  int nXsec = tXsec->GetSelectedRows(); assert(nXsec); 
  Double_t* gRho = tXsec->GetVal(0);
  Double_t* gPi =  tXsec->GetVal(1);
  Double_t* gXsec = tXsec->GetVal(2);
  TGraph2D *gXsec2D = new TGraph2D(nXsec, gRho, gPi, gXsec);
  cout<<"Done with 2D xsec graph. There were "<<nXsec<<endl;

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


  for(unsigned i=0; i<expPiLims.size(); ++i){///////
    exp->SetPoint(i, masses[i], expPiLims[i]);
    obs->SetPoint(i, masses[i], obsPiLims[i]);
    //cout<<"mass:expLim:obsLim:ExpPi:ObsPi = "<<masses[i]<<":"<<gExpLim->Eval(masses[i])<<":"<<gObsLim->Eval(masses[i])<<"\t"<<expPiLims[i]<<"\t\t"<<obsPiLims[i]<<endl;  
  }

  //Find Optimistic Limits
  for(unsigned i=0; i<expPiLims.size(); ++i){///////
    if(masses[i] < 350) continue;
    if(expPiLims[i] > masses[i] - 80.4){
      cout<<"Upper Optimistic Exp Limit is "<<masses[i]<<endl;
      break;
    }
  }
  for(unsigned i=expPiLims.size()-1; i>=0; --i){///////
    if(masses[i] > 350) continue;
    if(expPiLims[i] > masses[i] - 80.4){
      cout<<"Lower Optmistic Exp Limit is "<<masses[i]<<endl;
      break;
    }
  }
  
  for(unsigned i=0; i<obsPiLims.size(); ++i){///////
    if(masses[i] < 350) continue;
    if(obsPiLims[i] > masses[i] - 80.4){
      cout<<"Upper Optimistic Obs Limit is "<<masses[i]<<endl;
      break;
    }
  }
  for(unsigned i=obsPiLims.size()-1; i>=0; --i){///////
    if(masses[i] > 350) continue;
    if(obsPiLims[i] > masses[i] - 80.4){
      cout<<"Lower Optimistic Obs Limit is "<<masses[i]<<endl;
      break;
    }
  }

  /////Les Houches Param/////////////

  for(unsigned i=0; i<expPiLims.size(); ++i){///////
    if(masses[i] < 350) continue;
    if(expPiLims[i] > 0.75 * masses[i] - 25){
      cout<<"Upper Les Houches Exp Limit is "<<masses[i]<<endl;
      break;
    }
  }
  for(unsigned i=expPiLims.size()-1; i>=0; --i){///////
    if(masses[i] > 350) continue;
    if(expPiLims[i] > 0.75 * masses[i] - 25){
      cout<<"Lower Les Houches Exp Limit is "<<masses[i]<<endl;
      break;
    }
  }

  for(unsigned i=0; i<obsPiLims.size(); ++i){///////
    if(masses[i] < 350) continue;
    if(obsPiLims[i] > 0.75 * masses[i] - 25){
      cout<<"Upper Les Houches Obs Limit is "<<masses[i]<<endl;
      break;
    }
  }
  for(unsigned i=obsPiLims.size()-1; i>=0; --i){///////
    if(masses[i] > 350) continue;
    if(obsPiLims[i] > 0.75 * masses[i] - 25){
      cout<<"Lower Les Houches Obs Limit is "<<masses[i]<<endl;
      break;
    }
  }



  /////CDF Bump Param/////////

  for(unsigned i=0; i<obsPiLims.size(); ++i){///////
    if(masses[i] == 290) cout<<"For rho="<<masses[i]<<", pi limit is "<<obsPiLims[i]
                             <<" (obs) and "<<expPiLims[i]<<" (exp)"<<endl;
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
  
  TH2F* frame = new TH2F("frame", "", 100, minRHO, maxRHO, 100, minPI, maxPI);
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
  ax->SetTitle("M(#rho_{TC}) (GeV)");
  ay->SetTitle("M(#pi_{TC}) (GeV)");
  frame->Draw();
  obs->Draw("F");
  exp->Draw("F");
  exp->Draw("L");
  obs->Draw("L");

  TGraph* cdfPoint = new TGraph(1);
  cdfPoint->SetMarkerColor(kGreen);
  cdfPoint->SetMarkerStyle(20);
  cdfPoint->SetPoint(0,290, 160);
  //cdfPoint->Draw("p");

  TLegend* leg = new TLegend(0.30, 0.75, 0.60, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(exp, "Exp. Limit", "fl");
  leg->AddEntry(obs, "Obs. Limit", "fl");
  //leg->AddEntry(cdfPoint, "CDF Bump", "p");
  leg->Draw();

  TLatex* text0 = new TLatex(550, 200, "CMS Preliminary 2011");
  text0->SetTextSize(0.05);
  text0->Draw();

  TLatex* text1 = new TLatex(660, 130, "#sqrt{s} = 7 TeV");
  text1->SetTextSize(0.05);
  text1->Draw();

  TLatex* text2 = new TLatex(610, 285, Form("#int L dt = %.1f fb^{-1}",lumi[0]/1000));
  text2->SetTextSize(0.05);
  text2->Draw();

  TLine* line1 = new TLine(minRHO, minRHO-80.4, maxRHO, maxRHO-80.4);
  line1->SetLineStyle(2);
  line1->SetLineWidth(2);
  line1->Draw();

  TLine* line2 = new TLine(100, 50, 900, 650);
  line2->SetLineStyle(1);
  line2->SetLineWidth(2);
  //line2->Draw();

  TLatex* text3 = new TLatex(317, 265, "M(#pi_{TC}) = M(#rho_{TC}) - M(W)");
  text3->SetTextSize(0.06);
  text3->SetTextAngle(38);
  text3->Draw(); 

  c1->RedrawAxis();
  c1->Print("tcLimit-2D.pdf");
  c1->Print("tcLimit-2D.C");

  cout<<"Extra points:\n";
  for(int i=100; i<=900; i+=25){
    cout<<i<<" "<<i*3./4. - 25<<" "<<gXsec2D->Interpolate(i, i*3./4. - 25)<<endl;
  }

}
