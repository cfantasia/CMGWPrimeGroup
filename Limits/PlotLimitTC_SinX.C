#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "../root_macros/CMSStyle.C"
#include "../root_macros/setTDR_modified.C"
#include "consts.h"
#include "TGraph2D.h"

const float minSINX = 0.20;
const float maxSINX = 0.50;
const float maxRHO = 300;
const float minRHO = 280;
const float RHO_INC = 1;
const float SINX_INC = 0.005;

void PlotLimitTC_SinX(float MPi=160) {
  setTDRStyle();
  gROOT->ForceStyle();

  cout<<"Assuming Pi Mass = "<<MPi<<endl;

  TCanvas* c1 = new TCanvas("TC-RhoVsSinX","");

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
  tXsec->ReadFile("list.dat");
  //tXsec->ReadFile("xSec_TCWZ-RhoSinX.dat");
  tXsec->Draw("Rho:SinX:Xsec", Form("Pi==%.0f",MPi), "para goff");
  int nXsec = tXsec->GetSelectedRows(); assert(nXsec); 
  Double_t* gRho = tXsec->GetVal(0);
  Double_t* gSinX = tXsec->GetVal(1);
  Double_t* gXsec = tXsec->GetVal(2);
  TGraph2D *gXsec2D = new TGraph2D(nXsec, gRho, gSinX, gXsec);
  cout<<"Done with 2D xsec graph. There were "<<nXsec<<endl;

  vector<float> masses, expSinXLims, obsSinXLims;
  for(float rho=minRHO; rho<=maxRHO; rho+=RHO_INC){
    float limExp = maxSINX;
    float limObs = maxSINX;
    masses.push_back(rho);
    for(float sinx=minSINX; sinx<=maxSINX; sinx+=SINX_INC){
      if(gXsec2D->Interpolate(rho,sinx) > gExpLim->Eval(rho)){
        limExp = sinx;
        break;
      }
    }
    expSinXLims.push_back(limExp);

    for(float sinx=minSINX; sinx<=maxSINX; sinx+=SINX_INC){
      //cout<<"rho:sinX:Lim:xsec = "<<rho<<":"<<sinx<<":"<<gObsLim->Eval(rho)<<":"<<gXsec2D->Interpolate(rho,sinx)<<endl;
      if(gXsec2D->Interpolate(rho,sinx) > gObsLim->Eval(rho)){
        limObs = sinx;
        break;
      }
    }
    obsSinXLims.push_back(limObs);
  }
  cout<<"Done with looping "<<endl;
  TGraph* exp = new TGraph(expSinXLims.size()+2);//n+2 to close the curve
  TGraph* obs = new TGraph(expSinXLims.size()+2);


  for(unsigned i=0; i<expSinXLims.size(); ++i){///////
    exp->SetPoint(i, masses[i], expSinXLims[i]);
    obs->SetPoint(i, masses[i], obsSinXLims[i]);
    //cout<<"mass:expLim:obsLim:ExpSinX:ObsSinX = "<<masses[i]<<":"<<gExpLim->Eval(masses[i])<<":"<<gObsLim->Eval(masses[i])<<"\t"<<expSinXLims[i]<<"\t\t"<<obsSinXLims[i]<<"\t\t "<<gXsec2D->Interpolate(masses[i],obsSinXLims[i])<<"\t"<<gXsec2D->Interpolate(masses[i],obsSinXLims[i]-SINX_INC)<<endl;  
  }

  /*
  ///////Find Optimistic Limits///////////
  //Expected
  for(int i=expSinXLims.size()-1; i>=0; --i){///////
    if(masses[i] > 350) continue;
    if(expSinXLims[i] > masses[i] - 80.4){
      cout<<"Lower Optimistic Exp Limit is "<<masses[i]<<endl;
      break;
    }
  }
  for(unsigned i=0; i<expSinXLims.size(); ++i){///////
    if(masses[i] < 350) continue;
    if(expSinXLims[i] > masses[i] - 80.4){
      cout<<"Upper Optimistic Exp Limit is "<<masses[i]<<endl;
      break;
    }
  }
  //Observed
  for(int i=obsSinXLims.size()-1; i>=0; --i){///////
    if(masses[i] > 350) continue;
    if(obsSinXLims[i] > masses[i] - 80.4){
      cout<<"Lower Optimistic Obs Limit is "<<masses[i]<<endl;
      break;
    }
  }
  for(unsigned i=0; i<obsSinXLims.size(); ++i){///////
    if(masses[i] < 350) continue;
    if(obsSinXLims[i] > masses[i] - 80.4){
      cout<<"Upper Optimistic Obs Limit is "<<masses[i]<<endl;
      break;
    }
  }

  /////Les Houches Param/////////////
  //Expected
  for(int i=expSinXLims.size()-1; i>=0; --i){///////
    if(masses[i] > 350) continue;
    if(expSinXLims[i] > 0.75 * masses[i] - 25){
      cout<<"Lower Les Houches Exp Limit is "<<masses[i]<<endl;
      break;
    }
  }
  for(unsigned i=0; i<expSinXLims.size(); ++i){///////
    if(masses[i] < 350) continue;
    if(expSinXLims[i] > 0.75 * masses[i] - 25){
      cout<<"Upper Les Houches Exp Limit is "<<masses[i]<<endl;
      break;
    }
  }
  //Observed
  for(int i=obsSinXLims.size()-1; i>=0; --i){///////
    if(masses[i] > 350) continue;
    if(obsSinXLims[i] > 0.75 * masses[i] - 25){
      cout<<"Lower Les Houches Obs Limit is "<<masses[i]<<endl;
      break;
    }
  }
  for(unsigned i=0; i<obsSinXLims.size(); ++i){///////
    if(masses[i] < 350) continue;
    if(obsSinXLims[i] > 0.75 * masses[i] - 25){
      cout<<"Upper Les Houches Obs Limit is "<<masses[i]<<endl;
      break;
    }
  }

  */

  /////CDF Bump Param/////////

  for(unsigned i=0; i<obsSinXLims.size(); ++i){///////
    if(masses[i] == 280) cout<<"For rho="<<masses[i]<<", sinx limit is "<<obsSinXLims[i]
                             <<" (obs) and "<<expSinXLims[i]<<" (exp) "
                             <<gObsLim->Eval(masses[i])<<" pb (obs) and "<<gExpLim->Eval(masses[i])<<" pb (exp)"<<endl;
    if(masses[i] == 290) cout<<"For rho="<<masses[i]<<", sinx limit is "<<obsSinXLims[i]
                             <<" (obs) and "<<expSinXLims[i]<<" (exp) "
                             <<gObsLim->Eval(masses[i])<<" pb (obs) and "<<gExpLim->Eval(masses[i])<<" pb (exp)"<<endl;
  }
  

  exp->SetPoint(expSinXLims.size(), maxRHO, maxSINX);
  obs->SetPoint(expSinXLims.size(), maxRHO, maxSINX);
  exp->SetPoint(expSinXLims.size()+1, minRHO, maxSINX);
  obs->SetPoint(expSinXLims.size()+1, minRHO, maxSINX);

  obs->SetFillStyle(3003);
  exp->SetFillStyle(3004);
  obs->SetFillColor(kRed);
  exp->SetFillColor(4);
  obs->SetLineColor(kRed);
  exp->SetLineColor(4);

  TMultiGraph *mg = new TMultiGraph("mg", ";M(#rho_{TC}) (GeV);sin #chi");
  mg->Add(obs,"F");
  mg->Add(exp,"F");
  mg->Add(obs,"C");
  mg->Add(exp,"C");
  mg->SetMinimum(minSINX);
  mg->SetMaximum(maxSINX);
  mg->Draw("a");
  mg->GetXaxis()->SetNdivisions(505);
  mg->GetYaxis()->SetNdivisions(505);

  TGraph* cdfPoint = new TGraph(1);
  cdfPoint->SetMarkerColor(kGreen);
  cdfPoint->SetMarkerStyle(20);
  cdfPoint->SetPoint(0,290, 160);
  //cdfPoint->Draw("p");

  TLegend* leg = new TLegend(0.30, 0.75, 0.60, 0.90);
  leg->SetTextSize(0.05);
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(exp, "Exp. Limit", "f");
  leg->AddEntry(obs, "Obs. Limit", "f");
  //leg->AddEntry(cdfPoint, "CDF Bump", "p");
  leg->Draw();

  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextSize(0.05);
  latexLabel.SetTextFont(42);

  latexLabel.DrawLatex(0.33, 0.96, "CMS Preliminary 2011");
  latexLabel.DrawLatex(0.67, 0.8, "#sqrt{s} = 7 TeV");
  latexLabel.DrawLatex(0.62, 0.65, Form("#intL dt = %.1f fb^{-1}",lumi[0]/1000.));

  TLine* line1 = new TLine(minRHO, 0.3333, maxRHO, 0.3333);
  line1->SetLineStyle(kDashed);
  line1->SetLineWidth(2);
  line1->Draw();
  latexLabel.DrawLatex(0.25, 0.60, "Default: sin #chi = 1/3");


  latexLabel.DrawLatex(0.25, 0.45, Form("M(#pi_{TC}) = %.0f GeV", MPi));

  c1->RedrawAxis();
  c1->Print(Form("tcLimit-SinX-MPi%.0f.pdf",MPi));
  //c1->Print("tcLimit-SinX.png");

  TFile *fin = TFile::Open("WprimeWZ-Limits.root", "update"); assert(fin);
  c1->Write();
  fin->Close();
}
