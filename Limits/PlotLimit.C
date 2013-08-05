//Usage: root -l PlotLimit.C++

#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include "../root_macros/CMSStyle.C"
#include "../root_macros/setTDR_modified.C"
#include "consts.h"

void RemoveZeros(TGraph* g);
int GetEndPoints(int & start, int& end, Double_t* y, float line);
float FindLimit(const TGraph* xsec, const TGraph* limit, const bool findUpper);

struct SignalSample{
  string name;
  string sigXsec;
  string massString;
  string cutString;
  string legendString;
  TTree* tXsec;
  int lineColor;
  int lineStyle;
  int bandColor;
  TGraph* gXsec;
  TGraph* gXsecBand;

  SignalSample(){}
  SignalSample(string n, string xsec, string mass, string cut, string leg, int lc=kBlack, int ls=1, int bc=-1){
    name = n;
    sigXsec = xsec;
    massString = mass;
    cutString = cut;
    legendString = leg;
    tXsec = NULL;
    lineColor = lc;
    lineStyle = ls;
    bandColor = bc==-1 ? lc - 7 : bc;
    gXsec = NULL;
    gXsecBand = NULL;
  }
};

void
PlotLimit(string inName, string inFile="nLimit.txt", string outFile="", float xSecScale=1.){
  //gErrorIgnoreLevel = kWarning;
  //CMSstyle();
  setTDRStyle();
  gROOT->ForceStyle();

  vector<SignalSample> sigs;
  if(inName.find("WprimeWZ") != string::npos){
    sigs.push_back(SignalSample("W'","xSec_WZ.dat",  "Mass:Xsec:percentError",  "Mass>=172", "\\sigma_{W'}", kBlue));
    sigs.push_back(SignalSample("W'LO","test.dat",  "Mass:Xsec:percentError",  "Mass>=172", "\\sigma_{W'}^{LO}", kGreen));
    //sigs.push_back(SignalSample("W'","xSec_WZ.dat",  "Mass:Xsec:percentError",  "Mass>=172", "\\sigma_{W'}", kBlack, 1, kGray));
    //sigs.push_back(SignalSample("TC,sin(#chi)=#frac{1}{2}","xSec_TCWZ-sinchi1d2.dat",  "Rho:Xsec",  "Rho>=200",  "\\sigma_{TC sin(#chi)=#frac{1}{2}}", kRed, kDashed));
    sigs.push_back(SignalSample("TC",                      "xSec_TCWZ-sinchi1d3.dat",  "Rho:Xsec:percentError",  "Rho>=172",  "\\sigma_{TC sin(#chi)=#frac{1}{3}}", kRed));
    sigs.push_back(SignalSample("TCLO",                      "tc8TeVLO.dat",  "Rho:Xsec:percentError",  "Rho>=172",  "\\sigma_{TCLO sin(#chi)=#frac{1}{3}}", kOrange));
    //sigs.push_back(SignalSample("TC,sin(#chi)=#frac{1}{4}","xSec_TCWZ-sinchi1d4.dat",  "Rho:Xsec",  "Rho>=200",  "\\sigma_{TC sin(#chi)=#frac{1}{4}}", kRed, 3));
    //sigs.push_back(SignalSample("NoBkgObs",                 "../combined_limits/2013-05-20-NoBkg/nLimit_WprimeWZ_MarkovChainMC.txt",  "Mass:ObsLimit",  "",  "Obs_{NoBkg}", kBlue, 1));
    //sigs.push_back(SignalSample("NoBkg",                    "../combined_limits/2013-05-20-NoBkg/nLimit_WprimeWZ_MarkovChainMC.txt",  "Mass:ExpLimit",  "",  "Exp_{NoBkg}", kBlue, 2));
    //sigs.push_back(SignalSample("ZG50Percent",              "../combined_limits/2013-05-17-ZGamma50Percent/nLimit_WprimeWZ_MarkovChainMC.txt",  "Mass:ExpLimit",  "",  "Exp_{ZG50Percent}", kGreen, 2));
    //sigs.push_back(SignalSample("ZZ30Percent",              "../combined_limits/2013-05-28-ZZ30Percent/nLimit_WprimeWZ_MarkovChainMC.txt",  "Mass:ExpLimit",  "",  "Exp_{ZZ30Percent}", kGray, 2));
    //sigs.push_back(SignalSample("TT15DY30",                 "../combined_limits/2013-05-31-TT15DY30/nLimit_WprimeWZ_MarkovChainMC.txt",  "Mass:ExpLimit",  "",  "Exp_{TT15DY30}", kViolet, 2));
    //sigs.push_back(SignalSample("TT15DY50",                 "../combined_limits/2013-05-25-TT15DY50/nLimit_WprimeWZ_MarkovChainMC.txt",  "Mass:ExpLimit",  "",  "Exp_{TT15DY50}", kRed, 2));
    //sigs.push_back(SignalSample("TT15DY100",                "../combined_limits/2013-05-25-TT15DY100/nLimit_WprimeWZ_MarkovChainMC.txt",  "Mass:ExpLimit",  "",  "Exp_{TT15DY100}", kCyan, 2));
    //sigs.push_back(SignalSample("BinnedSys",                "../combined_limits/2013-05-28-BinnedSys/nLimit_WprimeWZ_MarkovChainMC.txt",  "Mass:ExpLimit",  "",  "Exp_{BinnedSys}", kOrange, 2));
    //sigs.push_back(SignalSample("LastApproval",                "../combined_limits/2013-05-10-Approval/nLimit_WprimeWZ_MarkovChainMC.txt",  "Mass:ExpLimit",  "",  "Exp_{Approval}", kBlue, 2));
    //sigs.push_back(SignalSample("5PercentMuPtUncert",         "../combined_limits/2013-06-12-5PercentMuPtScale/nLimit_WprimeWZ_MarkovChainMC.txt",  "Mass:ExpLimit",  "",  "Exp_{5PercentMuPtUncert}", kBlue, 2));
    //sigs.push_back(SignalSample("0p2PercentMuPtUncert",         "../combined_limits/2013-06-12-0p2PercentMuPtScale/nLimit_WprimeWZ_MarkovChainMC.txt",  "Mass:ExpLimit",  "",  "Exp_{0.2PercentMuPtUncert}", kRed, 2));

    if(outFile.empty()) outFile = "limitVsMass_WZ.pdf";
  }else if(inName.find("HadVZ") != string::npos){
    sigs.push_back(SignalSample("W'","xSec_WprimeVZ.dat", "Mass:Xsec", "", "\\sigma_{W'}"));
    sigs.push_back(SignalSample("RS","xSec_RSZZ.dat", "Mass:Xsec", "", "\\sigma_{RS}", kRed));
    if(outFile.empty()) outFile = "limitVsMass_VZ.pdf";
  }else if(inName.find("LNu") != string::npos){
    sigs.push_back(SignalSample("W'","xSec_SSMWprimeToLNu.dat", "Mass:Xsec", "", "\\sigma_{W'}"));
    if(outFile.empty()) outFile = "limitVsMass_LNu.pdf";
  }

  ////Load Trees////////
  TTree* tLimit = new TTree("tLimit", "Limits");
  tLimit->ReadFile(inFile.c_str());

  for(unsigned iSig=0; iSig<sigs.size(); ++iSig){
    sigs[iSig].tXsec = new TTree("tXsec", "Signal Cross Sections");
    sigs[iSig].tXsec->ReadFile(sigs[iSig].sigXsec.c_str());
  }
  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextSize(0.05);
  latexLabel.SetTextFont(42);

  //////Now Plot Limit vs Mass//////////////////////
  cout<<"Now Plot Limit vs Mass\n";
  TCanvas* c1 = new TCanvas("c1", "Exclusion Limit vs Mass");
  TMultiGraph *mg = new TMultiGraph("mg", ";M_{W', #rho_{TC}} (GeV);\\sigma #upoint BR (pb)");
  TLegend *leg = new TLegend(0.57, 0.45,0.85, 0.91,"");
  //TLegend *leg = new TLegend(0.6, 0.63,0.9, 0.89,"");
  c1->SetLogy();
  int n=0;
  float* x; 
  float* y;
  
  TGraph *gexp;
  TGraph *g1Sigma, *g2Sigma;
  TGraph* gData;

  //cout<<"Drawing Limit Curve\n";
  tLimit->Draw("Mass:Lumi:ObsLimit:ExpLimit:ExpLimitP1:ExpLimitM1:ExpLimitP2:ExpLimitM2", "", "para goff");
  n = tLimit->GetSelectedRows(); assert(n); 
  int treeIdx=0;
  const Double_t* mass = tLimit->GetVal(treeIdx++);
  const Double_t* lumi = tLimit->GetVal(treeIdx++);
  const Double_t* ObsLimit = tLimit->GetVal(treeIdx++);
  const Double_t* ExpLimit = tLimit->GetVal(treeIdx++);
  const Double_t* ExpLimitP1 = tLimit->GetVal(treeIdx++);
  const Double_t* ExpLimitM1 = tLimit->GetVal(treeIdx++);
  const Double_t* ExpLimitP2 = tLimit->GetVal(treeIdx++);
  const Double_t* ExpLimitM2 = tLimit->GetVal(treeIdx++);
  
  gexp = new TGraph(n, mass, ExpLimit);
  gexp->Sort();
  RemoveZeros(gexp);
  g2Sigma = new TGraph(2*n);
  g1Sigma = new TGraph(2*n);
  for (int i=0;i<n;i++) {
    g2Sigma->SetPoint(  i,mass[i]    ,ExpLimitP2[i]    );
    g2Sigma->SetPoint(n+i,mass[n-i-1],ExpLimitM2[n-i-1]);
    g1Sigma->SetPoint(  i,mass[i]    ,ExpLimitP1[i]    );
    g1Sigma->SetPoint(n+i,mass[n-i-1],ExpLimitM1[n-i-1]);
  }
  RemoveZeros(g2Sigma);
  RemoveZeros(g1Sigma);
  
  g2Sigma->SetFillColor(kYellow);
  mg->Add(g2Sigma, "F");//Cory
  g1Sigma->SetFillColor(kGreen);
  mg->Add(g1Sigma, "F");//Cory
  
  gexp->SetLineColor(kBlack);
  gexp->SetLineStyle(kDashed);//dashed
  gexp->SetLineWidth(3);
  mg->Add(gexp, "C");

  gData = new TGraph(n, mass, ObsLimit);
  gData->SetMarkerStyle(20);
  gData->SetLineWidth(3);
  mg->Add(gData, "CP");//Cory
  leg->AddEntry(gData,"Obs. 95% C.L.", "LP");//Cory

  leg->AddEntry(gexp,"Exp. 95% C.L.", "L");
  leg->AddEntry(g1Sigma,"Exp. #pm 1#sigma", "F");//Cory
  leg->AddEntry(g2Sigma,"Exp. #pm 2#sigma", "F");//Cory
   
  //////
  
  for(unsigned iSig=0; iSig<sigs.size(); ++iSig){
    bool drawBand = sigs[iSig].massString.find("percentError") != string::npos;
    sigs[iSig].tXsec->Draw( sigs[iSig].massString.c_str(), sigs[iSig].cutString.c_str(), "goff");
    n = sigs[iSig].tXsec->GetSelectedRows(); assert(n);
    x = new float[n]; y = new float[n];
    sigs[iSig].gXsecBand = new TGraph(2*n);
    float* sxSec = new float[n];
    for(int i=0; i<n; ++i){
      x[i] = sigs[iSig].tXsec->GetV1()[i];
      y[i] = sigs[iSig].tXsec->GetV2()[i];

      y[i] *= xSecScale;

      //Fill theory band around xsec curve
      if(drawBand){
        sxSec[i] = sigs[iSig].tXsec->GetV3()[i] / 100.;
        sigs[iSig].gXsecBand->SetPoint(    i  , x[i], y[i]*(1 + sxSec[i]));
        sigs[iSig].gXsecBand->SetPoint(2*n-i-1, x[i], y[i]*(1 - sxSec[i]));
      }
    }
    if(drawBand){
      sigs[iSig].gXsecBand->SetFillColor(sigs[iSig].bandColor);
      mg->Add(sigs[iSig].gXsecBand,"F");
    }

    sigs[iSig].gXsec = new TGraph(n, x, y);
    sigs[iSig].gXsec->Sort();
    sigs[iSig].gXsec->SetLineColor(sigs[iSig].lineColor);
    sigs[iSig].gXsec->SetLineStyle(sigs[iSig].lineStyle);
    if(drawBand) sigs[iSig].gXsec->SetFillColor(sigs[iSig].bandColor);
    mg->Add(sigs[iSig].gXsec,"C");
    if( 1 || sigs[iSig].lineStyle == 1 )//Only have legend entry for main band
      leg->AddEntry(sigs[iSig].gXsec,sigs[iSig].legendString.c_str(), drawBand ? "LF" : "L");

    
    delete[] x; delete[] y; delete sxSec;
  }

  cout<<"Drawing multigraph"<<endl;
  mg->SetMaximum(2.);
  mg->SetMinimum(0.00001);
  mg->Draw("a");
  //cout<<"ndiv is "<<mg->GetXaxis()->GetNdivisions()<<endl;
  mg->GetXaxis()->SetNdivisions(505);
  //cout<<"ndiv is "<<mg->GetXaxis()->GetNdivisions()<<endl;
  mg->Draw("a");

  for(unsigned iSig=0; iSig<sigs.size(); ++iSig){
    float lowObsLimit(-1), upObsLimit(-1);
    lowObsLimit = FindLimit(sigs[iSig].gXsec, gData, false); 
    upObsLimit  = FindLimit(sigs[iSig].gXsec, gData, true); 
    if(0){
      lowObsLimit = roundToNearest(lowObsLimit, 10);
      upObsLimit = roundToNearest(upObsLimit, 10);
    }
    cout<<sigs[iSig].name<<" Obs Limits are ["<<lowObsLimit<<", "<<upObsLimit<<"]"<<endl;
    if(upObsLimit > 0. && sigs[iSig].lineStyle == 1){
      //latexLabel.DrawLatex(0.20, 0.3-iSig*0.04, Form("#font[42]{Limit_{%s} = %.0f GeV}",sigs[iSig].name.c_str(),upObsLimit));
    }  
    float lowExpLimit(-1), upExpLimit(-1);
    lowExpLimit = FindLimit(sigs[iSig].gXsec, gexp, false); 
    upExpLimit  = FindLimit(sigs[iSig].gXsec, gexp, true); 
    if(0){
      lowExpLimit = roundToNearest(lowExpLimit, 10);
      upExpLimit = roundToNearest(upExpLimit, 10);
    }
    cout<<sigs[iSig].name<<" Exp Limits are ["<<lowExpLimit<<", "<<upExpLimit<<"]"<<endl;
  }

  latexLabel.DrawLatex(0.33, 0.96, "CMS Preliminary 2012");
  latexLabel.DrawLatex(0.19, 0.30, "#sqrt{s} = 8 TeV");
  latexLabel.DrawLatex(0.17, 0.20, Form("#intL dt = %.1f fb^{-1}",lumi[0]/1000.));
  
  //cout<<"Drawing legend"<<endl;
  leg->SetTextSize(0.05);
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  //leg->SetNColumns(2);
  //leg->SetColumnSeparation(0.05);
  leg->Draw();

  cout<<"Saving as "<<outFile<<endl;
  c1->SaveAs(outFile.c_str());
  string pngoutFile = outFile;
  pngoutFile.replace(pngoutFile.find(".pdf"), 4, ".png"); 
  c1->SaveAs(pngoutFile.c_str());

  return;
}

int
GetEndPoints(int & start, int& end, Double_t* y, float line){
  for(int i=start; i<end; ++i){
    if(y[i] >= line && y[i+1] >= line) start = i+1;
  }
  start = TMath::Min(start, end-1);
  end   = TMath::Min(start+2, end);
  return end-start + 1;
}

float 
FindLimit(const TGraph* xsec, const TGraph* limit, const bool findUpper){
/*  float lower = 300;///min(xsec->GetX()[0], limit->GetX()[0]);
  float upper = max(xsec->GetX()[xsec->GetN()-1], limit->GetX()[limit->GetN()-1]);
  const float INC = 1;
  for(float mass=lower; mass<upper; mass+=INC){
    if(xsec->Eval(mass) < limit->Eval(mass)) return mass;
  }
*/
  float upLimit = -1;
  float lower = 0;///min(xsec->GetX()[0], limit->GetX()[0]);
  float upper = max(xsec->GetX()[xsec->GetN()-1], limit->GetX()[limit->GetN()-1]);
  const float INC = 1;
  for(float mass=upper; mass>=lower; mass-=INC){
    if(xsec->Eval(mass) > limit->Eval(mass)){
      upLimit = mass+1;
      break;
    }
  }
  if(findUpper) return upLimit;
  for(float mass=upLimit-1; mass>=lower; mass-=INC){
    if(xsec->Eval(mass) < limit->Eval(mass)) return mass+1;
  }
  

  return -1;//should this be big now
}

             
void
RemoveZeros(TGraph* g){
  for(int i=0; i<g->GetN(); ++i){
    if(g->GetY()[i] < 1.e-9){
      cout<<"Removing point "<<i<<" with value "<<g->GetY()[i]<<endl;
      g->RemovePoint(i);
      i--;
    }
  }
}
