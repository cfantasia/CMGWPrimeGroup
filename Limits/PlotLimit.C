//Usage: root -l PlotLimit.C++

#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "../root_macros/CMSStyle.C"
#include "consts.h"

int GetEndPoints(int & start, int& end, Double_t* y, float line);
float FindLimit(const TGraph* xsec, const TGraph* limit);
float Slope(const float x1, const float y1,
            const float x2, const float y2);
float Intercept(const float x1, const float y1,
                const float x2, const float y2);
float IntersectionX(const float s1, const float i1,
                    const float s2, const float i2);
void
PlotLimit(string inName){
  gErrorIgnoreLevel = kWarning;
  CMSstyle();

  string sigXsec, legendString, massString, cutString, outFile;
  if(inName.find("WprimeWZ") != string::npos){
    sigXsec = "xSec_WZ.dat";
    legendString = "\\sigma_{W'}";
    massString = "Mass:Xsec";
    outFile = "limitVsMass_WZ.pdf";
  }else if(inName.find("TCWZ") != string::npos){
    sigXsec = "xSec_TCWZ.dat";
    legendString = "\\sigma_{TC}";
    massString = "Rho:Xsec";
    cutString = "Rho>=300 && Pi==Rho*3/4 - 25";
    outFile = "limitVsMass_TCWZ.pdf";
  }else if(inName.find("HadVZ") != string::npos){
    cout<<"Using HadVZ\n";
    sigXsec = "xSec_WprimeVZ.dat";
    legendString = "\\sigma_{W'}";
//    sigXsec = "xSec_RSZZ.dat";
//    legendString = "\\sigma_{RS}";
    massString = "Mass:Xsec";
    outFile = "limitVsMass_VZ.pdf";
  }

  ////Load Trees////////
  TTree* tLimit = new TTree("tLimit", "Limits");
  tLimit->ReadFile("nLimit.txt");

  TTree* tXsec = new TTree("tXsec", "Signal Cross Sections");
  tXsec->ReadFile(sigXsec.c_str());

  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextSize(0.045);

  //////Now Plot Limit vs Mass//////////////////////
  cout<<"Now Plot Limit vs Mass\n";
  TCanvas* c1 = new TCanvas("c1", "Exclusion Limit vs Mass");
  TMultiGraph *mg = new TMultiGraph("mg", ";M_{WZ} (GeV);\\sigma #upoint BR (pb)");
  TLegend *legend = new TLegend(0.18,0.18,0.52,0.42,"");
  c1->SetLogy();
  int n=0;
  float* x; 
  float* y;

  TGraph *glumi;
  TGraph *g1Sigma, *g2Sigma;
  TGraph* gData;

  //cout<<"Drawing Limit Curve\n";
  tLimit->Draw("SignalCode:Mass:Lumi:ObsLimit:ExpLimit:ExpLimitP1:ExpLimitM1:ExpLimitP2:ExpLimitM2", cutString.c_str(), "para goff");
  n = tLimit->GetSelectedRows(); assert(n); 
  
  const Double_t* mass = tLimit->GetVal(1);
  const Double_t* lumi = tLimit->GetVal(2);
  const Double_t* ObsLimit = tLimit->GetVal(3);
  const Double_t* ExpLimit = tLimit->GetVal(4);
  const Double_t* ExpLimitP1 = tLimit->GetVal(5);
  const Double_t* ExpLimitM1 = tLimit->GetVal(6);
  const Double_t* ExpLimitP2 = tLimit->GetVal(7);
  const Double_t* ExpLimitM2 = tLimit->GetVal(8);
  
  glumi = new TGraph(n, mass, ExpLimit);
  g2Sigma = new TGraph(2*n);
  g1Sigma = new TGraph(2*n);
  for (int i=0;i<n;i++) {
    g2Sigma->SetPoint(  i,mass[i]    ,ExpLimitP2[i]    );
    g2Sigma->SetPoint(n+i,mass[n-i-1],ExpLimitM2[n-i-1]);
    g1Sigma->SetPoint(  i,mass[i]    ,ExpLimitP1[i]    );
    g1Sigma->SetPoint(n+i,mass[n-i-1],ExpLimitM1[n-i-1]);
  }
  
//      g2Sigma->SetFillStyle(4004);
  g2Sigma->SetFillColor(kGreen);
  mg->Add(g2Sigma, "F");
//      g1Sigma->SetFillStyle(4004);
  g1Sigma->SetFillColor(kYellow);
  mg->Add(g1Sigma, "F");
  
  gData = new TGraph(n, mass, ObsLimit);
  gData->SetMarkerStyle(20);
  mg->Add(gData, "P");
  legend->AddEntry(gData,"Obs. Limit", "P");
  
  
  glumi->SetLineColor(kBlack);   
  mg->Add(glumi, "L");
  legend->AddEntry(glumi,"Exp. Limit", "L");
  legend->AddEntry(g1Sigma,"\\pm 1\\sigma", "F");
  legend->AddEntry(g2Sigma,"\\pm 2\\sigma", "F");
  
  
  float limit(-1);
  
  tXsec->Draw(massString.c_str(), cutString.c_str(), "goff");
  n = tXsec->GetSelectedRows(); assert(n);
  x = new float[n]; y = new float[n];
  for(int i=0; i<n; ++i){
    x[i] = tXsec->GetV1()[i];
    y[i] = tXsec->GetV2()[i];
  }
  TGraph* gxsec = new TGraph(n, x, y);
  mg->Add(gxsec,"l");
  legend->AddEntry(gxsec,legendString.c_str(), "L");
  delete x; delete y;

  limit = FindLimit(gxsec, gData); cout<<"Limit is "<<limit<<endl;
   

  //cout<<"Drawing multigraph"<<endl;
  mg->SetMinimum(0.001);
  mg->Draw("a");

  if(limit > 0.) 
    latexLabel.DrawLatex(0.69, 0.75, Form("#font[132]{Limit = %.0f GeV}",limit));
  
  latexLabel.DrawLatex(0.18, 0.96, "#font[132]{CMS Preliminary 2011}");
  latexLabel.DrawLatex(0.5, 0.955, "#font[132]{#sqrt{s} = 7 TeV}");
  latexLabel.DrawLatex(0.73, 0.87, Form("#font[132]{#intL dt = %.2f fb^{-1}}",lumi[0]/1000.));
  
  //cout<<"Drawing legend"<<endl;
  legend->SetTextSize(0.05);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  //legend->SetNColumns(2);
  legend->Draw();

  c1->SaveAs(outFile.c_str());

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
FindLimit(const TGraph* xsec, const TGraph* limit){
  int Nxsec = xsec->GetN();
  Double_t* Xxsec = xsec->GetX(); 
  Double_t* Yxsec = xsec->GetY(); 

  int Nlimit = limit->GetN();
  Double_t* Xlimit = limit->GetX(); 
  Double_t* Ylimit = limit->GetY(); 
  //assert(Nxsec == Nlimit);

  int idx1(-1), idx2(-1);
  int Npoints = min(Nxsec,Nlimit);
  for(int i=0; i<Npoints; ++i){
    if(Ylimit[i] > Yxsec[i]){
      idx2=i; 
      break;
    }
  }
  if(idx2 < 1) return -1;
  idx1 = max(0, idx2 -1);
  
  float s1(0.), i1(0.), s2(0.), i2(0.);
  s1 = Slope(Xxsec[idx1], Yxsec[idx1],
             Xxsec[idx2], Yxsec[idx2]);
  i1 = Intercept(Xxsec[idx1], Yxsec[idx1],
                 Xxsec[idx2], Yxsec[idx2]);
  s2 = Slope(Xlimit[idx1], Ylimit[idx1],
             Xlimit[idx2], Ylimit[idx2]);
  i2 = Intercept(Xlimit[idx1], Ylimit[idx1],
                 Xlimit[idx2], Ylimit[idx2]);
  return IntersectionX(s1, i1, s2, i2);
}

float Slope(const float x1, const float y1,
            const float x2, const float y2){
  return (y2-y1)/(x2-x1);
}
float Intercept(const float x1, const float y1,
                const float x2, const float y2){
  return y1 - x1*(y2-y1)/(x2-x1);
}

float
IntersectionX(const float s1, const float i1,
              const float s2, const float i2){
  return (i2-i1)/(s1-s2);
}
             