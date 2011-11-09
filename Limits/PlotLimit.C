//Usage: root -l PlotLimit.C++

#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "../root_macros/CMSStyle.C"
#include "consts.h"

int GetEndPoints(int & start, int& end, Double_t* y, float line);
float FindLimit(const TGraph* xsec, const TGraph* limit);

struct SignalSample{
  string name;
  string sigXsec;
  string massString;
  string cutString;
  string legendString;
  TTree* tXsec;
  int lineColor;
  TGraph* gXsec;

  SignalSample(){}
  SignalSample(string n, string xsec, string mass, string cut, string leg, int lc=kBlack){
    name = n;
    sigXsec = xsec;
    massString = mass;
    cutString = cut;
    legendString = leg;
    tXsec = NULL;
    lineColor = lc;
    gXsec = NULL;
  }
};

void
PlotLimit(string inName){
  gErrorIgnoreLevel = kWarning;
  CMSstyle();

  string outFile;
  vector<SignalSample> sigs;
  if(inName.find("WprimeWZ") != string::npos){
    sigs.push_back(SignalSample("W'","xSec_WZ.dat",  "Mass:Xsec",  "", "\\sigma_{W'}"));
    sigs.push_back(SignalSample("TC","xSec_TCWZ.dat",  "Rho:Xsec",  "Rho>=300 && Pi==Rho*3/4 - 25",  "\\sigma_{TC}", kRed));
    outFile = "limitVsMass_WZ.pdf";
  }else if(inName.find("HadVZ") != string::npos){
    sigs.push_back(SignalSample("W'","xSec_WprimeVZ.dat", "Mass:Xsec", "", "\\sigma_{W'}"));
    sigs.push_back(SignalSample("RS","xSec_RSZZ.dat", "Mass:Xsec", "", "\\sigma_{RS}", kRed));
    outFile = "limitVsMass_VZ.pdf";
  }

  ////Load Trees////////
  TTree* tLimit = new TTree("tLimit", "Limits");
  tLimit->ReadFile("nLimit.txt");

  for(unsigned iSig=0; iSig<sigs.size(); ++iSig){
    sigs[iSig].tXsec = new TTree("tXsec", "Signal Cross Sections");
    sigs[iSig].tXsec->ReadFile(sigs[iSig].sigXsec.c_str());
  }
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
  ///////////tLimit->Draw("SignalCode:Mass:Lumi:ObsLimit:ExpLimit:ExpLimitP1:ExpLimitM1:ExpLimitP2:ExpLimitM2", cutString.c_str(), "para goff");
  tLimit->Draw("SignalCode:Mass:Lumi:ObsLimit:ExpLimit:ExpLimitP1:ExpLimitM1:ExpLimitP2:ExpLimitM2", "", "para goff");
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
  
  glumi->SetLineColor(kBlack);
  glumi->SetLineStyle(2);//dashed
  mg->Add(glumi, "L");
  legend->AddEntry(glumi,"Exp. Limit", "L");
  legend->AddEntry(g1Sigma,"Exp. \\pm 1\\sigma", "F");
  legend->AddEntry(g2Sigma,"Exp. \\pm 2\\sigma", "F");

  gData = new TGraph(n, mass, ObsLimit);
  gData->SetMarkerStyle(20);
  mg->Add(gData, "P");
  legend->AddEntry(gData,"Obs. Limit", "P");
   
  //////
  
  for(unsigned iSig=0; iSig<sigs.size(); ++iSig){
    sigs[iSig].tXsec->Draw( sigs[iSig].massString.c_str(), sigs[iSig].cutString.c_str(), "goff");
    n = sigs[iSig].tXsec->GetSelectedRows(); assert(n);
    x = new float[n]; y = new float[n];
    for(int i=0; i<n; ++i){
      x[i] = sigs[iSig].tXsec->GetV1()[i];
      y[i] = sigs[iSig].tXsec->GetV2()[i];
    }
    sigs[iSig].gXsec = new TGraph(n, x, y);
    sigs[iSig].gXsec->SetLineColor(sigs[iSig].lineColor);
    mg->Add(sigs[iSig].gXsec,"l");
    legend->AddEntry(sigs[iSig].gXsec,sigs[iSig].legendString.c_str(), "L");
    delete x; delete y;
  }

  cout<<"Drawing multigraph"<<endl;
  mg->SetMinimum(0.0001);
  mg->Draw("a");

  for(unsigned iSig=0; iSig<sigs.size(); ++iSig){
    float limit(-1);
    limit = FindLimit(sigs[iSig].gXsec, gData); 
    cout<<sigs[iSig].name<<" Limit is "<<limit<<endl;
    if(limit > 0.) 
      latexLabel.DrawLatex(0.69, 0.75-iSig*0.05, Form("#font[132]{Limit_{%s} = %.0f GeV}",sigs[iSig].name.c_str(),limit));
  }  
  latexLabel.DrawLatex(0.18, 0.96, "#font[132]{CMS Preliminary 2011}");
  latexLabel.DrawLatex(0.5, 0.955, "#font[132]{#sqrt{s} = 7 TeV}");
  latexLabel.DrawLatex(0.73, 0.87, Form("#font[132]{#intL dt = %.2f fb^{-1}}",lumi[0]/1000.));
  
  //cout<<"Drawing legend"<<endl;
  legend->SetTextSize(0.05);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  //legend->SetNColumns(2);
  legend->Draw();

  cout<<"Saving as "<<outFile<<endl;
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
  float lower = min(xsec->GetX()[0], limit->GetX()[0]);
  float upper = max(xsec->GetX()[xsec->GetN()-1], limit->GetX()[limit->GetN()-1]);
  const float INC = 1;
  for(float mass=lower; mass<upper; mass+=INC){
    if(xsec->Eval(mass) < limit->Eval(mass)) return mass;
  }

  return -1;
}

             
