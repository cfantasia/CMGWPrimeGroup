#include "consts.h"
#include "TROOT.h"
#include "TTree.h"
#include "TGraph.h"
#include "TSystem.h"

void
ModXSec(string inFile, string outFile){
  TTree* tLimit = new TTree("tLimit", "Limits");
  tLimit->ReadFile(inFile.c_str());

  ofstream out(outFile.c_str());

  out<<"Mass/F:"
     <<"Xsec/F:"
     <<"percentError/F"
     <<endl;

  //kfactors
  TGraph* g = new TGraph(20);
  g->SetPoint(0,200,1.357);
  g->SetPoint(1,250,1.357);
  g->SetPoint(2,300,1.357);
  g->SetPoint(3,400,1.357);
  g->SetPoint(4,500,1.357);
  g->SetPoint(5,600,1.351);
  g->SetPoint(6,700,1.352);
  g->SetPoint(7,800,1.347);
  g->SetPoint(8,900,1.341);
  g->SetPoint(9,1000,1.332);
  g->SetPoint(10,1100,1.325);
  g->SetPoint(11,1200,1.311);
  g->SetPoint(12,1300,1.298);
  g->SetPoint(13,1400,1.279);
  g->SetPoint(14,1500,1.265);
  g->SetPoint(15,1600,1.255);
  g->SetPoint(16,1700,1.234);
  g->SetPoint(17,1800,1.211);
  g->SetPoint(18,1900,1.191);
  g->SetPoint(19,2000,1.184);

  TGraph* gPerErr = new TGraph(20);
  gPerErr->SetPoint(0,200,2.370);
  gPerErr->SetPoint(1,250,2.370);
  gPerErr->SetPoint(2,300,2.370);
  gPerErr->SetPoint(3,400,2.764);
  gPerErr->SetPoint(4,500,3.181);
  gPerErr->SetPoint(5,600,3.704);
  gPerErr->SetPoint(6,700,4.198);
  gPerErr->SetPoint(7,800,4.624);
  gPerErr->SetPoint(8,900,5.135);
  gPerErr->SetPoint(9,1000,5.695);
  gPerErr->SetPoint(10,1100,6.088);
  gPerErr->SetPoint(11,1200,6.516);
  gPerErr->SetPoint(12,1300,7.349);
  gPerErr->SetPoint(13,1400,7.760);
  gPerErr->SetPoint(14,1500,8.471);
  gPerErr->SetPoint(15,1600,8.471);
  gPerErr->SetPoint(16,1700,8.471);
  gPerErr->SetPoint(17,1800,8.471);
  gPerErr->SetPoint(18,1900,8.471);
  gPerErr->SetPoint(19,2000,8.471);
  
//  tLimit->Draw("Rho:Pi:Xsec:SinX", 
//  tLimit->Draw("Rho:Pi:Xsec", 
  tLimit->Draw("Mass:Xsec", 
               "", "para goff");
  float n = tLimit->GetSelectedRows(); 
  for(int isample=0; isample<n; ++isample){
    int idx=0;
    Double_t mass = tLimit->GetVal(idx++)[isample];
    Double_t xsec  = tLimit->GetVal(idx++)[isample];
    //Double_t perErr  = tLimit->GetVal(idx++)[isample];

    //cout<<"For rho: "<<rho<<" the kfactor is "<<g->Eval(rho)<<endl; 
    
    //xsec *= 1e9; //convert from mb to pb
    xsec *= 4./9.; //convert from emt to em
    xsec *= g->Eval(mass);//apply k factor
    //xsec *= (1.+0.75);//+75%
    //xsec *= (1.-0.40);//-40%%

    Double_t perErr  = gPerErr->Eval(mass);

    out.precision(0) ;
    out.setf ( ios::fixed, ios::floatfield);
    out  <<mass<<"\t";

    out.precision(4) ;
    out.setf ( ios::scientific, ios::floatfield);
    out  <<xsec<<"\t";
    
    out.precision(3);
    out.setf ( ios::fixed, ios::floatfield);
    out<<perErr<<"\t";
    out<<endl;

  }

}
