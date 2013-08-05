#include "consts.h"
#include "TROOT.h"
#include "TTree.h"
#include "TGraph.h"
#include "TSystem.h"

void
AddBand(string inFile, string outFile){
  TTree* tLimit = new TTree("tLimit", "Limits");
  tLimit->ReadFile(inFile.c_str());

  ofstream out(outFile.c_str());

  out<<"Rho/F:"
     <<"Xsec/F:"
     <<"percentError/F"
     <<endl;

  //kfactors
  TGraph* g = new TGraph(15);
  g->SetPoint(0,200,2.370);
  g->SetPoint(1,250,2.370);
  g->SetPoint(2,300,2.370);
  g->SetPoint(3,400,2.764);
  g->SetPoint(4,500,3.181);
  g->SetPoint(5,600,3.704);
  g->SetPoint(6,700,4.198);
  g->SetPoint(7,800,4.624);
  g->SetPoint(8,900,5.135);
  g->SetPoint(9,1000,5.695);
  g->SetPoint(10,1100,6.088);
  g->SetPoint(11,1200,6.516);
  g->SetPoint(12,1300,7.349);
  g->SetPoint(13,1400,7.760);
  g->SetPoint(14,1500,8.471);
  
  tLimit->Draw("Rho:Xsec", 
               "", "para goff");
  float n = tLimit->GetSelectedRows(); 
  for(int isample=0; isample<n; ++isample){
    int idx=0;
    Double_t rho = tLimit->GetVal(idx++)[isample];
    Double_t xsec  = tLimit->GetVal(idx++)[isample];

    //cout<<"For rho: "<<rho<<" the kfactor is "<<g->Eval(rho)<<endl; 
    
    //xsec *= 1e9; //convert from mb to pb
    //xsec *= 4./9.; //convert from emt to em
    double percentError = g->Eval(rho);
    //xsec *= (1.+0.75);//+75%
    //xsec *= (1.-0.40);//-40%%

    out.precision(0) ;
    out.setf ( ios::fixed, ios::floatfield);
    out  <<rho<<"\t";

    out.precision(4) ;
    out.setf ( ios::scientific, ios::floatfield);
    out  <<xsec<<"\t";
    out  <<percentError<<"\t"
         <<endl;

  }

}
