//Author: Cory Fantasia 2012
//Purpose: Modify input file cross sections to account for tau or k-factors
//Usage: root -b -l -q 'ModXSec.C+(string inFile, string outFile, bool removeTaus=true, bool applyKFactors=true)'

#include "consts.h"
#include "TROOT.h"
#include "TTree.h"
#include "TGraph.h"
#include "TSystem.h"

void
ModXSec(string inFile, string outFile, bool removeTaus=true, bool applyKFactors=true){
  bool doTC = inFile.find("TC") != string::npos;
  cout<<"doTC is "<<doTC<<endl;

  TTree* tLimit = new TTree("tLimit", "Limits");
  tLimit->ReadFile(inFile.c_str());

  ofstream out(outFile.c_str());

  if(doTC)  out<<"Rho/F:";
  else      out<<"Mass/F:";
  if(doTC){
    out<<"Pi/F:";
    out<<"SinX/F:";
  }
  out<<"Xsec/F:"
     <<"percentError/F"
     <<endl;

  //kfactors
  TGraph* gK = new TGraph(20);
  gK->SetPoint( 0, 200,1.347);
  gK->SetPoint( 1, 300,1.347);
  gK->SetPoint( 2, 500,1.363);
  gK->SetPoint( 3, 700,1.351);
  gK->SetPoint( 4, 900,1.347);
  gK->SetPoint( 5,1100,1.331);
  gK->SetPoint( 6,1300,1.317);
  gK->SetPoint( 7,1500,1.293);
  gK->SetPoint( 8,1700,1.257);
  gK->SetPoint( 9,1900,1.230);
  gK->SetPoint( 10,2000,1.214);

  //Signal pdf uncertainties
  TGraph* gPerErr = new TGraph(20);
  gPerErr->SetPoint(  0, 200,2.192);
  gPerErr->SetPoint(  1, 300,2.192);
  gPerErr->SetPoint(  2, 500,2.632);
  gPerErr->SetPoint(  3, 700,3.070);
  gPerErr->SetPoint(  4, 900,3.448);
  gPerErr->SetPoint(  5,1100,3.771);
  gPerErr->SetPoint(  6,1300,4.101);
  gPerErr->SetPoint(  7,1500,4.339);
  gPerErr->SetPoint(  8,1700,4.581);
  gPerErr->SetPoint(  9,1900,4.846);
  gPerErr->SetPoint( 10,2000,4.981);
  
  if(doTC) tLimit->Draw("Rho:Xsec:Pi:SinX", "SinX > 0.32 && SinX<0.34", "para goff");
  //if(doTC) tLimit->Draw("Rho:Xsec:Pi:SinX", "(280<=Rho && Rho<=300) && (150<=Pi && Pi<=160)", "para goff");
  else     tLimit->Draw("Mass:Xsec",        "", "para goff");
  float n = tLimit->GetSelectedRows(); 
  for(int isample=0; isample<n; ++isample){
    int idx=0;
    Double_t mass = tLimit->GetVal(idx++)[isample];
    Double_t xsec  = tLimit->GetVal(idx++)[isample];

    Double_t pi(-1), sinx(-1);
    if(doTC){
      pi   = tLimit->GetVal(idx++)[isample];
      sinx = tLimit->GetVal(idx++)[isample];
    }

    //cout<<"For rho: "<<rho<<" the kfactor is "<<gK->Eval(rho)<<endl; 
    
    //xsec *= 1e9; //convert from mb to pb
    if(removeTaus) xsec *= 4./9.; //convert from emt to em
    if(applyKFactors) xsec *= gK->Eval(mass);//apply k factor

    Double_t perErr  = gPerErr->Eval(mass);

    out.precision(0) ;
    out.setf ( ios::fixed, ios::floatfield);
    out  <<mass<<"\t";

    if(doTC){
      out.precision(0) ;
      out.setf ( ios::fixed, ios::floatfield);
      out  <<pi<<"\t";

      out.precision(4) ;
      out.setf ( ios::fixed, ios::floatfield);
      out  <<sinx<<"\t";
    }

    out.precision(4) ;
    out.setf ( ios::scientific, ios::floatfield);
    out  <<xsec<<"\t";
    
    out.precision(3);
    out.setf ( ios::fixed, ios::floatfield);
    out<<perErr<<"\t";
    out<<endl;

  }

}
