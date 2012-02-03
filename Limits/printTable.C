/* 
Author: Cory Fantasia
Script to print important numbers in LaTex form
To Use: 
root -b -q printTable.C+
*/
#include "consts.h"
#include "TROOT.h"
#include "TSystem.h"

int FindIdx(Double_t* arr, int size, double num){
  for(int i=0; i<size; ++i)
    if(arr[i] == num) return i;
  
  cout<<"Didn't find number: "<<num<<endl;
  abort();
}

void
printTable(int mode){
  //Need 3 files open:
  //Cuts File
  TTree* tCuts = new TTree("tCuts", "Cuts");
  tCuts->ReadFile("cutValues.wz.dat");
  tCuts->Draw("SignalCode:Mass:minWindow:maxWindow:HtCut", "", "para goff");
  float nCuts = tCuts->GetSelectedRows(); 
  
  //# of Events File
  TTree* tEvts = new TTree("tEvts", "Evts");
  tEvts->ReadFile("nEvents.txt");
  tEvts->Draw("SignalCode:DataEvts:MCEvts:statMCEvts:BkgEvts:statBkgEvts:Eff:sEff", 
             "", "para goff");
  float nEvts = tEvts->GetSelectedRows(); 

  //Limits File
  TTree* tLims = new TTree("tLims", "Limits");
  tLims->ReadFile("nLimit.txt");
  tLims->Draw("SignalCode:ExpLimit:ObsLimit", "", "para goff");
  float nLims = tLims->GetSelectedRows(); 
  
  switch(mode){

  case 0:
    cout<<"Mass point Window NBkgMC eSig NSig Data Exp. Limit Obs. Limit"<<endl;

    //loop over signals:
    for(int iLims=0; iLims<nLims; ++iLims){
      const Double_t  signalCode = tLims->GetVal(0)[iLims];
      int iCuts = FindIdx(tCuts->GetVal(0), nCuts, signalCode);
      int iEvts = FindIdx(tEvts->GetVal(0), nEvts, signalCode);

      const Double_t  mass = tCuts->GetVal(1)[iCuts];
      const Double_t  minWind = tCuts->GetVal(2)[iCuts];
      const Double_t  maxWind = tCuts->GetVal(3)[iCuts];

      const Double_t  DataEvts = tEvts->GetVal(1)[iCuts];
      const Double_t     MCEvts = tEvts->GetVal(2)[iCuts];
      const Double_t statMCEvts = tEvts->GetVal(3)[iCuts];
      const Double_t     BkgEvts = tEvts->GetVal(4)[iCuts];
      const Double_t statBkgEvts = tEvts->GetVal(5)[iCuts];
      const Double_t     SigEvts = MCEvts - BkgEvts;
      const Double_t statSigEvts = sqrt(statMCEvts*statMCEvts - statBkgEvts*statBkgEvts);
      const Double_t       Eff = tEvts->GetVal(6)[iCuts];
      const Double_t      sEff = tEvts->GetVal(7)[iCuts];

      const Double_t expLim = tLims->GetVal(1)[iLims];
      const Double_t obsLim = tLims->GetVal(2)[iLims];

      cout<<setiosflags(ios::fixed) << setprecision(0)
          <<"W' "<<mass<<" & "
          <<minWind<<"-"<<maxWind<<" & "
          <<std::fixed << std::setprecision(1)
          <<Value(BkgEvts,statBkgEvts)<<" $\\pm$ "<<Value(statBkgEvts)<<" & "
        //<<BkgEvts<<" $\\pm$ "<<statBkgEvts<<" & "
          <<(int)DataEvts<<" & "
          <<std::fixed << std::setprecision(1)
        //<<SigEvts<<" $\\pm$ "<<statSigEvts<<" & "
          <<Value(SigEvts,statSigEvts)<<" $\\pm$ "<<Value(statSigEvts)<<" & "
          <<std::fixed << std::setprecision(1)
        //<<Eff*100<<" $\\pm$ "<<sEff*100<<" & "
          <<Value(Eff*100,sEff*100)<<" $\\pm$ "<<Value(sEff*100)<<" & "
          <<setprecision(4)
          <<Value(expLim, -2)<<" & "
          <<Value(obsLim, -2)
          <<" \\\\ \\hline"
          <<endl;

    }//Table 0

    break;

  case 1:
    cout<<"\n\n ----After Ht Cut Table (remember to use noWind option!!!!---------------\n\n";
  
    //loop over signals:
    for(int iEvts=0; iEvts<nEvts; ++iEvts){
      const Double_t  signalCode = tEvts->GetVal(0)[iEvts];
      int iCuts = FindIdx(tCuts->GetVal(0), nCuts, signalCode);

      const Double_t  mass = tCuts->GetVal(1)[iCuts];
      const Double_t  minWind = tCuts->GetVal(2)[iCuts];
      const Double_t  maxWind = tCuts->GetVal(3)[iCuts];
      const Double_t  minHt   = tCuts->GetVal(4)[iCuts];
      const Double_t  DataEvts = tEvts->GetVal(1)[iCuts];
      const Double_t     MCEvts = tEvts->GetVal(2)[iCuts];
      const Double_t statMCEvts = tEvts->GetVal(3)[iCuts];
      const Double_t     BkgEvts = tEvts->GetVal(4)[iCuts];
      const Double_t statBkgEvts = tEvts->GetVal(5)[iCuts];
      const Double_t     SigEvts = MCEvts - BkgEvts;
      const Double_t statSigEvts = sqrt(statMCEvts*statMCEvts - statBkgEvts*statBkgEvts);

      cout<<setiosflags(ios::fixed) << setprecision(0)
          <<"W' "<<mass<<" & "
          <<minHt<<" & "
          <<setiosflags(ios::fixed) << setprecision(1)
        //<<BkgEvts<<" $\\pm$ "<<statBkgEvts<<" & "
          <<Value(BkgEvts,statBkgEvts)<<" $\\pm$ "<<Value(statBkgEvts)<<" & "
          <<(int)DataEvts<<" & "
          <<setiosflags(ios::fixed) << setprecision(1)
        //<<SigEvts<<" $\\pm$ "<<statSigEvts
          <<Value(SigEvts,statSigEvts)<<" $\\pm$ "<<Value(statSigEvts)
          <<" \\\\ \\hline"
          <<endl;
    }//Table 1
    break;
    
  default:
    cout<<" What mode did you want???"<<endl;
  }
}
