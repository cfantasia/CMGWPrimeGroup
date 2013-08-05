//Usage: root -b -q 'OptimizeWindows.C+()'
#include "consts.h"

void
OptimizeWindows(){
//Read in limits file and find best one for each mass
  

  TTree* tCodes = new TTree("tCodes", "Codes");
  tCodes->ReadFile("cutValues.wz.dat");
  tCodes->Draw("SignalCode", "", "para goff");
  float nCodes = tCodes->GetSelectedRows(); 
  const Double_t* codes = tCodes->GetVal(0);

  string outfile("cutValues.wz.dat.new");
  ofstream out(outfile.c_str());
  if(!out) { 
    cout << "Cannot open file " << outfile << endl; 
    abort();
  } 

  out<<"SignalCode/F:Mass/F:nGen/F:bkgSysErr/F:HtCut/F:ZptCut/F:WptCut/F:minWindow/F:maxWindow/F:ZJets/F:sZJets/F:WZKFactor/F"<<endl;
  out  << setiosflags(ios::fixed) << setprecision(0);

  TTree* tLims = new TTree("tLims", "Lims");
  tLims->ReadFile("nLimit.txt");

  TTree* tCuts = new TTree("tCuts", "Cuts");
  tCuts->ReadFile("cutValues.wzFull.dat");

  for(int icode=0; icode<nCodes; ++icode){
    tLims->Draw("ExpLimit", Form("SignalCode==%f", codes[icode]), "para goff");
    float nLims = tLims->GetSelectedRows(); 
    int bestIdx = 0;
    for(int iLims=0; iLims<nLims; ++iLims){
      const Double_t expLim     = tLims->GetVal(0)[iLims];
      if(expLim < tLims->GetVal(0)[bestIdx]) bestIdx = iLims;
    }

    if(bestIdx == 0) cout<<" For code "<<codes[icode]<<" the best window was the smaller, consider trying even smaller\n";
    if(bestIdx == nLims-1) cout<<" For code "<<codes[icode]<<" the best window was the biggest, consider trying even bigger\n";

    tCuts->Draw("SignalCode:Mass:nGen:bkgSysErr:HtCut:ZptCut:WptCut:minWindow:maxWindow:ZJets:sZJets:WZKFactor",
                Form("SignalCode==%f", codes[icode]), "para goff");
    int iVar=0;
    const Double_t  signalCode = tCuts->GetVal(iVar++)[bestIdx];
    const Double_t  mass      = tCuts->GetVal(iVar++)[bestIdx];
    const Double_t  nGen      = tCuts->GetVal(iVar++)[bestIdx];
    const Double_t  bkgSysErr = tCuts->GetVal(iVar++)[bestIdx];
    const Double_t  HtCut     = tCuts->GetVal(iVar++)[bestIdx];
    const Double_t  ZptCut    = tCuts->GetVal(iVar++)[bestIdx];
    const Double_t  WptCut    = tCuts->GetVal(iVar++)[bestIdx];
    const Double_t  minWind   = tCuts->GetVal(iVar++)[bestIdx];
    const Double_t  maxWind   = tCuts->GetVal(iVar++)[bestIdx];
    const Double_t  ZJets     = tCuts->GetVal(iVar++)[bestIdx];
    const Double_t  sZJets    = tCuts->GetVal(iVar++)[bestIdx];
    const Double_t  WZKFactor = tCuts->GetVal(iVar++)[bestIdx];

    out<<setprecision(1)
       <<signalCode<<" "
       <<setprecision(0)
       <<mass<<" "
       <<nGen<<" "
       <<setprecision(3)
       <<bkgSysErr<<" "
       <<setprecision(0)
       <<HtCut<<" "
       <<ZptCut<<" "
       <<WptCut<<" "
       <<minWind<<" "
       <<maxWind<<" "
       <<setprecision(2)
       <<ZJets<<" "
       <<sZJets<<" "
       <<setprecision(1)
       <<WZKFactor<<" "
       <<endl;
  }
}
