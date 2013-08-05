//Usage: root -b -q 'OptimizeWindows.C+()'
#include "consts.h"

void
MakeAllWindows(){

  TTree* tCuts = new TTree("tCuts", "Cuts");
  tCuts->ReadFile("cutValues.wz.dat");
  tCuts->Draw("SignalCode:Mass:nGen:bkgSysErr:HtCut:ZptCut:WptCut:minWindow:maxWindow:ZJets:sZJets:WZKFactor", "", "para goff");
  float nCuts = tCuts->GetSelectedRows(); 

  string outfile("cutValues.wzFull.dat");
  ofstream out(outfile.c_str());
  if(!out) { 
    cout << "Cannot open file " << outfile << endl; 
    abort();
  } 

  out<<"SignalCode/F:Mass/F:nGen/F:bkgSysErr/F:HtCut/F:ZptCut/F:WptCut/F:minWindow/F:maxWindow/F:ZJets/F:sZJets/F:WZKFactor/F"<<endl;
  out  << setiosflags(ios::fixed) << setprecision(0);

  //loop over signals:
  for(int iCuts=0; iCuts<nCuts; ++iCuts){
    int iVar=0;
    const Double_t  signalCode = tCuts->GetVal(iVar++)[iCuts];
    const Double_t  mass = tCuts->GetVal(iVar++)[iCuts];
    const Double_t  nGen = tCuts->GetVal(iVar++)[iCuts];
    const Double_t  bkgSysErr = tCuts->GetVal(iVar++)[iCuts];
    const Double_t  HtCut = tCuts->GetVal(iVar++)[iCuts];
    const Double_t  ZptCut = tCuts->GetVal(iVar++)[iCuts];
    const Double_t  WptCut = tCuts->GetVal(iVar++)[iCuts];
    const Double_t  minWind = tCuts->GetVal(iVar++)[iCuts];
    const Double_t  maxWind = tCuts->GetVal(iVar++)[iCuts];
    const Double_t  ZJets = tCuts->GetVal(iVar++)[iCuts];
    const Double_t  sZJets = tCuts->GetVal(iVar++)[iCuts];
    const Double_t  WZKFactor = tCuts->GetVal(iVar++)[iCuts];

    out  << setiosflags(ios::fixed) << setprecision(0);
    out<<"# Window Optmization for mass "<<mass<<" GeV"<<endl;
    
    const int htInc = 10;
    int minHtOffset(-70), maxHtOffset(70);
    for(int minHt=HtCut+minHtOffset; minHt<HtCut+maxHtOffset; minHt+=htInc){
      if(minHt < 0){
        maxHtOffset += htInc;
        continue;
      }
      const int windInc = 10;
      int minWindOffset(-6*windInc), maxWindOffset(6*windInc);
      for(int windOffset=minWindOffset; windOffset<=maxWindOffset; windOffset+=windInc){
        float lowerW = minWind - windOffset;
        float upperW = maxWind + windOffset;
        
        if(lowerW >= mass || upperW <= mass){
          maxWindOffset += windInc;//if we can't go narrower, let's try wider
          continue;
        }
        
        out<<setprecision(1)
           <<signalCode<<" "
           <<setprecision(0)
           <<mass<<" "
           <<nGen<<" "
           <<setprecision(3)
           <<bkgSysErr<<" "
           <<setprecision(0)
           <<minHt<<" "
           <<ZptCut<<" "
           <<WptCut<<" "
           <<lowerW<<" "//low mass window
           <<upperW<<" "//high mass window
           <<setprecision(2)
           <<ZJets<<" "
           <<sZJets<<" "
           <<setprecision(1)
           <<WZKFactor<<" "
           <<endl;
      }
    }
  }
}
