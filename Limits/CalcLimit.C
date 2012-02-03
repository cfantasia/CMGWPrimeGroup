/* 
Author: Cory Fantasia
Wrapper to calculate limits for each mass point
To Use: 
root -b -q -n CalcLimit.C+
*/
#include "consts.h"
#include "TROOT.h"
#include "TSystem.h"
#include "../../../StatisticalTools/RooStatsRoutines/root/roostats_cl95.C"

void
CalcLimit(bool useCLs=true, string inName="nEvents.txt", string outName="nLimit.txt"){
  gErrorIgnoreLevel = kWarning;
//  gSystem->SetIncludePath( "-I$ROOFITSYS/include" );
  gSystem->SetIncludePath( "-I/afs/hep.wisc.edu/cern/.root/root_v5.30.00.Linux-slc5_amd64-gcc4.3/include/RooStats" );
//  gROOT->ProcessLine(".L ../../../StatisticalTools/RooStatsRoutines/root/roostats_cl95.C+");
  
  string outfile(outName.c_str());
  ofstream out(outfile.c_str());
  if(!out) { 
    cout << "Cannot open file " << outfile << endl; 
    abort();
  } 
  
  out  << setiosflags(ios::fixed) << setprecision(4) << setiosflags(ios::left);
    
  out<<"SignalCode/F:"
     <<"Mass/F:"
     <<"Lumi/F:"
     <<"sLumi/F:"
     <<"Eff/F:"
     <<"sEff/F:"
     <<"DataEvts/F:"
     <<"BkgEvts/F:"
     <<"sBkgEvts/F:"
     <<"ObsLimit/F:"
     <<"ExpLimit/F:"
     <<"ExpLimitP1/F:"
     <<"ExpLimitM1/F:"
     <<"ExpLimitP2/F:"
     <<"ExpLimitM2/F"
     <<endl;

  TTree* tree = new TTree("tree", "Number of Events");
  tree->ReadFile(inName.c_str());
  tree->Draw("SignalCode:Mass:Lumi:DataEvts:BkgEvts:sBkgEvts:Eff:sEff", 
             "", "para goff");
  float n = tree->GetSelectedRows(); 
  for(int isample=0; isample<n; ++isample){
    const Double_t SignalCode = tree->GetVal(0)[isample];
    const Double_t  mass = tree->GetVal(1)[isample];
    const Double_t  lumi = tree->GetVal(2)[isample];
    const Double_t  DataEvts = tree->GetVal(3)[isample];
    const Double_t    BkgEvts = tree->GetVal(4)[isample];
    const Double_t   sBkgEvts = tree->GetVal(5)[isample];
    const Double_t       Eff = tree->GetVal(6)[isample];
    const Double_t      sEff = tree->GetVal(7)[isample];
    
    cout<<"Calculating limit for mass: "<<mass<<" and lumi: "<<lumi<<endl;
        
    float sLumi = sLumiFrac*lumi;
    
    Double_t obs_limit, exp_limit;
    Double_t exp_up, exp_down, exp_2up, exp_2down;
    
    if(useCLs){////CLs Limits
      //Does not work for bayesian, only works with cls    
      LimitResult limit = roostats_limit(lumi, sLumi, Eff, sEff, BkgEvts, sBkgEvts, DataEvts, false, 0, "cls", "", 12345);
      obs_limit = limit.GetObservedLimit();
      exp_limit = limit.GetExpectedLimit();
      exp_up    = limit.GetOneSigmaHighRange();
      exp_down  = limit.GetOneSigmaLowRange();
      exp_2up   = limit.GetTwoSigmaHighRange();
      exp_2down = limit.GetTwoSigmaLowRange();        
    }else{
      ////Bayesian Limits
      LimitResult limit  = roostats_clm (lumi, sLumi, Eff, sEff, BkgEvts, sBkgEvts);
      //obs_limit = limit.GetObservedLimit();
      obs_limit = roostats_cl95(lumi, sLumi, Eff, sEff, BkgEvts, sBkgEvts, DataEvts, false, 0, "bayesian", "");
      
      exp_limit = limit.GetExpectedLimit();
      exp_up    = limit.GetOneSigmaHighRange();
      exp_down  = limit.GetOneSigmaLowRange();
      exp_2up   = limit.GetTwoSigmaHighRange();
      exp_2down = limit.GetTwoSigmaLowRange();        
    }
  
    out<<setprecision(1)
       <<SignalCode<<"\t"
       <<setprecision(0)
       <<mass<<"\t"
       <<lumi<<"\t"
       <<sLumi<<"\t"
       <<setprecision(4)
       <<Eff<<"\t"
       <<sEff<<"\t"
       <<setprecision(0)
       <<DataEvts<<"\t"
       <<setprecision(4)
       <<BkgEvts<<"\t"
       <<sBkgEvts<<"\t"
//       <<setprecision(4)
       <<Value(obs_limit,-4)<<"\t"
       <<Value(exp_limit,-4)<<"\t"
       <<Value(exp_up,-4)<<"\t"
       <<Value(exp_down,-4)<<"\t"
       <<Value(exp_2up,-4)<<"\t"
       <<Value(exp_2down,-4)
       <<endl;
  }
  
  out.close(); 
  cout<<"Done\n";
  return;
}

