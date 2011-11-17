/* To Use: 
   root -l CalcLimit.C
*/
#include "consts.h"
#include "TROOT.h"
#include "TSystem.h"
//#include "../../../StatisticalTools/RooStatsRoutines/root/roostats_cl95.C"

void
CalcLimit(){
  gErrorIgnoreLevel = kWarning;
//  gSystem->SetIncludePath( "-I$ROOFITSYS/include" );
  gSystem->SetIncludePath( "-I/afs/hep.wisc.edu/cern/.root/root_v5.30.00.Linux-slc5_amd64-gcc4.3/include/RooStats" );
  gROOT->ProcessLine(".L ../../../StatisticalTools/RooStatsRoutines/root/roostats_cl95.C+");
  
  string outfile("nLimit.txt");
  ofstream out(outfile.c_str());
  if(!out) { 
    cout << "Cannot open file " << outfile << endl; 
    abort();
  } 
  
  out  << setiosflags(ios::fixed) << setprecision(4) << setiosflags(ios::left);
    
  out<<"SignalCode/I:"
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
  tree->ReadFile("nEvents.txt");
  tree->Draw("SignalCode:Mass:Lumi:DataEvts:BkgEvts:sBkgEvts:Eff:sEff", 
             "", "para goff");
  float n = tree->GetSelectedRows(); 
  for(int isample=0; isample<n; ++isample){
    const int SignalCode = tree->GetVal(0)[isample];
    const Double_t  mass = tree->GetVal(1)[isample];
    const Double_t  lumi = tree->GetVal(2)[isample];
    const Double_t  DataEvts = tree->GetVal(3)[isample];
    const Double_t    BkgEvts = tree->GetVal(4)[isample];
    const Double_t   sBkgEvts = tree->GetVal(5)[isample];
    const Double_t       Eff = tree->GetVal(6)[isample];
    const Double_t      sEff = tree->GetVal(7)[isample];
    
    cout<<"Calculating limit for mass: "<<mass<<" and lumi: "<<lumi<<endl;
        
    float sLumi = sLumiFrac*lumi;
    
    /*
    ////CLs Limits
    //Does not work for bayesian, only works with cls    
    LimitResult limit = roostats_limit(lumi, sLumi, Eff, sEff, BkgEvts, sBkgEvts, DataEvts, false, 0, "cls", "");
    Double_t obs_limit = limit.GetObservedLimit();
    Double_t exp_limit = limit.GetExpectedLimit();
    Double_t exp_up    = limit.GetOneSigmaHighRange();
    Double_t exp_down  = limit.GetOneSigmaLowRange();
    Double_t exp_2up   = limit.GetTwoSigmaHighRange();
    Double_t exp_2down = limit.GetTwoSigmaLowRange();        
    */

    ////Bayesian Limits
    LimitResult limit  = roostats_clm (lumi, sLumi, Eff, sEff, BkgEvts, sBkgEvts);
    //Double_t obs_limit = limit.GetObservedLimit();
    Double_t obs_limit = roostats_cl95(lumi, sLumi, Eff, sEff, BkgEvts, sBkgEvts, DataEvts, false, 0, "bayesian", "");
    Double_t exp_limit = limit.GetExpectedLimit();
    Double_t exp_up    = limit.GetOneSigmaHighRange();
    Double_t exp_down  = limit.GetOneSigmaLowRange();
    Double_t exp_2up   = limit.GetTwoSigmaHighRange();
    Double_t exp_2down = limit.GetTwoSigmaLowRange();        
    
  
    out<<setprecision(0)
       <<SignalCode<<"\t"
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
       <<obs_limit<<"\t"
       <<exp_limit<<"\t"
       <<exp_up<<"\t"
       <<exp_down<<"\t"
       <<exp_2up<<"\t"
       <<exp_2down
       <<endl;
  }
  
  out.close(); 
  cout<<"Done\n";
  return;
}

