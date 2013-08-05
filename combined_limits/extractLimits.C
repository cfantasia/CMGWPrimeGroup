//Author: Cory Fantasia 2012
//Purpose: Extract limits from higgs output for later plotting
//Usage: root -b -l -q 'extractLimits.C+(string inName, string modeName, string algoName, float xSecScale=1)'

#include "../root_macros/common.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include <fstream>

double getMedian(TTree* tree, const string & cuts);
double getMedian(TTree* tree, const string & cuts, double & lo95,double & lo68,double & hi68,double & hi95);

void
extractLimits(string inName, string modeName, string algoName, float xSecScale=1.){
  TFile *fLumi = TFile::Open(inName.c_str(), "read"); assert(fLumi);
  float lumi = GetLumiUsed(fLumi);


  string outfile(Form("nLimit_%s_%s.txt", modeName.c_str(), algoName.c_str()));
  ofstream out(outfile.c_str());
  if(!out) { 
    cout << "Cannot open file " << outfile << endl; 
    abort();
  } 
  
  out  << setiosflags(ios::fixed) << setprecision(4) << setiosflags(ios::left);
//out<<"SignalCode/F:"
  out<<"Mass/F:"
     <<"XSec/F:"
     <<"NToys/F:"
     <<"Lumi/F:"
//     <<"sLumi/F:"
//     <<"Eff/F:"
//     <<"sEff/F:"
//     <<"DataEvts/F:"
//     <<"BkgEvts/F:"
//     <<"sBkgEvts/F:"
     <<"ObsLimit/F:"
     <<"ExpLimit/F:"
     <<"ExpLimitP1/F:"
     <<"ExpLimitM1/F:"
     <<"ExpLimitP2/F:"
     <<"ExpLimitM2/F"
     <<endl;

  TGraph* gxsec = new TGraph("../Limits/xSec_WZ.dat", "%lg %lg");
  
  const int mass_inc = 1;
  for(int mass=100; mass<=2000; mass+=mass_inc){
    //string inputFile = Form("higgsCombineWprimeWZ.ProfileLikelihood.mH%i.root", mass);
    //string inputFile = Form("higgsCombineWprimeWZ.MarkovChainMC.mH%i.root", mass);
    string inputFile = Form("higgsCombine%s.%s.mH%i.root", modeName.c_str(), algoName.c_str(), mass);

    if(1 == gSystem->AccessPathName(inputFile.c_str())) continue;//ret val 1 means not found
    TFile *fIn = TFile::Open(inputFile.c_str(), "read"); if(!fIn) continue;
    
    //float mass = atof(inputFile.substr(inputFile.find(".mH")+3, inputFile.find(".root")).c_str());
    //cout<<"read mass of "<<mass<<" from "<<inputFile<<endl;

    TTree* tree = (TTree*) fIn->Get("limit"); if(!tree) continue;
    int nLimits = tree->Draw("limit", "quantileExpected==-1");
    double obs_limit,medianLimit,lo95,lo68,hi68,hi95;
    obs_limit = medianLimit = lo95 = lo68 = hi68 = hi95 = 0;

    if(algoName == "MarkovChainMC" || algoName == "Asymptotic"){
      obs_limit = getMedian(tree, "iToy==0");
      medianLimit = getMedian(tree, "iToy>0", lo95,lo68,hi68,hi95);
    }else{//Not MarkovChainMC
       obs_limit   = getMedian(tree, "quantileExpected==-1"); 
       medianLimit = getMedian(tree, "quantileExpected==0.5"); 
      
       lo95 = getMedian(tree, "0.0 < quantileExpected && quantileExpected < 0.1"); 
       lo68 = getMedian(tree, "0.1 < quantileExpected && quantileExpected < 0.2"); 
       hi68 = getMedian(tree, "0.8 < quantileExpected && quantileExpected < 0.9"); 
       hi95 = getMedian(tree, "0.9 < quantileExpected && quantileExpected < 1.0");      
    }
    double xsec = gxsec->Eval(mass);
    xsec *= xSecScale;
    out<<setprecision(0)
       <<mass<<"\t"
       <<setprecision(9)
       <<xsec<<"\t"
       <<setprecision(0)
       <<nLimits<<"\t"
       <<lumi<<"\t";
    out<<setprecision(8)
       <<xsec * obs_limit<<"\t"
       <<xsec * medianLimit<<"\t"
       <<xsec * hi68<<"\t"
       <<xsec * lo68<<"\t"
       <<xsec * hi95<<"\t"
       <<xsec * lo95
       <<endl;
    
    delete fIn;
  }//mass loop
 out.close();
}

double
getMedian(TTree* tree, const string & cuts){
  double lo95(0), lo68(0), hi68(0), hi95(0);
  return getMedian(tree,cuts,lo95,lo68,hi68,hi95);
}

double
getMedian(TTree* tree, const string & cuts, double & lo95,double & lo68,double & hi68,double & hi95){
  vector<float> limitHistory;
  
  int nLimits = tree->Draw("limit", cuts.c_str());
  if (nLimits <= 0){
    cout<<" Failed to find any lines in tree\n";
    return 0.;
  }
  Double_t*  limit = tree->GetVal(0);
  for(int ievt=0; ievt<nLimits; ++ievt){
    limitHistory.push_back(limit[ievt]);
  }
  
  sort(limitHistory.begin(), limitHistory.end());
  
  double medianLimit = (nLimits % 2 == 0 ? 0.5*(limitHistory[nLimits/2-1]+limitHistory[nLimits/2]) : limitHistory[nLimits/2]);
  hi95 = limitHistory[min<int>(nLimits-1,  ceil(0.975 * nLimits))];
  hi68 = limitHistory[min<int>(nLimits-1,  ceil(0.84  * nLimits))];
  lo68 = limitHistory[min<int>(nLimits-1, floor(0.16  * nLimits))];
  lo95 = limitHistory[min<int>(nLimits-1, floor(0.025 * nLimits))];
    //cout << "   68% expected band : " << lo68 << " < r < " << hi68 << endl;
    //cout << "   95% expected band : " << lo95 << " < r < " << hi95 << endl;
  return medianLimit;
}

/*
    vector<float> limitHistory;
  
    TTree* tree = (TTree*) fIn->Get("limit"); if(!tree) continue;
    tree->Draw("limit", "quantileExpected==0.5");
    int nLimits = tree->GetSelectedRows(); 
    Double_t*  limit = tree->GetVal(0);
    for(int ievt=0; ievt<nLimits; ++ievt){
      limitHistory.push_back(limit[ievt]);
    }

    sort(limitHistory.begin(), limitHistory.end());
    if (nLimits <= 0) continue;
    double obs_limit = 0.;//Fix..........
    
    double medianLimit = (nLimits % 2 == 0 ? 0.5*(limitHistory[nLimits/2-1]+limitHistory[nLimits/2]) : limitHistory[nLimits/2]);
    //cout << "median expected limit: r < " << medianLimit << " @ 95%CL (" <<nLimits << " toyMC)" << endl;
    double hi68 = limitHistory[min<int>(nLimits-1,  ceil(0.84  * nLimits))];
    double lo68 = limitHistory[min<int>(nLimits-1, floor(0.16  * nLimits))];
    double hi95 = limitHistory[min<int>(nLimits-1,  ceil(0.975 * nLimits))];
    double lo95 = limitHistory[min<int>(nLimits-1, floor(0.025 * nLimits))];
    //cout << "   68% expected band : " << lo68 << " < r < " << hi68 << endl;
    //cout << "   95% expected band : " << lo95 << " < r < " << hi95 << endl;
*/
