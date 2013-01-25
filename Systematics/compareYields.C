//Usage: root -b -l -q 'printSysTable.C+("../Systematics/SysBkgPDF.dat")'

#include <fstream>
#include "../root_macros/common.h"

void
compareYields(const string origName, const string modName, const string outName){
  TFile *fOrig = TFile::Open(origName.c_str(), "read"); assert(fOrig);
  TFile *fMod  = TFile::Open( modName.c_str(), "read"); assert(fMod);

  string outSigName = outName;
  outSigName.replace(outSigName.find("Sys"), 3, "SysSig"); 

  string outBkgName = outName;
  outBkgName.replace(outBkgName.find("Sys"), 3, "SysBkg"); 
  
  ofstream fSig(outSigName.c_str());
  if(!fSig) { 
    cout << "Cannot open file " << outSigName << endl; 
    abort();
  } 
  ofstream fBkg(outBkgName.c_str());
  if(!fBkg) { 
    cout << "Cannot open file " << outBkgName << endl; 
    abort();
  } 

  for(int mass=200; mass<=2000; mass+=100){
    string cuts = Form("weight*(%s)", AnalysisCuts(mass).c_str());
    string sample = "WZJetsTo3LNu";
    float origYield = GetNEvts(fOrig, sample, "tEvts_MET", cuts);
    float modYield  = GetNEvts(fMod , sample, "tEvts_MET", cuts);
    fSig<<mass<<"\t"<<fabs(modYield - origYield) / origYield<<endl;
  }

  for(int mass=200; mass<=2000; mass+=100){
    string cuts = Form("weight*(%s)", AnalysisCuts(mass).c_str());
    string sample = Form("WprimeToWZTo3LNu_M-%i", mass);
    float origYield = GetNEvts(fOrig, sample, "tEvts_MET", cuts);
    float modYield  = GetNEvts(fMod , sample, "tEvts_MET", cuts);
    fBkg<<mass<<"\t"<<fabs(modYield - origYield) / origYield<<endl;
  }
}
