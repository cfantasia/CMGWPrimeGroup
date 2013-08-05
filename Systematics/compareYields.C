//Usage: root -b -l -q 'printSysTable.C+("../Systematics/SysBkgPDF.dat")'

#include <fstream>
#include "../root_macros/common.h"

const int nch = 4;
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
    vector<string> bkgSamples = BkgSamples();
    fBkg<<mass;
    //cout<<mass<<endl;
    for(int ch=0; ch<nch; ++ch){
      string cuts = Form("weight*(%s)*(EvtType == %i)", AnalysisCuts(mass).c_str(), ch);
      float corrCoef = CalcCorrCoef(getTree(fOrig, bkgSamples, "tEvts_MET"), getTree(fMod, bkgSamples, "tEvts_MET"), cuts);

      Value origYield = GetNEvtsAndError(fOrig, bkgSamples, "tEvts_MET", cuts);
      Value  modYield = GetNEvtsAndError(fMod , bkgSamples, "tEvts_MET", cuts);
      Value sys = ShiftErr(origYield, modYield, corrCoef);
      fBkg<<"\t"<<Form("%7.4f %7.4f", sys.val, sys.err);
      origYield.setPrintErr(1); modYield.setPrintErr(1); sys.setPrintErr(1);
      //cout<<" : "<<origYield<<" : "<<modYield<<" : "<<sys<<endl;
    }
    fBkg<<endl;
  }
  for(int mass=200; mass<=2000; mass+=100){
    string sample = Form("WprimeToWZTo3LNu_M-%i", mass);
    fSig<<mass;
    //cout<<mass<<endl;
    for(int ch=0; ch<nch; ++ch){
      string cuts = Form("weight*(%s)*(EvtType == %i)", AnalysisCuts(mass).c_str(), ch);
      float corrCoef = CalcCorrCoef(getTree(fOrig, sample, "tEvts_MET"), getTree(fMod, sample, "tEvts_MET"), cuts);

      Value origYield = GetNEvtsAndError(fOrig, sample, "tEvts_MET", cuts);
      Value  modYield = GetNEvtsAndError(fMod , sample, "tEvts_MET", cuts);
      Value sys = ShiftErr(origYield, modYield, corrCoef);
      fSig<<"\t"<<Form("%7.4f %7.4f", sys.val, sys.err);
      origYield.setPrintErr(1); modYield.setPrintErr(1); sys.setPrintErr(1);
      //cout<<" : "<<origYield<<" : "<<modYield<<" : "<<sys<<endl;
    }
    fSig<<endl;
  }
}
