//Author: Cory Fantasia 2012
//Purpose: Print unweighted MC yields and errors
//Usage: root -b -l -q 'printSysTable.C+("../Systematics/SysBkgPDF.dat")'

#include <fstream>
#include "../root_macros/common.h"

const int nch = 4;
void
nEvtsUnweighted(){
  TFile *fOrig = TFile::Open("../../../WprimeWZ.root", "read"); assert(fOrig);

  for(int mass=200; mass<=2000; mass+=100){
    vector<string> bkgSamples = BkgSamples();
    cout<<mass;
    for(int ch=0; ch<nch; ++ch){
      string cuts = Form("(%s)*(EvtType == %i)", AnalysisCuts(mass).c_str(), ch);

      Value origYield = GetNEvtsAndError(fOrig, bkgSamples, "tEvts_MET", cuts);
      cout<<"\t"<<Form("%7.4f %7.4f", origYield.val, origYield.err);
    }
    cout<<endl;
  }

  for(int mass=200; mass<=2000; mass+=100){
    string sample = Form("WprimeToWZTo3LNu_M-%i", mass);
    cout<<mass;
    //cout<<mass<<endl;
    for(int ch=0; ch<nch; ++ch){
      string cuts = Form("(%s)*(EvtType == %i)", AnalysisCuts(mass).c_str(), ch);

      Value origYield = GetNEvtsAndError(fOrig, sample, "tEvts_MET", cuts);
      cout<<"\t"<<Form("%7.4f %7.4f", origYield.val, origYield.err);
    }
    cout<<endl;
  }
}
