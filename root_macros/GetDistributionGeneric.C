// ROOT stuff
#include <TROOT.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

// wprime stuff
#include "UserCode/CMGWPrimeGroup/root_macros/input_file.h"

// std stuff
#include <iostream>

extern void GetMuonPtDistribution(const vector<wprime::InputFile>& files, 
				  TFile *fout, string dir);

extern void GetMuonPtDistribution_JetIso(const vector<wprime::InputFile>& files, 
					 TFile *fout, string dir);

void GetDistributionGeneric(const vector<wprime::InputFile>& files, 
			    TFile *fout, string dir)
{
  GetMuonPtDistribution(files, fout, dir);
  //  GetMuonPtDistribution_JetIso(files, fout, dir);
  return;
}
