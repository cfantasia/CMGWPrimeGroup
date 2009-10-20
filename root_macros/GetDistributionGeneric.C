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

// option 1
extern void GetMuonPtDistribution_JetIso(const vector<wprime::InputFile>& files, 
					 TFile *fout, string dir);
// option 2
extern void GetMuonPtDistribution(const vector<wprime::InputFile>& files, 
				  TFile *fout, string dir, ofstream & out);

// option 3
extern void GetChargePtDistribution(const vector<wprime::InputFile>& files, 
				    TFile *fout, string dir);

const int option = 2;

void GetDistributionGeneric(const vector<wprime::InputFile>& files, 
			    TFile *fout, string dir, ofstream & out)
{
  if(option == 1)
    GetMuonPtDistribution_JetIso(files, fout, dir);
  else if (option == 2)
    GetMuonPtDistribution(files, fout, dir, out);
  else if (option == 3)
    GetChargePtDistribution(files, fout, dir);

  return;
}
