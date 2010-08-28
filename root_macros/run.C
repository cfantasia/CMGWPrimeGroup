#include <TROOT.h>
#include <TFile.h>

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "UserCode/CMGWPrimeGroup/root_macros/input_file.h"

using std::cout; using std::endl;
using std::vector; using std::string;

extern void GetDistributionGeneric(const wprime::InputFile & file, 
				   TFile *fout, ofstream & out, 
				   const int option = 1, 
				   const bool highestPtMuonOnly = true);

extern void loadInputFiles(vector<wprime::InputFile> & files, float lumiPb);

void run()
{
  //  float lumiPb = 100; // in pb^-1
  float lumiPb = 0.840; // in pb^-1

  TFile *fout = new TFile("Wprime_analysis.root","recreate");

  vector<wprime::InputFile> all_files; 

  string outfile("event_counts.txt");
  ofstream out(outfile.c_str());
  if(!out) { 
    cout << "Cannot open file " << outfile << endl; 
    abort();
  } 

  // option = 1 : GetMuonPtDistribution
  // option = 2 : GetChargePtDistribution
  const int option = 1;
  
  // switch on or off looping over muons;
  // if true, will only examine hardest muon in event (based on tracker-pt)
  bool highestPtMuonOnly = false; 

  loadInputFiles(all_files, lumiPb);
  
  vector<wprime::InputFile>::const_iterator it;
  for(it = all_files.begin(); it != all_files.end(); ++it)
    GetDistributionGeneric(*it, fout, out, option, highestPtMuonOnly);

  out.close(); 
  fout->Close();
}
