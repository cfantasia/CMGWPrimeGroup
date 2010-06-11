#include <TROOT.h>
#include <TFile.h>

#include <string>
#include <iostream>
#include <fstream>

#include "UserCode/CMGWPrimeGroup/root_macros/input_file.h"

using std::cout; using std::endl;


extern int loadCrossSections(vector<wprime::InputFile> & qcd_files, 
			     vector<wprime::InputFile> & z_files,
			     vector<wprime::InputFile> & w_files,
			     vector<wprime::InputFile> & top_files,
			     vector<wprime::InputFile> & wprime10_files,
			     vector<wprime::InputFile> & wprime15_files,
			     vector<wprime::InputFile> & wprime20_files,
			     unsigned int detector_conditions); 

extern void GetDistributionGeneric(const vector<wprime::InputFile>& files, 
				   TFile *fout, string dir, ofstream & out, 
				   const int option = 1, 
				   const bool highestPtMuonOnly = true);

extern void loadInputFiles(string file_desc, vector<wprime::InputFile> & files, 
			   float lumiPb, unsigned int detector_conditions);

void run()
{
  // flag for detector conditions 
  // options:
  // 1 -> 21x RECO, ideal conditions
  // 2 -> 2212 RECO, 50 pb-1 alignment for Wprime and W (all other bgd samples same as option 1)
  unsigned int detector_conditions = 2;
  
  float lumiPb = 100; // in pb^-1

  TFile *fout = new TFile("Wprime_analysis.root","recreate");

  vector<wprime::InputFile> qcd_files; vector<wprime::InputFile> z_files;
  vector<wprime::InputFile> w_files; vector<wprime::InputFile> top_files;
  vector<wprime::InputFile> wprime10_files;
  vector<wprime::InputFile> wprime15_files;
  vector<wprime::InputFile> wprime20_files;

  loadCrossSections(qcd_files, z_files, w_files,top_files, wprime10_files, wprime15_files, wprime20_files, detector_conditions);

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
  
  
  string dir = "QCD";
  loadInputFiles(dir, qcd_files, lumiPb, detector_conditions);
  GetDistributionGeneric(qcd_files, fout, dir, out, option, highestPtMuonOnly);
  
  dir = "Z";
  loadInputFiles(dir, z_files, lumiPb, detector_conditions);
  GetDistributionGeneric(z_files, fout, dir, out, option, highestPtMuonOnly);  

  dir = "W";
  loadInputFiles(dir, w_files, lumiPb, detector_conditions);
  GetDistributionGeneric(w_files, fout, dir, out, option, highestPtMuonOnly);    

  dir = "Top";
  loadInputFiles(dir, top_files, lumiPb, detector_conditions);
  GetDistributionGeneric(top_files, fout, dir, out, option, highestPtMuonOnly);    
  
  dir = "wprime10";
  loadInputFiles(dir, wprime10_files, lumiPb, detector_conditions);
  GetDistributionGeneric(wprime10_files, fout, dir, out, option, highestPtMuonOnly);    
  
  dir = "wprime15";
  loadInputFiles(dir, wprime15_files, lumiPb, detector_conditions);
  GetDistributionGeneric(wprime15_files, fout, dir, out, option, highestPtMuonOnly);    
  
  dir = "wprime20";
  loadInputFiles(dir, wprime20_files, lumiPb, detector_conditions);
  GetDistributionGeneric(wprime20_files, fout, dir, out, option, highestPtMuonOnly);    
  
  out.close(); 
  fout->Close();
}
