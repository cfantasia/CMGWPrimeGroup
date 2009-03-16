#include "UserCode/CMGWPrimeGroup/root_macros/input_file.h"
#include "UserCode/CMGWPrimeGroup/interface/wprimeEvent.h"

#include <string>
#include <vector>
#include <iostream>

#include <TTree.h>
#include <TFile.h>

using std::string; using std::vector;
using std::cout; using std::endl; using std::cerr;

// string top_level_dir = "UserCode/CMGWPrimeGroup/"; // Alessio
string top_level_dir = "/home/cleonido/wprime/v52/"; // Christos

int checkFiles(unsigned Nfiles, vector<wprime::InputFile> & files)
{
  if(Nfiles != files.size())
    {
      cerr << " *** Expected " << Nfiles << " files, got "
	   << files.size() << " instead..." << endl;
      return -1;
    }  
  return 0;
}

void loadInputFiles(string file_desc, vector<wprime::InputFile> & files, 
		   float lumiPb)
{
  int Nfiles = -1;
  cout << "\n Processing " << file_desc << " files " << endl << endl;

  // ==========================================================================
  if(file_desc == "QCD")
    {
      const int NfilesQCD = 9;
      Nfiles = NfilesQCD;
      if(checkFiles(Nfiles, files))
	return;

      string low[NfilesQCD]= 
	{"100", "150", "200", "300", "400", "600", "800", "1200", "1600"};
      string high[NfilesQCD]=
	{"150", "200", "300", "400", "600", "800","1200", "1600", "up"};
  
      for(int i = 0; i != Nfiles; ++i)
	{
	  files[i].pathname = top_level_dir + string("QCD_")+low[i]+string("_")+high[i]+string("_212_Ideal_Minv_ptGlobMu.root");
	}
      
    } // QCD

  // ==========================================================================
  else if(file_desc == "Z")
    {
      const int NfilesZ = 11;
      Nfiles = NfilesZ;        
      if(checkFiles(Nfiles, files))
	return;
      
      string lowZ[NfilesZ] = {"30", "110","200", "300", "400", "500", "600", "700", "800", "900", "1000"};
  
      for(int i = 0; i != Nfiles; ++i)
	{
	  files[i].pathname = top_level_dir + string("Z_212_Ideal_Minv_ptGlobalMu_pt")+lowZ[i]+string(".root");
	}
      
    } // Z

  // ==========================================================================
  else if(file_desc == "W")
    {
      const int NfilesW = 22;
      Nfiles = NfilesW;
      if(checkFiles(Nfiles, files))
	return;

      string lowW[NfilesW]={"0", "200", "250", "300", "350", "400", "450", "500", "600", "700", "800", "900", "1000", "1100", "1200", "1300", "1400", "1500", "1600", "1700", "1800", "1900"};
      string highW[NfilesW]={"200", "250", "300", "350", "400", "450", "500", "600", "700", "800", "900", "1000", "1100", "1200", "1300", "1400", "1500", "1600", "1700", "1800", "1900","2000"};
  
      for(int i = 0; i != Nfiles; ++i)
	files[i].pathname = top_level_dir + string("Wmunu_219_Ideal_Minv_")+lowW[i]+string("_")+highW[i]+string(".root");
      
    } // W
  // ==========================================================================
  else if(file_desc == "Top")
    {
      const int NfilesTop = 1;
      Nfiles = NfilesTop;
      if(checkFiles(Nfiles, files))
	return;

      files[0].pathname = top_level_dir + string("TTbar_Minv_PtMuGlobal.root");
      
    } // Top
  // ==========================================================================
  else if(file_desc == "wprime10"  || file_desc == "wprime15" ||
	  file_desc == "wprime20")
    {
      const int NfilesWprime = 1;
      Nfiles = NfilesWprime;
      if(checkFiles(Nfiles, files))
	return;

      if(file_desc == "wprime10")
	 files[0].pathname = top_level_dir + string("wprime_1TeV.root");
      else if(file_desc == "wprime15")
	files[0].pathname = top_level_dir + string("wprime_1.5TeV.root");
      else if(file_desc == "wprime20")
	files[0].pathname = top_level_dir + string("wprime_2TeV.root");
      
    } // Top
  // ==========================================================================
  else
    {
      cerr << " Unknown file description " << file_desc << endl;
      return;
    }

  // ==========================================================================

  for(int i = 0; i != Nfiles; ++i)
    { // loop over input files

      string pathname = files[i].pathname;
      files[i].weight = lumiPb*(files[i].x_sect)/(files[i].Nprod_evt);

      TFile * file = new TFile(pathname.c_str());
      if(!(file->IsOpen()))
	{
	  cerr <<" *** Missing file: "<< pathname << " !!! "<<endl; 
	  continue;
	}
      files[i].tree = (TTree *) file->Get("StdMu/wprime");
      if(! files[i].tree->GetBranch("wp"))
	{
	  cerr << " *** Can't find wp branch in file: " << pathname<< endl;
	  files[i].tree = 0;
	  continue;
	}
      TTree * jb = (TTree *) file->Get("StdMu/jobinfo");
      if(! jb->GetBranch("job"))
	{
	  cerr << " *** Can't find job branch in file: " << pathname<< endl;
	  continue;
	}
      wprime::JobInfo * job = new wprime::JobInfo();
      jb->SetBranchAddress("job", &job); jb->GetEntry(0);

      cout << " Opened file: " << job->sample;
      //      cout << " RECO version = " << job->RECOversion
      //	   << " HLT version = " << job->HLTversion << endl;
      cout << ",  # of entries = " << files[i].tree->GetEntries() 
	   << ", weight = " << files[i].weight << endl;

      files[i].description = job->sample;

      if(files[i].Nprod_evt != job->Nprod_evt)
	cout << " *** Oooops! ROOT file has Nprod_evt = " << job->Nprod_evt
	     << ", expected to find " << files[i].Nprod_evt << endl;
      

    } // loop over input files
      
  return;
}
