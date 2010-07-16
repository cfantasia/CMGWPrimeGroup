#include "UserCode/CMGWPrimeGroup/root_macros/input_file.h"
#include "UserCode/CMGWPrimeGroup/interface/wprimeEvent.h"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include <TTree.h>
#include <TFile.h>

using std::string; using std::vector;
using std::cout; using std::endl; using std::cerr;

void parseLine(const string & new_line, wprime::InputFile * in_file);

string top_level_dir = "";

void loadInputFiles(vector<wprime::InputFile> & files, float lumiPb)
{
  ifstream topdir_file("UserCode/CMGWPrimeGroup/config/top_directory.txt");
  getline(topdir_file, top_level_dir);
  if(top_level_dir.empty())
    {
      cerr << " *** Failed to load top level directory! " << endl;
      return;
    }

  cout << "\n Reading root-tuples from directory " << top_level_dir << endl;
  cout << " Applying weights to get distributions for " << lumiPb << " pb^-1"
       << endl << endl;

  ifstream in("UserCode/CMGWPrimeGroup/config/samples_cross_sections.txt");
  string new_line; wprime::InputFile * new_file = 0;
  while(getline(in, new_line))
    {
      // if DONE, we are done!
      if(new_line == "DONE")break;

      if(new_line.find("samplename = ") != string::npos)
	// new file found! create structure to put in info
	new_file = new wprime::InputFile();
      
      if(!new_line.empty())
	parseLine(new_line, new_file);
      else
	{
	  new_file->weight = lumiPb*(new_file->x_sect)
	    /(new_file->Nprod_evt);
	  cout << "; Weight to be applied on " << new_file->samplename
	       << " sample = " << new_file->weight << endl << endl;
	  // all info should now be logged in; check!
	  new_file->checkFile();
	  // if we made it here, everything looks good: add to vector of MC samples
	  files.push_back(*new_file);
	  // release memory
	  delete new_file;
	}
    }
}

void parseLine(const string & new_line, wprime::InputFile * in_file)
{
  unsigned int i = 0;
  i = new_line.find("samplename = ");
  if(i != string::npos)
    {
      in_file->samplename = new_line.substr(13, new_line.length() - 13);
      return;
    }

  i = new_line.find("description = ");
  if(i != string::npos)
    {
      in_file->description = new_line.substr(14, new_line.length() - 14);
      return;
    }

  i = new_line.find("pathname = ");
  if(i != string::npos)
    {
      string pathname = top_level_dir + new_line.substr(11, new_line.length() - 11);
      in_file->pathname = pathname;

      TFile * file = new TFile(pathname.c_str());
      if(!(file->IsOpen()))
	{
	  cerr <<" *** Missing file: "<< pathname << " !!! "<<endl; 
	  return;
	}
      in_file->tree = (TTree *) file->Get("StdMu/wprime");
      if(! in_file->tree->GetBranch("wp"))
	{
	  cerr << " *** Can't find wp branch in file: " << pathname<< endl;
	  in_file->tree = 0;
	  return;
	}
      
      cout << " Opened file: " << in_file->samplename;
      cout << " (" << in_file->description << ") " << endl;
      cout << " # of entries = " << in_file->tree->GetEntries();

      return;
    }

  i = new_line.find("x-section = ");
  if(i != string::npos)
    {
      in_file->x_sect = atof(new_line.substr(12, new_line.length() - 12).c_str());
      return;
    }

  i = new_line.find("Nprod_evt = ");
  if(i != string::npos)
    {
      in_file->Nprod_evt = atoi(new_line.substr(12, new_line.length() - 12).c_str());
      return;
    }
  
}
