#ifndef _wprime_input_file_h
#define _wprime_input_file_h

#include <vector>
#include <string>

#include "UserCode/CMGWPrimeGroup/root_macros/wprime_histo_constants.h"

class TTree;

namespace wprime{
  static std::string INVALID = "INVALID";

  struct InputFile
  {
    float x_sect; // cross-section in pb
    int Nprod_evt; // # of events produced (before any analysis cuts or other filtering)
    float weight; // cross-section * integrated luminosity / (# of events produced)
    // # of expected events after each selection cut, per reconstructor 
    // and for pt/Mt analysis
    float Nexp_evt_cut[Num_selection_cuts][Num_trkAlgos][Num_flavors]; 
    // efficiency for each selection cut, per reconstructor 
    // and for pt/Mt analysis wrt input sample
    float eff[Num_selection_cuts][Num_trkAlgos][Num_flavors]; 
    // efficiency uncertainty
    float deff[Num_selection_cuts][Num_trkAlgos][Num_flavors]; 
    std::string samplename;
    std::string pathname;
    std::string description; // sample description
    TTree * tree; // pointer to ROOT file
    //
    InputFile()
    {
      x_sect = -1; Nprod_evt = -1; weight = 0; tree = 0;
      samplename = pathname = description = INVALID;
      for(int i = 0; i != Num_selection_cuts; ++i)
	for(int f = 0; f != Num_flavors; ++f)
	  //	  for (int mual = MuAlgo_MIN; mual <= MuAlgo_MAX; ++mual)
	  for (int mual = 0; mual != Num_trkAlgos; ++mual)
	    {
	      Nexp_evt_cut[i][mual][f] = eff[i][mual][f] = 
		deff[i][mual][f] = 0;
	    }
    }
    void checkFile()
    {
      assert(x_sect > 0); assert(Nprod_evt > 0); assert(weight > 0);
      assert(pathname != INVALID); assert(description != INVALID); 
      assert(samplename != INVALID);
    }

  };
}



#endif // _wprime_input_file_h
