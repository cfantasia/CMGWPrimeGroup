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
    float N_aftercuts[Num_trkAlgos]; // # of expected events after all selection cuts
    float eff[Num_trkAlgos]; // efficiency of selection cuts wrt input sample
    float deff[Num_trkAlgos]; // efficiency uncertainty
    std::string samplename;
    std::string pathname;
    std::string description; // sample description
    TTree * tree; // pointer to ROOT file
    //
    InputFile()
    {
      x_sect = -1; Nprod_evt = -1; weight = 0; tree = 0;
      samplename = pathname = description = INVALID;
      for (int mual = 0; mual != Num_trkAlgos; ++mual)
	{
	  N_aftercuts[mual] = eff[mual] = deff[mual] = -1;
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
