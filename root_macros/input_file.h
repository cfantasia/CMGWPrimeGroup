#ifndef _wprime_input_file_h
#define _wprime_input_file_h

#include <vector>
#include <string>

class TTree;

namespace wprime{

  struct InputFile
  {
    float x_sect; // cross-section in pb
    int Nprod_evt; // # of events produced (before any analysis cuts or other filtering)
    float weight; // cross-section * integrated luminosity / (# of events produced)
    std::string pathname;
    std::string description; // sample description
    TTree * tree; // pointer to ROOT file
    //
    InputFile(float cross_sect, int nprod_evt)
    {
      x_sect = cross_sect; Nprod_evt = nprod_evt; weight = 0; tree = 0;
    }
  };
}



#endif // _wprime_input_file_h
