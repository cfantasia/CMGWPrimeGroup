#ifndef _util_h
#define _util_h

#include <vector>
#include <string>
#include <map>

namespace wprime{
  static std::string INVALID = "INVALID";

  struct InputFile
  {
    float x_sect; // cross-section in pb
    int Nprod_evt;// unweighted # of events produced (should correspond to x_sect!)
    int Nact_evt; // unweighted # of events surviving pre-selection/skimming and
    // actually contained in input file
    float weight; // cross-section * integrated luminosity / (# of events produced)
    // Nact_evt * weight = Nexp_evt for given integrated luminosity
    std::string samplename;
    std::vector<std::string> pathnames; // directory + filenames
    std::string description; // sample description
    //
    InputFile()
    {
      x_sect = -1; Nprod_evt = Nact_evt = -1; weight = 0;
      samplename = description = INVALID;
    }
    void checkFile()
    {
      assert(x_sect > 0); assert(Nprod_evt > 0); //assert(Nact_evt > 0);
      assert(weight > 0);
      assert(pathnames.size()); assert(description != INVALID); 
      assert(samplename != INVALID);
    }

  };



  struct FilterEff{
    // # of (unweighted!) events surviving after each selection cut
    int Nsurv_evt_cut;
    // # of (weighted!) events surviving after each selection cut
    float Nsurv_evt_cut_w;
    // efficiency for each selection cut
    float eff;
    // efficiency uncertainty
    float deff;

    FilterEff(){Nsurv_evt_cut = 0; Nsurv_evt_cut_w = eff = deff = 0.0;}
  };

  // key: samplename, value: vector<FilterEff> (ie. statistics for selection steps)
  typedef std::map<std::string, std::vector<FilterEff> > SampleStat;

  const float MUON_MASS = 0.105658366;      // GeV
  const float ELECTRON_MASS = 0.000511;     // GeV
}



#endif // _util_h
