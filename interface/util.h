#ifndef _util_h
#define _util_h

#include <vector>
#include <string>
#include <map>

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

////
#include <stdarg.h>
#include <boost/algorithm/string.hpp>

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"

#include "DataFormats/Common/interface/MergeableCounter.h"
////

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
    std::string subdir;
    std::vector<std::string> pathnames; // directory + subdir + filenames
    std::string description; // sample description
    //
    InputFile()
    {
      x_sect = -1; Nprod_evt = Nact_evt = -1; weight = 0;
      samplename = description = INVALID;
      subdir = "";
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

//// Global constants ////////////////////////////////////////////////////////

const float ZMASS = 91.188;
const float WMASS = 80.398;

const int PDGMUON = 13;
const int PDGELEC = 11;
const int PDGW = 24;
const int PDGZ = 23;
const int PDGWPRIME = 34;

const float PI    = 2.0 * TMath::ACos(0.);
const float TWOPI = 2.0 * PI;
const float NOCUT = 9e9;

//////// Comparators (for sorting vectors) ///////////////////////////////////

struct closestToZMass {
  bool operator() (const reco::Candidate & a, const reco::Candidate & b) {
    return fabs(ZMASS - a.mass()) < fabs(ZMASS - b.mass());
  }
};
struct highestPt {
  bool operator() (const reco::Candidate * a, const reco::Candidate * b) {
    return a->pt() > b->pt();
  }
};
struct highestPtLepton {
  bool operator() (const reco::CompositeCandidate a, 
                   reco::CompositeCandidate b) {
    return a.daughter(0)->pt() > b.daughter(0)->pt();
  }
};


inline bool areIdentical(const reco::Candidate & p1, 
                         const reco::Candidate & p2)
{
  float tolerance = 0.0001;
  if (p1.pdgId() == p2.pdgId() &&
      fabs(p1.eta() - p2.eta()) < tolerance &&
      fabs(p1.phi() - p2.phi()) < tolerance &&
      fabs(p1.pt () - p2.pt ()) < 0.1*p1.pt())
    return true;
  return false;
}

inline bool areOverlapping(const reco::Candidate & p1, 
                           const reco::Candidate & p2)
{
  if (p1.numberOfDaughters() == 0 && p2.numberOfDaughters() == 0)
    return areIdentical(p1, p2);
  for (size_t i = 0; i < p1.numberOfDaughters(); i++)
    if (areOverlapping(p2, * p1.daughter(i)))
      return true;
  for (size_t i = 0; i < p2.numberOfDaughters(); i++)
    if (areOverlapping(p1, * p2.daughter(i)))
      return true;
  return false;
}



inline const reco::Candidate * findMother(const reco::Candidate * p)
{
  const reco::Candidate * mother = p ? p->mother(0) : 0;
  if (mother) {
    if (mother->pdgId() == p->pdgId()) return findMother(mother);
    else return mother;
  }
  else return 0;
}

////
/// Return a pointer to product with productName
template <class T, class P>
  const T * getPointerFWLite(const P & ev, std::string productName)
{
  fwlite::Handle<T> handle;
  handle.getByLabel(ev, productName.c_str());
  if (handle.isValid()) {
    return handle.ptr();
  }
  return 0;
}

/// Update event counters if this is a new run/lumi block
template <class P>
void updateEventCounts(P & ev, std::vector<unsigned int> & nEvents, 
                       unsigned int & runNumber, unsigned int & lumiID,
                       std::vector<std::string> & names, bool verbose=false)
{
  if (ev.id().run() != runNumber ||
      ev.id().luminosityBlock() != lumiID) {
    if(verbose) std::cout<<" Changing Lumi block. Run: "<<runNumber<<" Lumi: "<<lumiID<<std::endl;
    fwlite::LuminosityBlock const & lumi = ev.getLuminosityBlock();
    runNumber = ev.id().run();
    lumiID = ev.id().luminosityBlock();
    if (nEvents.size() == 0)
      nEvents.assign(names.size(), 0);
    for (size_t i = 0; i < names.size(); i++) {
      nEvents[i] += (* getPointerFWLite<edm::MergeableCounter>(lumi, names[i])).value;
    }
    if (verbose)
      printf("New luminosity block: %i events present of %i total processed\n",
             nEvents[names.size() - 1], nEvents[0]);
  }
}


#endif // _util_h
