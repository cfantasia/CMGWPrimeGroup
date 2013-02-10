#ifndef _util_h
#define _util_h

#include <vector>
#include <string>
#include <map>

#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include <DataFormats/PatCandidates/interface/PFParticle.h>
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

////
#include <memory>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <TROOT.h>
#include <TH1F.h>
#include <TVector3.h>
#include <TFile.h>
#include <TSystem.h>

#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"

#include "UserCode/CMGWPrimeGroup/interface/TeVMuon.h"
#include "UserCode/CMGWPrimeGroup/interface/NewElec.h"
#include "SHarper/HEEPAnalyzer/interface/HEEPEle.h"


typedef unsigned int uint;
typedef std::vector<int> vint;
typedef std::vector<std::string> vstring;
typedef math::XYZPoint Point;
typedef math::XYZTLorentzVector LorentzVector;
typedef edm::ParameterSet Pset;
typedef std::vector<Pset> VPset;
typedef edm::MergeableCounter counter;

typedef std::vector<pat::Electron> PatElectronV;
typedef std::vector<pat::Muon    > PatMuonV;
typedef std::vector<pat::PFParticle > PFCandidateV;
typedef std::vector<heep::Ele > ElectronV;
typedef std::vector<TeVMuon  > MuonV;

typedef std::vector<pat::Jet     > JetV;
//typedef std::vector<pat::MET     > METV;
typedef std::vector<reco::PFMET   > METV;
typedef std::vector<reco::Track  > TrackV;
typedef std::vector<edm::InputTag> VInputTag;
typedef std::vector<reco::Candidate> CandV;
typedef std::vector<reco::GenParticle> GenParticleV;

typedef edm::Handle<PatElectronV > PatElectronVH;
typedef edm::Handle<PatMuonV > PatMuonVH;
typedef edm::Handle<JetV > JetVH;
typedef edm::Handle<PFCandidateV > PFCandidateVH;
typedef edm::Handle<METV > METVH;


//////

namespace wprime{
  struct FilterEff{
    // # of (unweighted!) events surviving after each selection cut
    int Nsurv_evt_cut;
    // # of (weighted!) events surviving after each selection cut
    float Nsurv_evt_cut_w;
    // efficiency for each selection cut
    float eff;
    // efficiency uncertainty
    float deff;
    // absolute efficiency for all selection cuts
    float eff_abs;
    // absolute efficiency uncertainty
    float deff_abs;

    FilterEff(){Nsurv_evt_cut = 0; Nsurv_evt_cut_w = eff = deff = eff_abs = deff_abs = 0.0;}
  };

  // key: samplename, value: vector<FilterEff> (ie. statistics for selection steps)
  typedef std::vector<FilterEff> EffV;
  typedef std::map<std::string, EffV > SampleStat;

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
    float signalMass; // in <x>.<y> TeV for signal samples (e.g. 1.2), -1 otherwise (default)
    bool isSignal() const {return signalMass > 0;}
    int splitInto;
    EffV results;
    //
    InputFile()
    {
      x_sect = signalMass = -1; Nprod_evt = Nact_evt = -1; weight = 0; splitInto = 1;
      samplename = description = INVALID;
      subdir = "";
    }
    void checkFile()
    {
      if(samplename.find("data") == std::string::npos){
        //Only require these if running on MC
        assert(x_sect > 0); 
        assert(Nprod_evt > 0); 
      }
      assert(weight > 0);
      assert(splitInto > 0);
      assert(pathnames.size()); 
      assert(samplename != INVALID);
      if(description.find(INVALID) == 0) description.replace(0, INVALID.size(),samplename); 
      assert(description != INVALID);
    }

  };

  struct AnalysisCut{
    std::string Name;
    std::string Desc;
    typedef boost::function<bool()> fnCut;
    std::map<std::string,fnCut > mFnPtrs_;
    std::vector<fnCut > CutFn;
  };
  typedef std::map<std::string, AnalysisCut> mCuts;

  const float MUON_MASS = 0.105658366;      // GeV
  const float ELECTRON_MASS = 0.000511;     // GeV
}

//// Global constants ////////////////////////////////////////////////////////

const float ZMASS = 91.188;
const float WMASS = 80.398;
const float VMASS = (ZMASS+WMASS)/2.;
const float TMASS = 172.9;
const float MAX_MBL = 153.1; //sqrt(TMASS*TMASS - WMASS*WMASS)

const int PDG_ID_ELEC = 11;
const int PDG_ID_ELECNEU = 12;
const int PDG_ID_MUON = 13;
const int PDG_ID_MUONNEU = 14;
const int PDG_ID_TAU = 15;
const int PDG_ID_TAUNEU = 16;

const int PDG_ID_W = 24;
const int PDG_ID_Z = 23;
const int PDG_ID_WPRIME = 34;

const float PI    = 2.0 * TMath::ACos(0.);
//const float TWOPI = 2.0 * PI;
const float NOCUT = 9e9;

//////// Comparators (for sorting vectors) ///////////////////////////////////
struct closestToZMass {
  bool operator() (const reco::Candidate & a, const reco::Candidate & b) {
    return fabs(ZMASS - a.mass()) < fabs(ZMASS - b.mass());
  }
};
struct closestToWTMass {
  bool operator() (const reco::Candidate & a, const reco::Candidate & b) {
    return fabs(WMASS - a.mt()) < fabs(WMASS - b.mt());
  }
};
struct closestToVMass {
  bool operator() (const reco::Candidate & a, const reco::Candidate & b) {
    return fabs(VMASS - a.mass()) < fabs(VMASS - b.mass());
  }
};

template<class C, class D>
struct highestPt {
  bool operator() (const C * a, const D * b) {
    return a->pt() > b->pt();
  }
  bool operator() (const C & a, const D & b) {
    return a.pt() > b.pt();
  }
};
struct highestPtLepton {
  bool operator() (const reco::CompositeCandidate a, 
                   const reco::CompositeCandidate b) {
    return a.daughter(0)->pt() > b.daughter(0)->pt();
  }
};
struct highestMuonPt {                                                                                                                                     
  bool operator() (const TeVMuon & a, const TeVMuon & b){                                                                                  
    return a.pt() > b.pt();                                                                                                                                
  }                                                                                                                                                        
};

struct highestElectronPt{                                                                                                                                   
  bool operator() (const heep::Ele & a, const heep::Ele & b){                                                                                              
    return a.patEle().pt() > b.patEle().pt();                                                                                                           
  }                                                                                                                                                        
};

struct highestJetPt {
  bool operator() (const pat::Jet & a, const pat::Jet & b){
    return a.pt() > b.pt();
  }
};

inline bool areIdentical(const heep::Ele & a, const heep::Ele & b)
{
  return areIdentical(a.patEle(), b.patEle());
}

template<class C, class D>
inline bool areIdentical(const C & p1, const D & p2)
{
  float tolerance = 0.0001;
  if (p1.pdgId() == p2.pdgId() &&
      fabs(p1.eta() - p2.eta()) < tolerance &&
      fabs(p1.phi() - p2.phi()) < tolerance &&
      fabs(p1.pt () - p2.pt ()) < 0.1*p1.pt())
    return true;
  return false;
}

template<class C, class D>
inline bool areOverlapping(const C & p1, const D & p2)
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
void updateEventcounts(P & ev, std::vector<unsigned int> & nEvents, 
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
/*
/// Return the DBS dataset name stored in the patTuple
template <class P>
std::string getDataSetName(const P & event, const std::string dataSetName) {
  if (dataSetName == "")
    return getProduct<std::string>(event, "wzPreselectionProducer:dataSetName", 
                              "Not found");
  return dataSetName;
}
*/
template <class P>
bool passTriggerMatch(const P & p, float cut, std::vector<std::string>& triggers){
  for (size_t i=0; i < triggers.size(); ++i){
    if (p.triggerObjectMatchesByPath(triggers[i], true, false).size() > 0){
      const pat::TriggerObjectStandAlone * trigRef = p.triggerObjectMatchByPath(triggers[i], true, false);
      if(trigRef->et() > cut) return true;;
    }
  }
  return false;
}


/// Functions to get products from the event
template <class T, class P, class R>
  T getProduct(const P & event, R productName) {
  edm::Handle<T> handle;
  event.getByLabel(edm::InputTag(productName), handle);
  return * handle;
}

template <class T, class P>
  T getProduct(const P & event, std::string productName, T defaultVal) {
  edm::Handle<T> handle;
  event.getByLabel(edm::InputTag(productName), handle);
  if (handle.isValid()) return * handle;
  return defaultVal;
}

template <class T, class P>
  const T * getPointer(const P & event, std::string productName) {
  edm::Handle<T> handle;
  event.getByLabel(edm::InputTag(productName), handle);
  if (handle.isValid()) return handle.product();
  return 0;
}


/////////////////////
////Others Fns///////
/////////////////////

template<class T1, class T2>
  bool Overlap(const T1 & p, const std::vector<T2>& vec, const float minDR=0.01, const size_t maxToCheck=0 ){
  uint max = maxToCheck == 0 ? vec.size() : std::min(vec.size(), maxToCheck);
  for(uint i=0; i<max; ++i){
    if(reco::deltaR(p, vec[i]) < minDR) return true;
  }
  return false;
}

#endif // _util_h
