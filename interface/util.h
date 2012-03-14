#ifndef _util_h
#define _util_h

#include <vector>
#include <string>
#include <map>

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

#include "PhysicsTools/SelectorUtils/interface/Selector.h"
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
typedef std::vector<pat::MET     > METV;
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

////Selectors//////////

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


//////// Selectors ///////////////////////////////////////////////////////////
/// Class derived from Selector that implements some convenience functions
template<class T>
class MinMaxSelector : public Selector<T> {
 public:
  bool cutOnMin(std::string param) {
    bool useMin = param.find("min") != std::string::npos;
    bool useMax = param.find("max") != std::string::npos;
    if (!useMin && !useMax) {
      std::cout << param << " has neither 'min' nor 'max' in its name!" << std::endl;
      exit(1);
    }
    return useMin;
  }
  template<class C>
    void loadFromPset(Pset params, std::string param, bool shouldset = true) {
    C defaultValue = 0;
    if (!this->cutOnMin(param))
      defaultValue = std::numeric_limits<C>::max();
    C val = params.getUntrackedParameter<C>(param, defaultValue);
    if(!params.exists(param)){
      shouldset = false;//If you don't list it, I'll ignore it
      std::cout<<" You didn't specify param "<<param<<", so it will be ignored\n";
    }
    this->push_back(param, val);
    this->set(param, shouldset);
  }
  template<class C>
    void setpassCut(std::string param, C value, pat::strbitset & ret) {
    bool useMin = this->cutOnMin(param);
    if (this->ignoreCut(param) || 
        ( useMin && value >= this->cut(param, C())) ||
        (!useMin && value <= this->cut(param, C())) )
      this->passCut(ret, param);
  }

  pat::strbitset bitmask;
};

/// Selector for electrons (either barrel or endcap) based on params
class ElectronSelectorBase : public MinMaxSelector<pat::Electron> {
public:
  ElectronSelectorBase() {}
  ElectronSelectorBase(Pset const params) {
    // set the last parameter to false to turn off the cut
    loadFromPset<double>(params, "minPt", true);
    loadFromPset<double>(params, "minConv", true);
    loadFromPset<double>(params, "maxSigmaIEtaIEta", true);
    loadFromPset<double>(params, "maxDeltaEta", true);
    loadFromPset<double>(params, "maxDeltaPhi", true);
    loadFromPset<double>(params, "maxCombRelIso", true);
    loadFromPset<int>(params, "maxMissingHits", true);

    bitmask = getBitTemplate();
  }
  virtual bool operator()(const pat::Electron & p, pat::strbitset & ret) {
    bool val = (*this)(p,0.);
    ret = bitmask;
    return val;
  }
  virtual bool operator()(const pat::Electron & p, const float pu) {
    bitmask.set(false);
    setpassCut("minPt", p.pt(), bitmask);
    if(ignoreCut("minConv") || 
       fabs(p.convDist()) >= cut("minConv", double()) || fabs(p.convDcot()) >= cut("minConv", double()))
      passCut(bitmask, "minConv");
    setpassCut("maxSigmaIEtaIEta", p.sigmaIetaIeta(), bitmask);
    setpassCut("maxDeltaEta", fabs(p.deltaEtaSuperClusterTrackAtVtx()), bitmask);
    setpassCut("maxDeltaPhi", fabs(p.deltaPhiSuperClusterTrackAtVtx()), bitmask);   
    setpassCut("maxCombRelIso",calcCombRelIso(p, pu), bitmask);
    setpassCut("maxMissingHits", 
               p.gsfTrack()->trackerExpectedHitsInner().numberOfHits(), bitmask);
    setIgnored(bitmask);
    return (bool) bitmask;
  }
  virtual bool operator()(const heep::Ele & p, pat::strbitset & ret) {
    bool val = (*this)(p,0.);
    ret = bitmask;
    return val;
  }
  virtual bool operator()(const heep::Ele & p, const float pu) {
    bitmask.set(false);
    const pat::Electron & e = p.patEle();
    setpassCut("minPt", p.et(), bitmask);
    if(ignoreCut("minConv") || 
       fabs(e.convDist()) >= cut("minConv", double()) || fabs(e.convDcot()) >= cut("minConv", double()))
      passCut(bitmask, "minConv");
    setpassCut("maxSigmaIEtaIEta", e.sigmaIetaIeta(), bitmask);
    setpassCut("maxDeltaEta", fabs(e.deltaEtaSuperClusterTrackAtVtx()), bitmask);
    setpassCut("maxDeltaPhi", fabs(e.deltaPhiSuperClusterTrackAtVtx()), bitmask);   
    setpassCut("maxCombRelIso",calcCombRelIso(e, pu), bitmask);
    setpassCut("maxMissingHits", 
               e.gsfTrack()->trackerExpectedHitsInner().numberOfHits(), bitmask);
    setIgnored(bitmask);
    return (bool) bitmask;
  }
};



/// A wrapper to handle the barrel/endcap split for electrons
class ElectronSelector {
 public:
  ElectronSelector() {}
  ElectronSelector(Pset pset, std::string selectorName) {
    Pset const params = pset.getParameter<Pset>(selectorName);
    barrelSelector_ = ElectronSelectorBase(params.getParameter<Pset>("barrel"));
    endcapSelector_ = ElectronSelectorBase(params.getParameter<Pset>("endcap"));
  }
  bool operator()(const pat::Electron & p, const float pu=0.) {
    if     (p.isEB()) return barrelSelector_(p, pu);
    else if(p.isEE()) return endcapSelector_(p, pu);
    return false;
  }
  bool operator()(const heep::Ele & p, const float pu=0.) {
    if     (p.isEB()) return barrelSelector_(p, pu);
    else if(p.isEE()) return endcapSelector_(p, pu);
    return false;
  }
  pat::strbitset getBitTemplate() { return barrelSelector_.getBitTemplate(); }
 private:
  ElectronSelectorBase barrelSelector_;
  ElectronSelectorBase endcapSelector_;
};



/// Selector for muons based on params
class MuonSelector : public MinMaxSelector<TeVMuon> {
public:
  MuonSelector() {}
  MuonSelector(Pset pset, std::string selectorName) {
    Pset const params = pset.getParameter<Pset>(selectorName);
    // set the last parameter to false to turn off the cut
    loadFromPset<double>(params, "minPt", true);
    loadFromPset<double>(params, "maxEta", true);
    loadFromPset<double>(params, "maxDxy", true);
    loadFromPset<double>(params, "maxIso", true);
    loadFromPset<double>(params, "maxIso03", true);
    loadFromPset<int>(params, "minIsGlobal", true);
    loadFromPset<int>(params, "minIsTracker", true);
    loadFromPset<int>(params, "minNMatches", true);
    loadFromPset<int>(params, "minNMatchedStations", true);
    loadFromPset<double>(params, "maxNormalizedChi2", true);
    loadFromPset<int>(params, "minNTrackerHits", true);
    loadFromPset<int>(params, "minNPixelHits", true);
    loadFromPset<int>(params, "minNMuonHits", true);
    loadFromPset<int>(params, "minNTrackerLayers", true);
    loadFromPset<double>(params, "minTrackerValidFrac", true);

    bitmask = getBitTemplate();
  }
  virtual bool operator()(const TeVMuon & p, pat::strbitset & ret) {
    bool val = (*this)(p,0.);
    ret = bitmask;
    return val;
  }
  bool operator()(const TeVMuon & p, const float pu=0) {
    bitmask.set(false);
    setpassCut("minPt", p.pt(), bitmask);
    setpassCut("maxEta", fabs(p.eta()), bitmask);
    setpassCut("maxDxy", fabs(p.dB()), bitmask);
    setpassCut("maxIso", p.combRelIsolation(), bitmask);
    setpassCut("maxIso03", p.combRelIsolation03(pu), bitmask);
    setpassCut("minIsGlobal", p.isGlobalMuon(), bitmask);
    setpassCut("minIsTracker", p.isTrackerMuon(), bitmask);
    setpassCut("minNMatches", p.numberOfMatches(), bitmask);
    setpassCut("minNMatchedStations", p.numberOfMatchedStations(), bitmask);
    
    reco::TrackRef global = p.globalTrack();
    if(!global.isNull()){
      const reco::HitPattern& gtHP = global->hitPattern();
      setpassCut("maxNormalizedChi2", global->normalizedChi2(), bitmask);
      setpassCut("minNTrackerHits", gtHP.numberOfValidTrackerHits(), bitmask);
      setpassCut("minNPixelHits", gtHP.numberOfValidPixelHits(), bitmask);
      setpassCut("minNMuonHits", gtHP.numberOfValidMuonHits(), bitmask);
      setpassCut("minNTrackerLayers", gtHP.trackerLayersWithMeasurement(), bitmask);
      setpassCut("minTrackerValidFrac", global->validFraction(), bitmask);
    }
    setIgnored(bitmask);
    return (bool) bitmask;
  }
};



/// Selector for jets based on params
class JetSelector : public MinMaxSelector<pat::Jet> {
public:
  JetSelector() {}
  JetSelector(Pset pset, std::string selectorName) {
    Pset const params = pset.getParameter<Pset>(selectorName);
    // set the last parameter to false to turn off the cut
    loadFromPset<double>(params, "minPt", true);
    loadFromPset<double>(params, "maxEta", true);
    loadFromPset<double>(params, "maxNHF", true);
    loadFromPset<double>(params, "maxNEF", true);
    loadFromPset<int>(params, "minNDaughters", true);
    loadFromPset<double>(params, "minCHF", true);
    loadFromPset<double>(params, "maxCEF", true);
    loadFromPset<int>(params, "minCMult", true);

    bitmask = getBitTemplate();
  }
  virtual bool operator()(const pat::Jet & p, pat::strbitset & ret) {
    bool val = (*this)(p);
    ret = bitmask;
    return val;
  }
  virtual bool operator()(const pat::Jet & p) {
    bitmask.set(false);
    bool inTracking = fabs(p.eta()) < 2.4;
    if(ignoreCut("minPt")  || p.pt() > cut("minPt", double())) passCut(bitmask, "minPt");
    if(ignoreCut("maxEta") || fabs(p.eta()) < cut("maxEta", double())) passCut(bitmask, "maxEta");
    if(ignoreCut("maxNHF") || p.neutralHadronEnergyFraction() < cut("maxNHF", double())) passCut(bitmask, "maxNHF");
    if(ignoreCut("maxNEF") || p.neutralEmEnergyFraction()     < cut("maxNEF", double())) passCut(bitmask, "maxNEF");
    if(ignoreCut("minNDaughters") || (int)p.numberOfDaughters() > cut("minNDaughters", int())) passCut(bitmask, "minNDaughters");
    //Below are only used for fabs(eta) < 2.4 b/c tracking needed
    if(ignoreCut("minCHF")   || !inTracking || p.chargedHadronEnergyFraction() > cut("minCHF",   double())) passCut(bitmask, "minCHF");
    if(ignoreCut("maxCEF")   || !inTracking || p.chargedEmEnergyFraction()     < cut("maxCEF",   double())) passCut(bitmask, "maxCEF");
    if(ignoreCut("minCMult") || !inTracking || (int)p.chargedMultiplicity()    > cut("minCMult",    int())) passCut(bitmask, "minCMult");
    setIgnored(bitmask);
    return (bool) bitmask;
  }
  
};

class PhotonSelector : public MinMaxSelector<pat::Photon> {
public:
  PhotonSelector(){}
  PhotonSelector(Pset pset, std::string selectorName){
    Pset const params = pset.getParameter<Pset>(selectorName);
    // set the last parameter to false to turn off the cut
    loadFromPset<double>(params, "minPt", true);
    loadFromPset<double>(params, "maxEta", true);
    loadFromPset<double>(params, "maxECalIso", true);
    loadFromPset<double>(params, "maxHCalIso", true);
    loadFromPset<double>(params, "maxTrkIso", true);
    loadFromPset<double>(params, "maxHoE", true);
    loadFromPset<double>(params, "maxSigmaee", true);
    loadFromPset<bool>  (params, "minHasSeed", true);

    bitmask = getBitTemplate();
  }
  
  virtual bool operator()(const pat::Photon & p, pat::strbitset & ret) {
    bool val = (*this)(p);
    ret = bitmask;
    return val;
  }
  virtual bool operator()(const pat::Photon & p) {
    bitmask.set(false);
    if(ignoreCut("minPt")       || p.pt() > cut("minPt", double())) passCut(bitmask, "minPt");
    if(ignoreCut("maxEta")      || fabs(p.superCluster()->eta()) < cut("maxEta", double())) passCut(bitmask, "maxEta");
    if(ignoreCut("maxECalIso")  || p.ecalRecHitSumEtConeDR04() < cut("maxECalIso", double())) passCut(bitmask, "maxECalIso");
    if(ignoreCut("maxHCalIso")  || p.hcalTowerSumEtConeDR04() < cut("maxHCalIso", double())) passCut(bitmask, "maxHCalIso");
    if(ignoreCut("maxTrkIso")   || p.trkSumPtHollowConeDR04() < cut("maxTrkIso", double())) passCut(bitmask, "maxTrkIso");
    if(ignoreCut("maxHoE")      || p.hadronicOverEm() < cut("maxHoE", double())) passCut(bitmask, "maxHoE");
    if(ignoreCut("maxSigmaee")  || p.sigmaIetaIeta() < cut("maxSigmaee", double())) passCut(bitmask, "maxSigmaee");
    if(ignoreCut("minHasSeed")  || p.hasPixelSeed() >= cut("minHasSeed", bool())) passCut(bitmask, "minHasSeed");
    setIgnored(bitmask);
    return (bool) bitmask;
  }
 
};

/////////////////////
////Others Fns///////
/////////////////////

template<class T1, class T2>
  bool Overlap(const T1 & p, const std::vector<T2>& vec, const float minDR=0.01, const size_t maxToCheck=0 ){
  uint max = maxToCheck == 0 ? vec.size() : min(vec.size(), maxToCheck);
  for(uint i=0; i<max; ++i){
    if(reco::deltaR(p, vec[i]) < minDR) return true;
  }
  return false;
}

#endif // _util_h
