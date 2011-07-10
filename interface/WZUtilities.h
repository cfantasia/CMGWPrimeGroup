#ifndef UserAnalysis_WZAnalysis_WZUtilities_H
#define UserAnalysis_WZAnalysis_WZUtilities_H

#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <TROOT.h>
#include <TH1F.h>
#include <TVector3.h>
#include <TFile.h>
#include <TSystem.h>

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"

#include "UserCode/CMGWPrimeGroup/interface/TeVMuon.h"
#include "UserCode/CMGWPrimeGroup/interface/NewElec.h"

typedef edm::ParameterSet PSet;

/// Functions to get products from the event
template <class T, class P>
  T getProduct(const P & event, std::string productName) {
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


//////// Classes /////////////////////////////////////////////////////////////



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
    void loadFromPset(PSet params, std::string param, bool shouldSet = true) {
    C defaultValue = 0;
    if (!this->cutOnMin(param))
      defaultValue = std::numeric_limits<C>::max();
    C val = params.getUntrackedParameter<C>(param, defaultValue);
    this->push_back(param, val);
    this->set(param, shouldSet);
  }
  template<class C>
    void setPassCut(std::string param, C value, pat::strbitset & ret) {
    bool useMin = this->cutOnMin(param);
    if (this->ignoreCut(param) || 
        ( useMin && value >= this->cut(param, C())) ||
        (!useMin && value <= this->cut(param, C())) )
      this->passCut(ret, param);
  }
};

/// Selector for electrons (either barrel or endcap) based on params
class ElectronSelectorBase : public MinMaxSelector<pat::Electron> {
public:
  ElectronSelectorBase() {}
  ElectronSelectorBase(PSet const params) {
    // Set the last parameter to false to turn off the cut
    loadFromPset<double>(params, "minPt", true);
    loadFromPset<double>(params, "minConv", true);
    loadFromPset<double>(params, "maxSigmaIEtaIEta", true);
    loadFromPset<double>(params, "maxDeltaEta", true);
    loadFromPset<double>(params, "maxDeltaPhi", true);
    loadFromPset<double>(params, "maxCombRelIso", true);
    loadFromPset<int>(params, "maxMissingHits", true);
  }
  virtual bool operator()(const pat::Electron & p, pat::strbitset & ret) {
    return (*this)(p,ret,0.);
  }
  virtual bool operator()(const pat::Electron & p, pat::strbitset & ret, const float pu) {
    ret.set(false);
    setPassCut("minPt", p.pt(), ret);
    if(ignoreCut("minConv") || 
       fabs(p.convDist()) >= cut("minConv", double()) || fabs(p.convDcot()) >= cut("minConv", double()))
      passCut(ret, "minConv");
    setPassCut("maxSigmaIEtaIEta", p.sigmaIetaIeta(), ret);
    setPassCut("maxDeltaEta", fabs(p.deltaEtaSuperClusterTrackAtVtx()), ret);
    setPassCut("maxDeltaPhi", fabs(p.deltaPhiSuperClusterTrackAtVtx()), ret);   
    setPassCut("maxCombRelIso",CalcCombRelIso(p, pu), ret);
    setPassCut("maxMissingHits", 
               p.gsfTrack()->trackerExpectedHitsInner().numberOfHits(), ret);
    setIgnored(ret);
    return (bool) ret;
  }
};



/// A wrapper to handle the barrel/endcap split for electrons
class ElectronSelector {
 public:
  ElectronSelector() {}
  ElectronSelector(PSet pset, std::string selectorName) {
    PSet const params = pset.getParameter<PSet>(selectorName);
    barrelSelector_ = ElectronSelectorBase(params.getParameter<PSet>("barrel"));
    endcapSelector_ = ElectronSelectorBase(params.getParameter<PSet>("endcap"));
  }
  bool operator()(const pat::Electron & p, pat::strbitset & ret, const float pu=0.) {
    if     (p.isEB()) return barrelSelector_(p, ret, pu);
    else if(p.isEE()) return endcapSelector_(p, ret, pu);
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
  MuonSelector(PSet pset, std::string selectorName) {
    PSet const params = pset.getParameter<PSet>(selectorName);
    // Set the last parameter to false to turn off the cut
    loadFromPset<double>(params, "minPt", true);
    loadFromPset<double>(params, "maxEta", true);
    loadFromPset<double>(params, "maxDxy", true);
    loadFromPset<double>(params, "maxNormalizedChi2", true);
    loadFromPset<double>(params, "maxIso", true);
    loadFromPset<double>(params, "maxIso03", true);
    loadFromPset<int>(params, "minIsGlobal", true);
    loadFromPset<int>(params, "minIsTracker", true);
    loadFromPset<int>(params, "minNTrackerHits", true);
    loadFromPset<int>(params, "minNPixelHits", true);
    loadFromPset<int>(params, "minNMuonHits", true);
    loadFromPset<int>(params, "minNMatches", true);

  }
  virtual bool operator()(const TeVMuon & p, pat::strbitset & ret) {
    return (*this)(p,ret,0.);
  }
  bool operator()(const TeVMuon & p, pat::strbitset & ret, const float pu) {
    ret.set(false);
    reco::TrackRef global = p.globalTrack();
    setPassCut("minPt", p.pt(), ret);
    setPassCut("maxEta", p.eta(), ret);
    setPassCut("maxDxy", fabs(p.dB()), ret);
    setPassCut("maxNormalizedChi2", global.isNull() ? 0. : 
               global->normalizedChi2(), ret);
    setPassCut("maxIso", p.combRelIsolation(), ret);
    setPassCut("maxIso03", p.combRelIsolation03(pu), ret);
    setPassCut("minIsGlobal", p.isGlobalMuon(), ret);
    setPassCut("minIsTracker", p.isTrackerMuon(), ret);
    setPassCut("minNTrackerHits", global.isNull() ? 0 :
               global->hitPattern().numberOfValidTrackerHits(), ret);
    setPassCut("minNPixelHits", global.isNull() ? 0 :
               global->hitPattern().numberOfValidPixelHits(), ret);
    setPassCut("minNMuonHits", global.isNull() ? 0 :
               global->hitPattern().numberOfValidMuonHits(), ret);
    setPassCut("minNMatches", p.numberOfMatches(), ret);
    setIgnored(ret);
    return (bool) ret;
  }
};



/// Selector for jets based on params
class JetSelector : public MinMaxSelector<pat::Jet> {
public:
  JetSelector() {}
  JetSelector(PSet pset, std::string selectorName) {
    PSet const params = pset.getParameter<PSet>(selectorName);
    // Set the last parameter to false to turn off the cut
    loadFromPset<double>(params, "minPt", true);
    loadFromPset<double>(params, "maxEta", true);
  }
  virtual bool operator()(const pat::Jet & p, pat::strbitset & ret) {
    ret.set(false);
    setPassCut("minPt", p.pt(), ret);
    setPassCut("maxEta", p.eta(), ret);
    setIgnored(ret);
    return (bool) ret;
  }
};

//////// Functions ///////////////////////////////////////////////////////////

/// Return the DBS dataset name stored in the patTuple
template <class P>
std::string getDatasetName(const P & event, const std::string datasetName) {
  if (datasetName == "")
    return getProduct<std::string>(event, "wzPreselectionProducer:datasetName", 
                              "Not found");
  return datasetName;
}




/*
TVector2
Heep::getPtDiff(){
  
}

TVector2
adjustPt(const ElectronV & electrons){
  TVector2 diff(0.,0.);
  for (ElectronV::const_iterator i = electrons.begin(); 
       i != electrons.end(); ++i){
//Cory: Under development    
    diff += i.getPtDiff();
  }
  return diff;
}
*/


#endif
