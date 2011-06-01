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

using namespace std;

typedef edm::ParameterSet PSet;

/// Functions to get products from the event
template <class T, class P>
T getProduct(const P & event, string productName) {
  edm::Handle<T> handle;
  event.getByLabel(edm::InputTag(productName), handle);
  return * handle;
}

template <class T, class P>
T getProduct(const P & event, string productName, T defaultVal) {
  edm::Handle<T> handle;
  event.getByLabel(edm::InputTag(productName), handle);
  if (handle.isValid()) return * handle;
  return defaultVal;
}

template <class T, class P>
const T * getPointer(const P & event, string productName) {
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
  bool cutOnMin(string param) {
    bool useMin = param.find("min") != string::npos;
    bool useMax = param.find("max") != string::npos;
    if (!useMin && !useMax) {
      cout << param << " has neither 'min' nor 'max' in its name!" << endl;
      exit(1);
    }
    return useMin;
  }
  template<class C>
  void loadFromPset(PSet params, string param, bool shouldSet = true) {
    C defaultValue = 0;
    if (!this->cutOnMin(param))
      defaultValue = numeric_limits<C>::max();
    C val = params.getUntrackedParameter<C>(param, defaultValue);
    this->push_back(param, val);
    this->set(param, shouldSet);
  }
  template<class C>
  void setPassCut(string param, C value, pat::strbitset & ret) {
    bool useMin = this->cutOnMin(param);
    bool passMin = useMin && value >= this->cut(param, C());
    bool passMax = !useMin && value <= this->cut(param, C());
    if (passMin || passMax || this->ignoreCut(param))
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
    loadFromPset<double>(params, "maxDeltaEta", true);
    loadFromPset<double>(params, "maxDeltaPhi", true);
    loadFromPset<double>(params, "maxSigmaEtaEta", true);
    loadFromPset<double>(params, "maxSigmaIEtaIEta", true);
    loadFromPset<double>(params, "minEoverP", true);
    loadFromPset<double>(params, "maxHoverE", true);
    loadFromPset<double>(params, "maxCaloIso", true);
    loadFromPset<double>(params, "maxTrackIso", true);
    loadFromPset<double>(params, "maxSwissCross", true);
    loadFromPset<int>(params, "maxMissingHits", true);
    loadFromPset<int>(params, "minSimpleEleId60relIso", true);
    loadFromPset<int>(params, "minSimpleEleId70relIso", true);
    loadFromPset<int>(params, "minSimpleEleId80relIso", true);
    loadFromPset<int>(params, "minSimpleEleId85relIso", true);
    loadFromPset<int>(params, "minSimpleEleId90relIso", true);
    loadFromPset<int>(params, "minSimpleEleId95relIso", true);
  }
  virtual bool operator()(const pat::Electron & p, pat::strbitset & ret) {
    ret.set(false);
    setPassCut("minPt", p.pt(), ret);
    setPassCut("maxDeltaEta", p.deltaEtaSuperClusterTrackAtVtx(), ret);
    setPassCut("maxDeltaPhi", p.deltaPhiSuperClusterTrackAtVtx(), ret);
    setPassCut("maxSigmaEtaEta", p.scSigmaEtaEta(), ret);
    setPassCut("maxSigmaIEtaIEta", p.sigmaIetaIeta(), ret);
    setPassCut("minEoverP", p.eSuperClusterOverP(), ret);
    setPassCut("maxHoverE", p.hadronicOverEm(), ret);
    setPassCut("maxCaloIso", p.userFloat("relCaloIso"), ret);
    setPassCut("maxTrackIso", p.userFloat("relTrackIso"), ret);
    setPassCut("maxSwissCross", p.userFloat("swissCross"), ret);
    setPassCut("maxMissingHits", 
               p.gsfTrack()->trackerExpectedHitsInner().numberOfHits(), ret);
    setPassCut("minSimpleEleId60relIso", p.electronID("simpleEleId60relIso"), ret);
    setPassCut("minSimpleEleId70relIso", p.electronID("simpleEleId70relIso"), ret);
    setPassCut("minSimpleEleId80relIso", p.electronID("simpleEleId80relIso"), ret);
    setPassCut("minSimpleEleId85relIso", p.electronID("simpleEleId85relIso"), ret);
    setPassCut("minSimpleEleId90relIso", p.electronID("simpleEleId90relIso"), ret);
    setPassCut("minSimpleEleId95relIso", p.electronID("simpleEleId95relIso"), ret);
    setIgnored(ret);
    return (bool) ret;
  }
};



/// A wrapper to handle the barrel/endcap split for electrons
class ElectronSelector {
 public:
  ElectronSelector() {}
  ElectronSelector(PSet pset, string selectorName) {
    PSet const params = pset.getParameter<PSet>(selectorName);
    barrelSelector_ = ElectronSelectorBase(params.getParameter<PSet>("barrel"));
    endcapSelector_ = ElectronSelectorBase(params.getParameter<PSet>("endcap"));
  }
  bool operator()(const pat::Electron & p, pat::strbitset & ret) {
    double fabsEta = fabs(p.eta());
    if (fabsEta < 1.479) return barrelSelector_(p, ret);
    if (fabsEta > 1.550) return endcapSelector_(p, ret);
    return false;
  }
  pat::strbitset getBitTemplate() { return barrelSelector_.getBitTemplate(); }
 private:
  ElectronSelectorBase barrelSelector_;
  ElectronSelectorBase endcapSelector_;
};



/// Selector for muons based on params
class MuonSelector : public MinMaxSelector<pat::Muon> {
public:
  MuonSelector() {}
  MuonSelector(PSet pset, string selectorName) {
    PSet const params = pset.getParameter<PSet>(selectorName);
    // Set the last parameter to false to turn off the cut
    loadFromPset<double>(params, "minPt", true);
    loadFromPset<double>(params, "maxSip", true);
    loadFromPset<double>(params, "maxIso", true);
    loadFromPset<double>(params, "maxNormalizedChi2", true);
    loadFromPset<int>(params, "minNTrackerHits", true);
    loadFromPset<int>(params, "minNPixelHits", true);
    loadFromPset<int>(params, "minNMuonHits", true);
    loadFromPset<int>(params, "minNMatches", true);
  }
  virtual bool operator()(const pat::Muon & p, pat::strbitset & ret) {
    ret.set(false);
    reco::TrackRef inner = p.innerTrack();
    reco::TrackRef global = p.globalTrack();
    setPassCut("minPt", p.pt(), ret);
    setPassCut("maxSip", p.userFloat("sip"), ret);
    setPassCut("maxIso", p.userFloat("relIso"), ret);
    setPassCut("maxNormalizedChi2", global.isNull() ? 0. : 
               global->normalizedChi2(), ret);
    setPassCut("minNTrackerHits", inner.isNull() ? 0 :
               inner->hitPattern().numberOfValidTrackerHits(), ret);
    setPassCut("minNPixelHits", inner.isNull() ? 0 :
               inner->hitPattern().numberOfValidPixelHits(), ret);
    setPassCut("minNMuonHits", global.isNull() ? 0 :
               global->hitPattern().numberOfValidMuonHits(), ret);
    setPassCut("minNMatches", p.numberOfMatches(), ret);
    setIgnored(ret);
    return (bool) ret;
  }
};




//////// Functions ///////////////////////////////////////////////////////////

/// Return the DBS dataset name stored in the patTuple
template <class P>
string getDatasetName(const P & event, const string datasetName) {
  if (datasetName == "")
    return getProduct<string>(event, "wzPreselectionProducer:datasetName", 
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
