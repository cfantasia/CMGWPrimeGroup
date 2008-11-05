// -*- C++ -*-
//
// Package:    Wprime_muonreco
// Class:      Wprime_muonreco
// 
/**\class Wprime_muonreco Wprime_muonreco.cc WPrime/MyAnalysisMuons/src/Wprime_muonreco.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Alessio Ghezzi, Christos Leonidopoulos, Silvia Goy Lopez
//
//
#include <vector>
#include <string>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/TriggerNames.h"

#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/EDMException.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"


//
// class declaration
//
class TFile;
class TH1F;
class TH2F;

class Wprime_muonreco : public edm::EDAnalyzer 
{
  public:
  explicit Wprime_muonreco(const edm::ParameterSet&);
  ~Wprime_muonreco();
  
  
  private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  edm::InputTag muonTag_;
  edm::InputTag metTag_;
  edm::InputTag HLTTag_; edm::InputTag L1Tag_;
  //      edm::InputTag isoTag_;
  edm::InputTag jetTag_;
  double eJetMin_;
  unsigned int nBinNMu;
  double minNMu;
  double maxNMu;
  unsigned int nBinPtMu;
  double minPtMu;
  double maxPtMu;
  unsigned int nBinEtaMu;
  double minEtaMu;
  double maxEtaMu;
  unsigned int nBinMET;
  double minMET;
  double maxMET;
  unsigned int nBinTMass;
  double minTMass;
  double maxTMass;
  unsigned int nBinAcop;
  double minAcop;
  double maxAcop;
  unsigned int nBinNJets;
  double minNJets;
  double maxNJets;
  
  struct trigEff {
    // key: (muon) trigger name, value: # of accepted events
    std::map<std::string, unsigned> trigger_count;
    // # of events processed
    unsigned Nev; 
    // initialize
    trigEff(){trigger_count.clear(); Nev = 0;}
  };

  typedef std::map<std::string, unsigned>::const_iterator It;

  trigEff genMuTrig; // muon-trigger efficiency wrt generator-muons
  trigEff MuTrig; // muon-trigger efficiency for all processed events

  // muon-detector acceptance
  const float detmu_acceptance;

  // (at least one) generated muon within detector acceptance 
  // (defined as |eta|< detmu_acceptance)
  bool genmu_acceptance;
  // whether event has been accepted by single non-isolated muon (L1 and HL) trigger
  bool muL1_acceptance;
  bool muHLT_acceptance;

  // "golden" single-muon trigger names
  std::string muHLT_20x;
  std::string muHLT_21x;
  std::string muL1;

  // whether this is real-data
  bool realData;

  edm::TriggerNames triggerNames;

  // generator-level muons
  std::vector<HepMC::GenParticle> gen_muons;

  // software versions used to produce HLT and RECO
  std::string HLTversion;
  std::string RECOversion;
  static const std::string INVALID_RELEASE;

  static bool is21x(const std::string & release_string);
  static bool is20x(const std::string & release_string);

  //    Histograms   
  
  TH1F *hPtGen,*hPtRecoOverPtGen;
  TH2F *hPtRecoVsPtGen;
  TH1F *hPtMuUnMatched;
  
  TH1F *hNMu;
  TH1F *hPtMu;
  TH1F *hEtaMu;
  TH1F *hMET;
  TH1F *hTMass;
  TH1F *hAcop;
  TH1F *hNjets;

  TH1F *h_mcmu_pt, *h_mcmu_pt_hlt;
  TH1F *h_mcmu_eta, *h_mcmu_eta_hlt, *h_mcmu_eta_l1;

  //    Root output file
  std::string outputFileName;
  TFile* outputRootFile;

  edm::Handle<reco::TrackCollection> muonCollection;
  
  // # of trigger paths in HLT configuration
  unsigned N_triggers;

  // # of reconstructed muons per event
  unsigned N_muons;

  // initialize run info
  void init_run();

  // initialize event info
  void init_event();

  // initialize histograms
  void init_histograms();

  // initialize trigger structure
  void init_trigger(const edm::Handle<edm::TriggerResults> & hltresults);
  
  // print summary info over full job
  void printSummary() const;
  // print summary info for real
  void printSummary2(const trigEff & trig, const std::string & description) const;

  // get the generator info, populate gen_muons, set genmu_acceptance flag
  void getGenMuons(const edm::Event & iEvent);

  // get trigger info, update muTrig/genMuTrig, set muL1/HLT_acceptance flag
  void getTriggers(const edm::Event & iEvent);

  double met_x, met_y, met;

  // get Calo-MET, initialize MET
  void getCaloMET(const edm::Event & iEvent);

  // get Jets
  void getJets(const edm::Event & iEvent);

  // get muons, update MET
  void getMuons(const edm::Event & iEvent);

  // do MC matching
  void doMCmatching();
};
