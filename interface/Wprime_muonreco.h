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

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

//
// class declaration
//
class TFile;
class TH1F;
class TH2F;

namespace Wprime_muonreco_histo
{
  // histogram parameters (energy/pt units: GeV)

  // Muon pt distribution
  static const unsigned  nBinPtMu = 1500;
  static const float minPtMu = 0;
  static const float  maxPtMu = 3000;

  // jet distributions (# of jets and Et)
  static const unsigned nBinNJets = 50;
  static const float minNJets = -0.5;
  static const float maxNJets = 49.5;
  static const unsigned nBinEtJets = 150;
  static const float minEtJets = 0;
  static const float maxEtJets = 300;
  // jet-activity veto (# of jets and Et)
  static const unsigned nBinEtJets_veto = 10;
  static const float minEtJets_veto = 0;
  static const float maxEtJets_veto = 200;
  static const unsigned nBinNJets_veto = 10;
  static const float minNJets_veto = -0.5;
  static const float maxNJets_veto = 9.5;

  // cone size and SumPt for isolation
  static const unsigned nBinCone = 20;
  static const float minCone = 0.05;
  static const float maxCone = 1.05;
  static const unsigned nBinSumPt = 600;
  static const float minSumPt = 0;
  static const float maxSumPt = 600;
}

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
  edm::InputTag tkIsoMapTag_;
  edm::InputTag ecalIsoMapTag_;
  edm::InputTag hcalIsoMapTag_;

  struct muonTrack {
    // reference to standard tracking
    reco::MuonRef mu;
    // reference to default-TeV, 1st-Hit, picky and cocktail/optimized tracking
    reco::TrackRef TeVMuons[4];
  };

  // reference to different tracking algorithms for all good muons
  std::vector<muonTrack> good_muons;

  typedef std::vector<muonTrack>::const_iterator mIt;

  const float eJetMin_;
  unsigned NJetsAboveThres; // # of jets in event above eJetMin_
  
  struct trigEff {
    // key: (muon) trigger name, value: # of accepted events
    std::map<std::string, unsigned> trigger_count;
    // # of events processed
    unsigned Nev; 
    // # of events processed with a single reconstructed muon
    unsigned Nev_1mu; 
    // initialize
    trigEff(){trigger_count.clear(); Nev = Nev_1mu = 0;}
  };

  typedef std::map<std::string, unsigned>::const_iterator tIt;

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
  TH1F *hPtGenJetVeto,*hPtRecoOverPtGenJetVeto;
  TH1F *hPtGenOne, *hPtGenOneMatch, *hPtGenOnePtlt5Match;
  TH1F *h1PtGen1PtReco;
  TH2F *h1PtGen1PtRecoVsPtGen;
  TH2F *hPtRecoVsPtGen;
  TH1F *hPtMuUnMatched;
  TH1F *hPtMuUnMatchedJetVeto;

  TH1F *hNMu;
  TH1F *hPtMu,*hPtMaxMu,*hPtOneMu, *hPtMuOneMatch;
  TH1F *hPtMuJetVeto,*hPtMaxMuJetVeto,*hPtOneMuJetVeto, *hPtMuOneMatchJetVeto;
  TH1F *hHitTrack,*hHitMuon,*hHitComb,*hHitCheckMu;
  TH1F *hEtaMu, *hEtaOnePtlt5Match;
  TH1F *hMET;
  TH1F *hTMass;
  TH1F *hAcop;
  TH1F *hNAlljets;
  TH1F *hNjetsGtEJetMin;
  TH1F *hJetEt;
  TH2F *hEthresNjet,*hEthresNjet_norm;

  TH1F *hPtSumR03,*hPtSumNR03;
  TH2F *hTMass_PtSumR03;
  TH2F *hPtTkIso_ConeSize,*hPtTkIso_mmupt_ConeSize,*hNTkIso_ConeSize,*hEcalIso_ConeSize,*hHcalIso_ConeSize;
  TH2F *hPtTkIsoThresh_ConeSize,*hPtTkIsoThresh_ConeSize_norm; 

  TH1F *hMuChi2,*hTrackerMuChi2,*hSAMuChi2;
  TH1F *hMuNdf,*hTrackerMuNdf,*hSAMuNdf;

  TH1F *hMuChi2UnMatched,*hTrackerMuChi2UnMatched,*hSAMuChi2UnMatched;
  TH1F *hMuNdfUnMatched,*hTrackerMuNdfUnMatched,*hSAMuNdfUnMatched;

  TH1F *hPtTevMu[4];
  TH1F *hPtTevMuJetVeto[4];
  TH1F *hPtTevMuOverPtMu[4];

  TH1F *h1PtGen1PtRecoTevMu[4];
  TH2F *h1PtGen1PtRecoVsPtGenTevMu[4];
  TH1F *hPtTevMuOverPtGen[4];


  TH1F *h_mcmu_pt, *h_mcmu_pt_hlt;
  TH1F *h_mcmu_eta, *h_mcmu_eta_hlt, *h_mcmu_eta_l1;

  //    Root output file
  edm::Service<TFileService> fs;


  edm::Handle<reco::MuonCollection> muonCollection;
  edm::Handle<reco::IsoDepositMap> tkMapH;
  edm::Handle<reco::IsoDepositMap> ecalMapH;
  edm::Handle<reco::IsoDepositMap> hcalMapH;

  const reco::TrackToTrackMap * tevMap_default;
  const reco::TrackToTrackMap * tevMap_1stHit;
  const reco::TrackToTrackMap * tevMap_picky;

  // # of trigger paths in HLT configuration
  unsigned N_triggers;

  // # of reconstructed muons per event
  unsigned N_muons;
  // Total # of reconstructed muons in job
  unsigned N_muons_tot;

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

  // get isolation
  void getIsolation(const edm::Event & iEvent);

  // get TeV muons
  void getTeVMuons(const edm::Event & iEvent);

  // get muons, update MET
  void getMuons(const edm::Event & iEvent);

  // do muon analysis
  void doMuons();
  // do TeV-muon analysis
  void doTeVanalysis(reco::MuonRef mu);

  // do isolation
  void doIsolation(reco::MuonRef mu, double massT);

  // do MC matching
  void doMCmatching();
};
