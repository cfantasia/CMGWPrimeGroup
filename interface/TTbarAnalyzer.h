#ifndef _TTbarAnalyzer_h_
#define _TTbarAnalyzer_h_

#include "UserCode/CMGWPrimeGroup/interface/AnalyzerBase.h"
#include "TH2F.h"
#include "TTree.h"
#include "UserCode/CMGWPrimeGroup/interface/JetCorrectionUncertainty.h"

/// The class TTbarAnalyzer will analyze X -> VZ -> (heavy) jet + dilepton.
/// For now, we do only simple cuts and simple histograms - fancy things come later.
/// We also have only the dimuon channel for the time being. 

class TTbarAnalyzer : public AnalyzerBase {

  edm::InputTag genLabel_;
  edm::Handle<reco::GenParticleCollection> genParticles;


public:
  TTbarAnalyzer();                         // constructor; initialize the list to be empty
  TTbarAnalyzer(const edm::ParameterSet & cfg, int fileToRun);
  ~TTbarAnalyzer();

  // +++++++++++++++++++Event Characteristics

  // +++++++++++++++++++General Cut values
  uint minNLeptons_;
  uint minNJets_;
  float maxAngleBetweenJets;
  
  // +++++++++++++++++++Z Cuts
  float minZmass_;
  float maxZmass_;
  float minZpt_;

  // +++++++++++++++++++Hadronic Boson Cuts
  float minVmass_;
  float maxVmass_;
  float minVpt_;

  //Selectors
  ElectronSelector looseElectron_;
  ElectronSelector tightElectron_;

  MuonSelector looseMuon_;
  MuonSelector tightMuon_;

  JetSelector looseJet_;

  // This is actually the equivalent of the analyze() method in CMSSW.
  void eventLoop(edm::EventBase const & event);

//methods for utilities


//methods for modifiers

  //methods for histograms 
  void defineHistos(const TFileDirectory& dir);
  void fillHistos(const int& index, const float& weight=1.);
  void fillJetMultiplicityHists();
  void fillPOGMuonHists();
  //void fillTightMuonHists();
  void fillGoodZHistos();
  void fillGoodHadVHistos();
  void fillValidVZHistos();
  void fillGoodVZHistos();
  void fillJetMergingHistos();


  //clean stuff
  void clearEvtVariables();

  //printers
  void printEventDetails() const;

  //methods for the cuts
  bool passValidVZCandCut();

///My calculated qualities//////////////////
  int runNumber_;
  int lumiNumber_;
  int evtNumber_;
  float VZMass_;
  float Zpt_;
  float ZMass_;
  float Vpt_;
  float VMass_;
  float Q_;
  float DeltaRll_, DeltaRVZ_;
  uint evtType_;


//POG hists

  TH1F* h_dptpt2; 
  TH2F* h_dptpt_vs_pt;
  TH2F* h_dptpt2_vs_pt;
  TH2F* h_dptpt_vs_invpt;
  TH2F* h_dptpt2_vs_invpt;

  TH2F* h_dptpt_vs_genpt;
  TH2F* h_dptpt2_vs_genpt;
  TH2F* h_dptpt_vs_invgenpt;
  TH2F* h_dptpt2_vs_invgenpt;



  std::vector<reco::GenParticle> genMuons;
  std::vector<reco::GenParticle> genElectrons;
  std::vector<reco::GenJet> ak7GenJet;

  //Handles
  edm::Handle<std::vector<reco::Vertex> > verticesH_;
  JetVH patJetsH_;


//////Chosen Vector Boson Candidates 

  ElectronV allElectrons_, looseElectrons_, tightElectrons_;
  MuonV allMuons_, looseMuons_, tightMuons_;
  JetV looseJets_, looseBJets_, tightBJets_;
  pat::MET met_;

  pat::Jet bCand1_;
  pat::Jet bCand2_;
  WCandidate wCand_;
  XWLeptonic tCand_;
  ZCandidate zCand_, vCand_;
  VZCandidate hadVZ_;
  VZCandidate hadTop_;
  pat::Jet wJet_;

  double gravMass_;

// +++++++++++++++++++ Histogram Definitions - loose
  TH1F* h_genWMass;
  TH1F* h_genZMass;

  TH2F* h_ptJetCut_nJets;
  TH2F* h_VWMass_nJets;
  TH2F* h_VZMass_nJets;

  TH1F* h_HadVWMass;
  TH1F* h_HadVWgenMass;
  TH1F* h_HadVZgenMass;


  TH1F* h_MET_AllCuts;
  TH1F* h_WMass;

  TH1F* h_HadVZ_res;

  TH1F* h_HadVZMass;
  TH1F* h_Zelec1_pt;	
  TH1F* h_Zelec1_eta;
  TH1F* h_Zelec1_phi;
  TH1F* h_Zelec2_pt;
  TH1F* h_Zelec2_eta;
  TH1F* h_Zelec2_phi;
  TH1F* h_Zmuon1_pt;	
  TH1F* h_Zmuon1_eta;
  TH1F* h_Zmuon1_phi;
  TH1F* h_Zmuon2_pt;
  TH1F* h_Zmuon2_eta;
  TH1F* h_Zmuon2_phi;
  TH1F* h_jet1_pt;
  TH1F* h_jet1_eta;
  TH1F* h_jet1_phi;
  TH1F* h_jet2_pt;
  TH1F* h_jet2_eta;
  TH1F* h_jet2_phi;
  TH1F* h_jet1_mass;
  TH1F* h_jet2_mass;
  // TH1F* h_jet1_Vwindow;
  TH1F* h_deltaR_elec1elec2;
  TH1F* h_deltaR_HadVelec1;
  TH1F* h_deltaR_HadVelec2;
  TH1F* h_deltaR_muon1muon2;
  TH1F* h_deltaR_HadVmuon1;
  TH1F* h_deltaR_HadVmuon2;
  TH1F* h_deltaR_jet1muon1;
  TH1F* h_deltaR_jet1muon2;
  TH1F* h_deltaR_jet2muon1;
  TH1F* h_deltaR_jet2muon2;	
  TH1F* h_HadVZpt;
  TH1F* h_HadVZeta;
  TH1F* h_HadVZphi;		
  TH1F* h_muons_pt;
  TH1F* h_muons_eta;
  TH1F* h_muons_phi;
  TH1F* h_jets_pt;
  TH1F* h_jets_eta;
  TH1F* h_jets_phi;
  TH1F* h_jet_HadV_pt;
  TH1F* h_jet_HadV_eta;
  TH1F* h_jet_HadV_phi;
  TH1F* h_jet_mult;
  TH1F* h_jet_mult_inc;
  TH1F* h_jet1jet2_mass;

  //Jet merging histos
  TH2F* h_jet1mass_jet2mass;
  TH1F* h_HadVZmass_Cory;
  TH1F* h_jet1jet2_mass_Restricted;
  TH1F* h_HadVZmass_Flavia;
  TH1F* h_HadV_mass_Cory;
  TH1F* h_HadV_mass_Flavia;
  TH2F* h_m1_vs_m12;
  TH1F* h_bestmass;
  TH1F* h_deltaR_jet1jet2;
  TH1F* h_deltaR_jet1jet2_R1_cut40;
  TH1F* h_deltaR_jet1jet2_R1_cut50;
  TH1F* h_deltaR_jet1jet2_R1_cut60;
  TH1F* h_deltaR_jet1jet2_R1_cut65;
  TH1F* h_deltaR_jet1jet2_R1_cut70;
  TH1F* h_deltaR_jet1jet2_R1_cut80;
  TH1F* h_deltaR_jet1jet2_R1_cut90;
  TH1F* h_deltaR_jet1jet2_R2_cut40;
  TH1F* h_deltaR_jet1jet2_R2_cut50;
  TH1F* h_deltaR_jet1jet2_R2_cut60;
  TH1F* h_deltaR_jet1jet2_R2_cut65;
  TH1F* h_deltaR_jet1jet2_R2_cut70;
  TH1F* h_deltaR_jet1jet2_R2_cut80;
  TH1F* h_deltaR_jet1jet2_R2_cut90;
  TH1F* h_deltaR_jet1Z_R1;
  TH1F* h_deltaR_jet2Z_R1;
  TH1F* h_deltaR_jet3Z_R1;
  TH1F* h_deltaR_jet1Z_R2;
  TH1F* h_deltaR_jet2Z_R2;
  TH1F* h_deltaR_jet3Z_R2;

  /////////////////////
  std::vector<TH1F*> hTBMass, hTBenuMass, hTBmnuMass;
  std::vector<TH1F*> hTMass, hTenuMass, hTmnuMass, hQ;
  std::vector<TH1F*> hT2Mass, hT2pt;
  std::vector<TH1F*> hEvtType;
  std::vector<TH1F*> hTpt, hTenupt, hTmnupt;

  std::vector<TH1F*> hMET, hMETSig;
  std::vector<TH1F*> hWTransMass, hWenuTransMass, hWmnuTransMass;
  std::vector<TH1F*> hWpt;
  std::vector<TH1F*> hW2TransMass, hW2Mass, hW2pt;
  std::vector<TH1F*> hB1Mass, hB1Disc, hB1pt;
  std::vector<TH1F*> hB2Mass, hB2Disc, hB2pt;
  std::vector<TH1F*> hMbl;
  std::vector<TH1F*> hNLElec, hNLMuon, hNLLeps;
  std::vector<TH1F*> hNLJets, hNLBJets, hNTBJets;
  std::vector<TH1F*> hWeight, hNVtxs;
  TH2F* hW2PtVsW2Mass;

  float Mbl_;
  ////////////////////

  // http://www.parashift.com/c++-faq-lite/pointers-to-members.html#faq-33.5
  // Good manners!
  
};




#endif//#define _TTbarAnalyzer_h_
