#ifndef _HadronicVZAnalyzer_h_
#define _HadronicVZAnalyzer_h_

#include "UserCode/CMGWPrimeGroup/interface/AnalyzerBase.h"
#include "TH2F.h"
#include "TTree.h"

/// The class HadronicVZAnalyzer will analyze X -> VZ -> (heavy) jet + dilepton.
/// For now, we do only simple cuts and simple histograms - fancy things come later.
/// We also have only the dimuon channel for the time being. 

class HadronicVZAnalyzer : public AnalyzerBase {

  edm::InputTag genLabel_;
  edm::Handle<reco::GenParticleCollection> genParticles;


public:
  HadronicVZAnalyzer();                         // constructor; initialize the list to be empty
  HadronicVZAnalyzer(const edm::ParameterSet & cfg, int fileToRun);
  ~HadronicVZAnalyzer();

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

  METVH metH_;
  XWLeptonic wzCand_;


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


//////Chosen Vector Boson Candidates 
  ZCandidate zCand_, vCand_;
  VZCandidate hadVZ_;

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
  std::vector<TH1F*> hVZMass, hVZeeMass, hVZmmMass;
  std::vector<TH1F*> hVZpt, hQ;

  std::vector<TH1F*> hZMass, hZeeMass, hZmmMass;
  std::vector<TH1F*> hZpt;
  std::vector<TH1F*> hEvtType;
  std::vector<TH1F*> hDRee, hDRmm;

  std::vector<TH1F*> hMET, hMETee, hMETmm;


  std::vector<TH1F*> heeVMass;
  std::vector<TH1F*> hmmVMass;
  std::vector<TH1F*> heeVpt;
  std::vector<TH1F*> hmmVpt;

  std::vector<TH1F*> hNVtxs;
  std::vector<TH1F*> hNLLeps;
  std::vector<TH1F*> hNLJets;

  TTree* tVZCand;
  ////////////////////

  // http://www.parashift.com/c++-faq-lite/pointers-to-members.html#faq-33.5
  // Good manners!
  
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



#endif//#define _HadronicVZAnalyzer_h_
