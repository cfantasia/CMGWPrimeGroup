#ifndef _HadronicVZAnalyzer_h_
#define _HadronicVZAnalyzer_h_

#include "UserCode/CMGWPrimeGroup/interface/AnalyzerBase.h"
#include "TH2F.h"

/// The class HadronicVZAnalyzer will analyze X -> VZ -> (heavy) jet + dilepton.
/// For now, we do only simple cuts and simple histograms - fancy things come later.
/// We also have only the dimuon channel for the time being. 

class HadronicVZAnalyzer : public AnalyzerBase {

public:
  HadronicVZAnalyzer();                         // constructor; initialize the list to be empty
  HadronicVZAnalyzer(const edm::ParameterSet & cfg, WPrimeUtil * wprimeUtil);
  ~HadronicVZAnalyzer();

  // +++++++++++++++++++Event Characteristics

  // +++++++++++++++++++General Cut values
  uint maxNJets;
  float maxAngleBetweenJets;
  
  // +++++++++++++++++++Z Cuts

  // +++++++++++++++++++Hadronic Boson Cuts

// +++++++++++++++++++Muon General Cuts
  float maxMuonEta;
  float minMuonLoosePt;
  float minMuonTightPt;
//VBTF Recommended Cuts
  float maxMuonDxy;
  float maxMuonNormChi2;
  int minMuonNPixHit;
  int minMuonNTrkHit;
  int minMuonStations;
  int minMuonHitsUsed;

// +++++++++++++++++++Jet General Cuts
  float minJetPt;
  float maxJetEta;
  //Jet ID Cuts
  float maxJetNHF;
  float maxJetNEF;
  float minJetCHF;
  float maxJetCEF;
  size_t minJetnumConst;
  size_t minJetcMult;

  // This is actually the equivalent of the analyze() method in CMSSW.
  void eventLoop(edm::EventBase const & event);

//methods for utilities

//methods for modifiers

  //methods for histograms 
  void Declare_Histos(const TFileDirectory& dir);
  void Fill_Histos(const int& index, const float& weight=1.);
  void FillJetMultiplicityHists();
  void FillLooseMuonHists();
  void FillTightMuonHists();
  void FillGoodZHistos();
  void FillGoodZTightHistos();
  void FillGoodHadVHistos();
  void FillValidVZHistos();
  void FillGoodVZHistos();
  void FillJetMergingHistos();


  //clean stuff
  void ClearEvtVariables();

//methods for the cuts
  bool PassMinNTightLeptonsCut();
  bool PassValidZCutTight();
  bool PassValidVZCandCut();
  bool PassValidVZTightCandCut();

  bool PassZMassCutTight();
  bool PassZptCutTight();

  bool PassMuonCut(const TeVMuon* mu);
  bool PassMuonLooseCut(const TeVMuon* mu);
  bool PassMuonTightCut(const TeVMuon* mu);
  bool PassMuonLoosePtCut(const TeVMuon* mu);
  bool PassMuonTightPtCut(const TeVMuon* mu);
  bool PassMuonGlobalCut(const TeVMuon* mu);
  bool PassMuonDxyCut(const TeVMuon* mu);
  bool PassMuonNpixhitCut(const TeVMuon* mu);
  bool PassMuonNtrkhitCut(const TeVMuon* mu);
  bool PassMuonNormChi2Cut(const TeVMuon* mu);
  bool PassMuonHitsUsedCut(const TeVMuon* mu);
  bool PassMuonStationsCut(const TeVMuon* mu);
  bool PassMuonEtaCut(const TeVMuon* mu);
  bool PassMuonCombRelIsoCut(const TeVMuon* mu);

  bool PassJetCut(const pat::Jet* jet);
  bool PassJetPtCut(const pat::Jet* jet);
  bool PassJetEtaCut(const pat::Jet* jet);
  bool PassJetNHFCut(const pat::Jet* jet);
  bool PassJetNEFCut(const pat::Jet* jet);
  bool PassJetNConstCut(const pat::Jet* jet);
  bool PassJetCHFCut(const pat::Jet* jet);
  bool PassJetCMultCut(const pat::Jet* jet);
  bool PassJetCEFCut(const pat::Jet* jet);
  bool PassJetIDCut(const pat::Jet* jet); 


// +++++++++++++++++++useful constants

  float weight_;
  
// +++++++++++++++++++location of data files and samples info

//////Chosen Vector Boson Candidates 
  ZCandidate zCandTight_;
  VZCandidate hadVZ_, hadVZTight_;

// +++++++++++++++++++ Histogram Definitions - loose
  TH1F* h_HadVZMass;
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
  TH1F* h_Zmuon1_VZCut_pt;
  TH1F* h_Zmuon1_VZCut_eta;
  TH1F* h_Zmuon1_VZCut_phi;
  TH1F* h_Zmuon2_VZCut_pt;
  TH1F* h_Zmuon2_VZCut_eta;
  TH1F* h_Zmuon2_VZCut_phi;
  TH1F* h_jet_VZCut_pt;
  TH1F* h_jet_VZCut_eta;
  TH1F* h_jet_VZCut_phi;
  TH1F* h_jet_mult;
  TH1F* h_jet_mult_inc;
  TH1F* h_jet1jet2_mass;



// +++++++++++++++++++ Histogram Definitions - tight
  TH1F* h_tight_HadVZMass;
  TH1F* h_tight_Zmuon1_pt;	
  TH1F* h_tight_Zmuon1_eta;
  TH1F* h_tight_Zmuon1_phi;
  TH1F* h_tight_Zmuon2_pt;
  TH1F* h_tight_Zmuon2_eta;
  TH1F* h_tight_Zmuon2_phi;
  /*  TH1F* h_tight_jet1_pt;
  TH1F* h_tight_jet1_eta;
  TH1F* h_tight_jet1_phi;
  TH1F* h_tight_jet2_pt;
  TH1F* h_tight_jet2_eta;
  TH1F* h_tight_jet2_phi;
  TH1F* h_tight_jet1_mass;
  TH1F* h_tight_jet2_mass;
  */
  TH1F* h_tight_deltaR_muon1muon2;
  TH1F* h_tight_deltaR_HadVmuon1;
  TH1F* h_tight_deltaR_HadVmuon2;
  /*  TH1F* h_tight_deltaR_jet1muon1;
  TH1F* h_tight_deltaR_jet1muon2;
  TH1F* h_tight_deltaR_jet2muon1;
  TH1F* h_tight_deltaR_jet2muon2;*/	
  TH1F* h_tight_HadVZpt;
  TH1F* h_tight_HadVZeta;
  TH1F* h_tight_HadVZphi;		
  TH1F* h_tight_muons_pt;
  TH1F* h_tight_muons_eta;
  TH1F* h_tight_muons_phi;
  /*  TH1F* h_tight_jets_pt;
  TH1F* h_tight_jets_eta;
  TH1F* h_tight_jets_phi;
  TH1F* h_tight_jet_HadV_pt;
  TH1F* h_tight_jet_HadV_eta;
  TH1F* h_tight_jet_HadV_phi;*/
  TH1F* h_tight_Zmuon1_VZCut_pt;
  TH1F* h_tight_Zmuon1_VZCut_eta;
  TH1F* h_tight_Zmuon1_VZCut_phi;
  TH1F* h_tight_Zmuon2_VZCut_pt;
  TH1F* h_tight_Zmuon2_VZCut_eta;
  TH1F* h_tight_Zmuon2_VZCut_phi;
  TH1F* h_tight_jet_VZCut_pt;
  TH1F* h_tight_jet_VZCut_eta;
  TH1F* h_tight_jet_VZCut_phi;
  //  TH1F* h_tight_jet_mult;
  //  TH1F* h_tight_jet_mult_inc;
  //  TH1F* h_tight_jet1jet2_mass;

  //Jet merging histos
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

  TH1F* h_deltaR_jet1Z;
  TH1F* h_deltaR_jet2Z;
  TH1F* h_deltaR_jet3Z;


  //Muon work histos
  TH1F* h_dptpt2; 
  TH2F* h_dptpt_vs_pt;
  TH2F* h_dptpt2_vs_pt;
  TH2F* h_dptpt_vs_invpt;
  TH2F* h_dptpt2_vs_invpt;
  TH2F* h_dptpt_vs_eta;
  TH2F* h_dptpt2_vs_eta;
  TH2F* h_dptpt_vs_pt_meta09;
  TH2F* h_dptpt2_vs_pt_meta09;
  TH2F* h_dptpt_vs_invpt_meta09;
  TH2F* h_dptpt2_vs_invpt_meta09;
  TH2F* h_dptpt_vs_pt_meta0912;
  TH2F* h_dptpt2_vs_pt_meta0912;
  TH2F* h_dptpt_vs_invpt_meta0912;
  TH2F* h_dptpt2_vs_invpt_meta0912;
  TH2F* h_dptpt_vs_pt_meta1225;
  TH2F* h_dptpt2_vs_pt_meta1225;
  TH2F* h_dptpt_vs_invpt_meta1225;
  TH2F* h_dptpt2_vs_invpt_meta1225;
  TH2F* h_dptpt_vs_pt_teta09;
  TH2F* h_dptpt2_vs_pt_teta09;
  TH2F* h_dptpt_vs_invpt_teta09;
  TH2F* h_dptpt2_vs_invpt_teta09;
  TH2F* h_dptpt_vs_pt_teta0915;
  TH2F* h_dptpt2_vs_pt_teta0915;
  TH2F* h_dptpt_vs_invpt_teta0915;
  TH2F* h_dptpt2_vs_invpt_teta0915;
  TH2F* h_dptpt_vs_pt_teta1524;
  TH2F* h_dptpt2_vs_pt_teta1524;
  TH2F* h_dptpt_vs_invpt_teta1524;
  TH2F* h_dptpt2_vs_invpt_teta1524;


  // http://www.parashift.com/c++-faq-lite/pointers-to-members.html#faq-33.5
  // Good manners!
  
  // MuonCutFnPtr type is: "pointer to member function of HadronicVZAnalyzer"
  // It takes a const TeVMuon* as argument, and returns a bool.
  typedef bool (HadronicVZAnalyzer::*MuonCutFnPtr)(const TeVMuon*); 
  // Vector of strings which define cuts.
  std::vector<std::string> MuonCuts_;
  int NMuonCuts_;
  // Vector of member function pointers which APPLY those cuts
  std::vector<MuonCutFnPtr> MuonCutFns_;
  // Map between strings and member function pointers
  std::map<std::string, MuonCutFnPtr> mMuonFnPtrs_;

};

struct highestMuonPt {                                                                                                                                     
  bool operator() (const TeVMuon & a, const TeVMuon & b){                                                                                  
                                                                                                                                                             
    return a.pt() > b.pt();                                                                                                                                
  }                                                                                                                                                        
};


struct highestJetPt {
  bool operator() (const pat::Jet & a, const pat::Jet & b){
    return a.pt() > b.pt();
  }
};



#endif//#define _HadronicVZAnalyzer_h_
