#ifndef _HadronicVZAnalyzer_h_
#define _HadronicVZAnalyzer_h_

#include <vector>
#include <string>
#include <TH2F.h>
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "UserCode/CMGWPrimeGroup/interface/WPrimeUtil.h"
#include "UserCode/CMGWPrimeGroup/interface/BosonFinder.h"
#include "UserCode/CMGWPrimeGroup/interface/TeVMuon.h"

/// The class HadronicVZAnalyzer will analyze X -> VZ -> (heavy) jet + dilepton.
/// For now, we do only simple cuts and simple histograms - fancy things come later.
/// We also have only the dimuon channel for the time being. 

class HadronicVZAnalyzer  {

public:
  HadronicVZAnalyzer();                         // constructor; initialize the list to be empty
  HadronicVZAnalyzer(const edm::ParameterSet & cfg, WPrimeUtil * wprimeUtil);
  ~HadronicVZAnalyzer();

  // +++++++++++++++++++Event Characteristics
  uint numZs;
  float LeadPt;
  
  // +++++++++++++++++++General Cut values
  uint maxNumZs;
  uint minNLeptons;
  uint minNJets;
  uint maxNJets;
  float maxAngleBetweenJets;
  float minLeadPt;
  
  // +++++++++++++++++++Z Cuts
  float minZpt;
  float minZmass;
  float maxZmass;

  // +++++++++++++++++++Hadronic Boson Cuts
  float minHadVpt;
  float minHadVmass;
  float maxHadVmass;

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

  //void PrintElectron(const heep::Ele* elec, int parent);
  void PrintMuon(const TeVMuon* mu, int parent);
  void PrintJet(const pat::Jet* jet, int parent);
    
  bool SameTrigger(std::string & A, std::string & B);

//methods for utilities
//  void CheckStream(ofstream& stream, std::string s);

//methods for modifiers
  void SetCandEvtFile(std::string s);
  void SetLogFile(std::string s);
  void SetOutputFile(std::string s);

  //methods for histograms 
  void Declare_Histos(TFileDirectory& dir);
  void Fill_Histos(int index, float weight=1.);

  //clean stuff
  void ClearEvtVariables();

//methods for the cuts
  bool PassNoCut();
  bool PassTriggersCut();
  bool PassNLeptonsCut();
  bool PassNJetsCut();
  bool PassValidHadVCut();
  bool PassValidZCut();  
  bool PassValidZCutTight();
  bool PassValidHadVZCandCut();
  bool PassLeadingLeptonPtCut();
  bool PassNumberOfZsCut();
  bool PassZMassCut();
  bool PassZMassCutTight();
  bool PassZptCut();
  bool PassZptCutTight();
  bool PassHadVMassCut();
  bool PassHadVptCut();

  //  bool PassZLepPtCut();
  
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

//////////
  void reportProgress(int eventNum);
  void verbose(const char *string, ...);

  void beginFile(std::vector<wprime::InputFile>::const_iterator fi);
  void endFile(std::vector<wprime::InputFile>::const_iterator fi, ofstream & out);
  void endAnalysis(ofstream & out);

// +++++++++++++++++++useful constants
  bool debugme;//print stuff if active
  bool doPreselect_;

  std::string electronsLabel_;
  std::string muonsLabel_;
  std::string jetsLabel_;
  std::string metLabel_;
  std::string hltEventLabel_;
  std::string pileupLabel_;
  vstring triggersToUse_;

  int PDGMUON;
  int PDGELEC;
  int PDGW;
  int PDGZ;
  int PDGZPRIME;
  int PDGWPRIME;

  double PI;
  double TWOPI;
  float NOCUT;
  float weight_;
  float PU_NumInteractions3BX_;
  float PU_NumInteractions1BX_;

// +++++++++++++++++++
  std::vector<uint> nEvents;
  int eventNum;
  uint runNumber;
  uint lumiID;

  typedef std::pair< std::map<std::string, bool>, std::vector<std::string> > OptArgPair;
  std::map<std::string, int> intOptions_;
  std::map<std::string, std::string> stringOptions_;
  std::vector<std::string> filenames_;
  std::string datasetName;
  
// +++++++++++++++++++location of data files and samples info
  std::string top_level_dir;
  ofstream outCandEvt;
  ofstream outLogFile;

  WPrimeUtil * wprimeUtil_;
  
  uint muonAlgo_;

//////Chosen Candidates
  ElectronV electrons_, looseElectrons_, tightElectrons_;
  MuonV muons_, looseMuons_, tightMuons_;
  JetV jets_;
  JetV jetsDR_;
  // Not using now
  //pat::MET met;

//////Chosen Vector Boson Candidates 
  ZCandidate zCand;
  WCandidate wCand;
  WZCandidate wzCand;

  pat::TriggerEvent triggerEvent_; 
  std::vector< PileupSummaryInfo > PupInfo_; 
  std::vector<float> Num_surv_cut_;

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
  
  TeVMuon & FindMuon(reco::Candidate & p);
  bool Match(TeVMuon & p1, reco::Candidate & p2);


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

  ZCandidate zCand_;
  ZCandidate zCandTight_;


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
