#ifndef _WZAnalyzer_h_
#define _WZAnalyzer_h_

#include <vector>
#include <string>

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
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "UserCode/CMGWPrimeGroup/interface/util.h"
#include "UserCode/CMGWPrimeGroup/interface/WPrimeUtil.h"
#include "UserCode/CMGWPrimeGroup/interface/BosonFinder.h"
#include "UserCode/CMGWPrimeGroup/interface/TeVMuon.h"
#include "TH3F.h"
#include "TTree.h"

class WZAnalyzer  {
public:
  WZAnalyzer();                         // constructor; initialize the list to be empty
  WZAnalyzer(const edm::ParameterSet & cfg, WPrimeUtil * wprimeUtil);
  ~WZAnalyzer();

  //methods for stuff to be once per job
  void FillCutFns();
  void endAnalysis(ofstream & out);


  //methods for stuff to be done for each sample
  void beginFile(std::vector<wprime::InputFile>::const_iterator fi);
  void endFile(std::vector<wprime::InputFile>::const_iterator fi, ofstream & out);
  void ResetCounters();
  void Declare_Histos(TFileDirectory& dir);
  void DeclareHistoSet(std::string n, std::string t, std::string xtitle,
                       int nbins, float min, float max, std::string units,
                       std::vector<TH1F*>& h, TFileDirectory& d);
  void DeclareHisto(std::string n, std::string t, std::string xtitle,
                    int nbins, float min, float max,
                    TH1F* h, TFileDirectory& d);
  void tabulateSummary();
  void printSummary(const std::string& dir, ofstream & out) const;

  //methods for stuff to be done for each event
  void eventLoop(edm::EventBase const & event);
  bool PassCuts(const float& weight=1.);
  void ClearEvtVariables();
  void Fill_Histos(int index, float weight=1.);
  void Tabulate_Me(int& cut_index,const float& weight);

  //methods for printers
  void PrintEventFull(edm::EventBase const & event) const;
  void PrintPassingEvent(edm::EventBase const & event);
  void PrintDebugEvent() const;
  void PrintEventToFile(edm::EventBase const & event);
  void PrintEvent(edm::EventBase const & event) const;
  void PrintEventDetails() const;
  void PrintEventLeptons() const;
  void PrintTrigger() const;
  void PrintLeptons() const;
  void PrintElectron(const heep::Ele& elec, int parent=0) const;
  void PrintMuon(const TeVMuon& mu, int parent=0) const;
  
//methods for utilities
  bool SameTrigger(const std::string & A, const std::string & B) const;
  int CountZCands(ZCandV & Zs) const;
  void CalcZVariables();
  void CalcWVariables();
  void CalcWElecVariables();
  void CalcWMuonVariables();
  void CalcWZVariables();
  void CalcEventVariables();

  int   Calc_EvtType() const;
  float CalcLeadPt(int type=0) const;
  float Calc_Q() const;
  float Calc_Ht() const;
  float CalcTriLepMass() const;
  float Calc_GenWZInvMass() const;
  bool inEE(const TeVMuon& mu) const;

  bool EMuOverlap(const pat::Electron & e,
                  const std::vector<pat::Muon > & ms) const;
  
  const heep::Ele & FindElectron(const reco::Candidate & p) const;
  const TeVMuon & FindMuon(const reco::Candidate & p) const ;

  bool Match(const heep::Ele & p1, const reco::Candidate & p2) const;
  bool Match(const TeVMuon & p1, const reco::Candidate & p2) const;
  
  float WLepPt() const;
  float ZLepPt(int idx) const;
  
  float ElecPU(const heep::Ele & e) const;
  float MuonPU(const TeVMuon & m) const;

//methods for modifiers
  void SetCandEvtFile(const std::string& s);

//methods for the cuts
  bool PassNoCut();
  bool PassTriggersCut();
  bool PassMinNLeptonsCut();
  bool PassMaxNLeptonsCut();
  bool PassValidWCut();
  bool PassValidWElecCut();
  bool PassValidWMuonCut();
  bool PassValidZCut();
  bool PassValidWandZCut();
  bool PassValidWZCandCut();
  bool PassLeadingLeptonPtCut();
  bool PassNumberOfZsCut();
  bool PassWptCut();
  bool PassZptCut();
  bool PassHtCut();
  bool PassMETCut();

  bool PassZMassCut();
  bool PassZLepPtCut();
  bool PassZLepTriggerMatchCut();
  bool PassZeePtCut();
  bool PassZmumuPtCut();

  bool PassWtransMassCut();
  bool PassWLepPtCut();
  bool PassWLepIsoCut();

  bool PassWLepTightCut();
  bool PassWFlavorElecCut();
  bool PassWFlavorMuonCut();
  bool PassFakeEvtCut();
  bool PassFakeLeptonTagCut();
  bool PassFakeLeptonProbeLooseCut();
  bool PassFakeLeptonProbeTightCut();

  bool PassTriggerMatch(const heep::Ele& e1, const heep::Ele& e2) const;
  bool PassTriggerMatch(const TeVMuon& m1, const TeVMuon& m2) const;
  bool PassTriggerMatch(const pat::Electron& p, const float cut, const vstring& triggers) const;
  bool PassTriggerMatch(const TeVMuon& p, const float cut, const vstring& triggers) const;
  bool PassTriggerEmulation(const heep::Ele& elec, const float minPt=0.) const;


//////////////////
/////Variables////
//////////////////

  bool debugme;//print stuff if active
  bool doPreselect_;

  std::string looseElectronType_, tightElectronType_;
  std::string looseMuonType_, tightMuonType_;

  std::string electronsLabel_;
  std::string muonsLabel_;
  std::string metLabel_;
  std::string hltEventLabel_;
  std::string pileupLabel_;
  vstring triggersToUse_;

// +++++++++++++++++++location of data files and samples info
  WPrimeUtil * wprimeUtil_;
  ofstream outCandEvt_;

  uint muonAlgo_;
  bool useAdjustedMET_;
  double rhoFastJet_;
  std::vector<double> effectiveElecArea_;
  std::vector<double> effectiveMuonArea_;

///My calculated qualities//////////////////
  float WZMass_;
  float Zpt_;
  float Wpt_;
  float weight_;
  float Ht_;
  float TriLepMass_;
  float Q_;
  uint evtType_;
  uint numZs_;
  float LeadPt_;
  float LeadElecPt_;
  float LeadMuonPt_;
  bool TT, TF;
  float PU_NumInteractions_;

// +++++++++++++++++++General Cut values
  uint maxNumZs_;
  uint minNLeptons_;
  uint maxNLeptons_;
  float minLeadPt_;
  float minMET_;

// +++++++++++++++++++Ht Cuts
  float minHt_;

// +++++++++++++++++++W Cuts
  float minWtransMass_;
  float minWpt_;

  int   cutWenuWPRelIsoMask_;
  std::string cutElecWPTightType_;
  float minDeltaR_;

  float minWlepPt_;
// +++++++++++++++++++Z Cuts
  float maxZMassDiff_;

  float minZpt_;
  float minZmass_;
  float maxZmass_;
  float minZeePt1_;
  float minZeePt2_;
  float minZmmPt1_;
  float minZmmPt2_;

//////Chosen Candidates
  ElectronV electrons_, looseElectrons_, tightElectrons_;
  MuonV muons_, looseMuons_, tightMuons_;
  JetV  jets_;
  pat::MET met_;
  ZCandidate zCand_;
  WCandidate wCand_;
  WZCandidate wzCand_;
  pat::TriggerEvent triggerEvent_; 
  std::vector< PileupSummaryInfo > PupInfo_; 
  std::vector<wprime::FilterEff> results_;

// +++++++++++++++++++ Histogram Definitions
  TH1F * hEffRel;
  TH1F * hEffAbs;
  TH1F * hNumEvts;

  std::vector<TH1F*> hEvtType;
  std::vector<TH1F*> hEvtTypeP;
  std::vector<TH1F*> hEvtTypeM;

  std::vector<TH1F*> hWZMass     ;
  std::vector<TH1F*> hWZ3e0muMass;
  std::vector<TH1F*> hWZ2e1muMass;
  std::vector<TH1F*> hWZ1e2muMass;
  std::vector<TH1F*> hWZ0e3muMass;

  std::vector<TH1F*> hWZTransMass;
  std::vector<TH1F*> hHt;
  std::vector<TH1F*> hWpt;
  std::vector<TH1F*> hZpt;

  std::vector<TH1F*> hMET;
  std::vector<TH1F*> hMETee;
  std::vector<TH1F*> hMETmumu;
  std::vector<TH1F*> hMET3e0mu;
  std::vector<TH1F*> hMET2e1mu;
  std::vector<TH1F*> hMET1e2mu;
  std::vector<TH1F*> hMET0e3mu;

  std::vector<TH1F*> hZMass     ;
  std::vector<TH1F*> hZeeMass   ;
  std::vector<TH1F*> hZmumuMass ;
  std::vector<TH1F*> hZ3e0muMass;
  std::vector<TH1F*> hZ2e1muMass;
  std::vector<TH1F*> hZ1e2muMass;
  std::vector<TH1F*> hZ0e3muMass;
  std::vector<TH1F*> hZeeMassTT;
  std::vector<TH1F*> hZeeMassTF;
  std::vector<TH1F*> hZmumuMassTT ;
  std::vector<TH1F*> hZmumuMassTF ;

  std::vector<TH1F*> hWTransMass     ;
  std::vector<TH1F*> hWenuTransMass  ;
  std::vector<TH1F*> hWmunuTransMass ;
  std::vector<TH1F*> hW3e0muTransMass;
  std::vector<TH1F*> hW2e1muTransMass;
  std::vector<TH1F*> hW1e2muTransMass;
  std::vector<TH1F*> hW0e3muTransMass;

  std::vector<TH1F*> hWQ     ;
  std::vector<TH1F*> hWenuQ  ;
  std::vector<TH1F*> hWmunuQ ;
  std::vector<TH1F*> hW3e0muQ;
  std::vector<TH1F*> hW2e1muQ;
  std::vector<TH1F*> hW1e2muQ;
  std::vector<TH1F*> hW0e3muQ;

  std::vector<TH1F*> hQ;

  std::vector<TH1F*> hNLElec;
  std::vector<TH1F*> hNLMuon;
  std::vector<TH1F*> hNLLeps;

  std::vector<TH1F*> hNTElec;
  std::vector<TH1F*> hNTMuon;
  std::vector<TH1F*> hNTLeps;

  std::vector<TH1F*> hNJets;
  std::vector<TH1F*> hNVtxs;

  std::vector<TH1F*> hLeadPt;
  std::vector<TH1F*> hLeadElecPt;
  std::vector<TH1F*> hLeadMuonPt;

  std::vector<TH1F*> hWenuCombRelIso;
  std::vector<TH1F*> hWmunuCombRelIso;

  std::vector<TH1F*> hElecPt;
  std::vector<TH1F*> hElecEt;
  std::vector<TH1F*> hElecdEta;
  std::vector<TH1F*> hElecdPhi;
  std::vector<TH1F*> hElecSigmann;
  std::vector<TH1F*> hElecEP;
  std::vector<TH1F*> hElecHE;
  std::vector<TH1F*> hElecTrkRelIso;
  std::vector<TH1F*> hElecECalRelIso;
  std::vector<TH1F*> hElecHCalRelIso;

  std::vector<TH1F*> hMuonPt;
  std::vector<TH1F*> hMuonDxy;
  std::vector<TH1F*> hMuonNormChi2;
  std::vector<TH1F*> hMuonNPix;
  std::vector<TH1F*> hMuonNTrk;
  std::vector<TH1F*> hMuonRelIso;
  std::vector<TH1F*> hMuonStation;
  std::vector<TH1F*> hMuonSip;
  std::vector<TH1F*> hMuonTightCombIso;

  std::vector<TH1F*> hTriLepMass;

  TTree* tWZCand;
//Cuts
  
  int NCuts_;
  std::vector<std::string> Cuts_;
  typedef bool (WZAnalyzer::*CutFnPtr)(); 
#ifndef __CINT__
  std::map<std::string,CutFnPtr> mFnPtrs_;
  std::vector<CutFnPtr> CutFns_;
#endif

  PSet eSelectorPset_;
  ElectronSelector looseElectron_;
  ElectronSelector tightElectron_;
  pat::strbitset electronResult_;

  PSet mSelectorPset_;
  MuonSelector looseMuon_;
  MuonSelector tightMuon_;
  pat::strbitset muonResult_;

};

#endif//#define _WZAnalyzer_h_
