#ifndef _HadronicVWAnalyzer_h_
#define _HadronicVWAnalyzer_h_

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

class HadronicVWAnalyzer  {
public:
  HadronicVWAnalyzer();                         // constructor; initialize the list to be empty
  HadronicVWAnalyzer(const edm::ParameterSet & cfg, WPrimeUtil * wprimeUtil);
  ~HadronicVWAnalyzer();

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
  void PrintElectron(const pat::Electron& elec, int parent=0) const;
  void PrintElectron(const heep::Ele& elec, int parent=0) const;
  void PrintMuon(const TeVMuon& mu, int parent=0) const;
  
//methods for utilities
  bool SameTrigger(const std::string & A, const std::string & B) const;
  int CountZCands(ZCandV & Zs) const;
  void CalcVVariables();
  void CalcWVariables();
  void CalcWElecVariables();
  void CalcWMuonVariables();
  void CalcWVVariables();
  void CalcEventVariables();

  int   Calc_EvtType() const;
  float CalcLeadPt(int type=0) const;
  float Calc_Q() const;
  float Calc_Ht() const;
  float CalcTriLepMass() const;
  float Calc_GenWVInvMass() const;
  bool inEE(const TeVMuon& mu) const;

  bool EMuOverlap(const pat::Electron & e,
                  const MuonV & ms) const;
  
  const heep::Ele & FindElectron(const reco::Candidate & p) const;
  const TeVMuon & FindMuon(const reco::Candidate & p) const ;

//  bool Match(const heep::Ele & p1, const reco::Candidate & p2) const;
//  bool Match(const TeVMuon & p1, const reco::Candidate & p2) const;
  
  float WLepPt() const;
  float VLepPt(int idx) const;
  
  float ElecPU(const heep::Ele & e) const;
  float MuonPU(const TeVMuon & m) const;

//methods for modifiers
  void SetCandEvtFile(const std::string& s);

//methods for the cuts
  bool PassNoCut();
  bool PassTriggersCut();
  bool PassMinNLeptonsCut();
  bool PassMaxNLeptonsCut();
  bool PassMinNJetsCut();
  bool PassValidWCut();
  bool PassValidWElecCut();
  bool PassValidWMuonCut();
  bool PassValidVCut();
  bool PassValidWandVCut();
  bool PassValidWVCandCut();
  bool PassWptCut();
  bool PassVptCut();
  bool PassHtCut();
  bool PassMETCut();

  bool PassVMassCut();

  bool PassWtransMassCut();

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
  std::string looseJetType_;

  edm::InputTag electronsLabel_;
  edm::InputTag muonsLabel_;
  edm::InputTag pfCandsLabel_;
  edm::InputTag metLabel_;
  edm::InputTag hltEventLabel_;
  edm::InputTag pileupLabel_;
  edm::InputTag vertexLabel_;
  edm::InputTag jetsLabel_;
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
  float WVMass_;
  float Vpt_;
  float Wpt_;
  float weight_;
  float Ht_;
  float TriLepMass_;
  float Q_;
  uint evtType_;
  float LeadPt_;
  float LeadElecPt_;
  float LeadMuonPt_;
  bool TT, TF;

// +++++++++++++++++++General Cut values
  uint minNLeptons_;
  uint maxNLeptons_;
  uint minNJets_;
  float minMET_;

// +++++++++++++++++++Ht Cuts
  float minHt_;

// +++++++++++++++++++W Cuts
  float minWtransMass_;
  float minWpt_;

// +++++++++++++++++++V Cuts
  float minVpt_;
  float minVmass_;
  float maxVmass_;

  //Handles
  PatElectronVH patElectronsH_;
  PatMuonVH patMuonsH_;
  JetVH patJetsH_;
  METVH metH_;
  PFCandidateVH pfCandidatesH_;

//////Chosen Candidates
  ElectronV electrons_, looseElectrons_, tightElectrons_;
  MuonV muons_, looseMuons_, tightMuons_;
  JetV  looseJets_;
  pat::MET met_;
  WCandidate vCand_;
  WCandidate wCand_;
  WVCandidate wvCand_;
  pat::TriggerEvent triggerEvent_; 
  std::vector<reco::Vertex>  vertices_;
  std::vector< PileupSummaryInfo > PupInfo_; 
  wprime::EffV results_;

// +++++++++++++++++++ Histogram Definitions
  TH1F * hNumEvts;

  std::vector<TH1F*> hWVMass;
  std::vector<TH1F*> hWV3e0muMass, hWV2e1muMass, hWV1e2muMass, hWV0e3muMass;

  std::vector<TH1F*> hQ;
  std::vector<TH1F*> hWVTransMass;
  std::vector<TH1F*> hWVpt;
  std::vector<TH1F*> hHt;
  std::vector<TH1F*> hTriLepMass;
  std::vector<TH1F*> hEvtType, hEvtTypeP, hEvtTypeM;
  std::vector<TH1F*> hLeadPt, hLeadPtVee, hLeadPtVmm;
  std::vector<TH1F*> hLeadElecPt, hLeadMuonPt;

  std::vector<TH1F*> hVMass, hVeeMass, hVmmMass ;
  std::vector<TH1F*> hV3e0muMass, hV2e1muMass, hV1e2muMass, hV0e3muMass;
  std::vector<TH1F*> hVeeMassTT, hVeeMassTF, hVmmMassTT, hVmmMassTF;
  std::vector<TH1F*> hVpt,hVeept,hVmmpt;

  std::vector<TH1F*> hMET, hMETee, hMETmm;
  std::vector<TH1F*> hMET3e0mu, hMET2e1mu, hMET1e2mu, hMET0e3mu;

  std::vector<TH1F*> hWTransMass, hWenuTransMass, hWmnuTransMass;
  std::vector<TH1F*> hW3e0muTransMass, hW2e1muTransMass, hW1e2muTransMass, hW0e3muTransMass;
  std::vector<TH1F*> hWpt,hWptVee,hWptVmm;
  std::vector<TH1F*> hWQ, hWenuQ, hWmnuQ;
  std::vector<TH1F*> hW3e0muQ, hW2e1muQ, hW1e2muQ, hW0e3muQ;

  std::vector<TH1F*> hNLElec;
  std::vector<TH1F*> hNLMuon;
  std::vector<TH1F*> hNLLeps, hNLLepsVee, hNLLepsVmm;

  std::vector<TH1F*> hNTElec;
  std::vector<TH1F*> hNTMuon;
  std::vector<TH1F*> hNTLeps;

  std::vector<TH1F*> hNJets,hNJetsVee,hNJetsVmm;
  std::vector<TH1F*> hNVtxs, hNVtxsVee, hNVtxsVmm;

  std::vector<TH1F*> hWenuCombRelIso, hWmnuCombRelIso;

  TTree* tWVCand;
//Cuts
  
  int NCuts_;
  std::vector<std::string> Cuts_;
  typedef bool (HadronicVWAnalyzer::*CutFnPtr)(); 
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

  PSet jSelectorPset_;
  JetSelector looseJet_;
  pat::strbitset jetResult_;

};

#endif//#define _HadronicVWAnalyzer_h_
