#ifndef _WZAnalyzer_h_
#define _WZAnalyzer_h_

#include "UserCode/CMGWPrimeGroup/interface/AnalyzerBase.h"
#include "TTree.h"
#include <boost/signals.hpp>

class WZAnalyzer : public AnalyzerBase, public boost::signals::trackable {
public:
  WZAnalyzer();                         // constructor; initialize the list to be empty
  WZAnalyzer(const edm::ParameterSet & cfg, int fileToRun);
  ~WZAnalyzer();

  //methods for stuff to be once per job
  void setupCutOrder();

  //methods for stuff to be done for each sample
  void defineHistos(const TFileDirectory& dir);

  //methods for stuff to be done for each event
  void eventLoop(edm::EventBase const & event);
  void clearEvtVariables();
  void fillHistos(const int& index, const float& weight=1.);

  //methods for printers
  void printDebugEvent() const;
  void printEventDetails() const;
  void printEventLeptons() const;
  
//methods for utilities
  int countZCands(ZCandV & Zs) const;
  void calcZVariables();
  void calcWVariables();
  void calcWZVariables();
  void calcEventVariables();

  int   calcEvtType() const;
  float calcLeadPt(int type=0) const;
  float calcQ() const;
  float calcHt() const;
  float calcTriLepMass() const;
  float calcGenWZInvMass() const;
  bool inEE(const TeVMuon& mu) const;

  float WLepPt() const;
  float ZLepPt(int idx) const;
  
  float ElecPU(const heep::Ele & e) const;
  float MuonPU(const TeVMuon & m) const;

//methods for modifiers

//methods for the cuts
  bool passLeadingLeptonPtCut() const;
  bool passNumberOfZsCut() const;
  bool passValidWCut(WCandidate& w);
  bool passValidZCut(ZCandidate& z);
  bool passValidWZCut(XWLeptonic& xw);
  bool passHtCut() const;

  bool passZLepPtCut() const;
  bool passZeePtCut() const;
  bool passZmumuPtCut() const;

  bool passWLepPtCut() const;
  bool passWLepIsoCut() const;

  bool passWLepTightCut() const;
  bool passWFlavorElecCut() const;
  bool passWFlavorMuonCut() const;
  bool passFakeEvtCut() const;
  bool passFakeLeptonProbeCut() const;

  bool passTriggerMatch(const heep::Ele& e, const float cut, const vstring& triggers) const;
  bool passTriggerMatch(const TeVMuon& p, const float cut, const vstring& triggers) const;
  bool passTriggerEmulation(const heep::Ele& elec, const float minPt=0.) const;

//////////////////
/////Variables////
//////////////////

  double rhoFastJet_;
  std::vector<double> effectiveElecArea_;
  std::vector<double> effectiveMuonArea_;

///My calculated qualities//////////////////
  int runNumber_;
  int lumiNumber_;
  int evtNumber_;
  float WZMass_;
  float Zpt_;
  float ZMass_;
  float Wpt_;
  float WTransMass_;
  float Ht_;
  float TriLepMass_;
  float Q_;
  float MET_;
  float METSig_;
  uint evtType_;
  uint numZs_; 
  uint NVtxs_;
  float LeadPt_;
  float LeadElecPt_;
  float LeadMuonPt_;
  bool TT, TF;

// +++++++++++++++++++General Cut values
  uint minNLeptons_;
  uint maxNumZs_;
  float minLeadPt_;
  float minMET_;

// +++++++++++++++++++Ht Cuts
  float minHt_;

// +++++++++++++++++++W Cuts
  float minDeltaR_;

  float minWlepPt_;
  float minWtransMass_;
  float minWpt_;

// +++++++++++++++++++Z Cuts
  float minZmass_;
  float maxZmass_;
  float maxZMassDiff_;
  float minZpt_;

  float minZeePt1_;
  float minZeePt2_;
  float minZmmPt1_;
  float minZmmPt2_;

  //Selectors
  ElectronSelector looseElectron_;
  ElectronSelector tightElectron_;

  MuonSelector looseMuon_;
  MuonSelector tightMuon_;

  //Handles
  PatElectronVH patElectronsH_;
  PatMuonVH patMuonsH_;
  METVH metH_;
  PFCandidateVH pfCandidatesH_;

//////Chosen Candidates
  ZCandidate zCand_;
  WCandidate wCand_;
  XWLeptonic wzCand_;
  NuAlgos wzAlgo_;

// +++++++++++++++++++ Histogram Definitions
  std::vector<TH1F*> hWZMass;
  std::vector<TH1F*> hWZ3e0muMass, hWZ2e1muMass, hWZ1e2muMass, hWZ0e3muMass;

  std::vector<TH1F*> hQ;
  std::vector<TH1F*> hWZTransMass;
  std::vector<TH1F*> hWZpt;
  std::vector<TH1F*> hHt;
  std::vector<TH1F*> hTriLepMass;
  std::vector<TH1F*> hEvtType, hEvtTypeP, hEvtTypeM;
  std::vector<TH1F*> hLeadPt, hLeadPtZee, hLeadPtZmm;
  std::vector<TH1F*> hLeadElecPt, hLeadMuonPt;

  std::vector<TH1F*> hZMass, hZeeMass, hZmmMass ;
  std::vector<TH1F*> hZeeMassTT, hZeeMassTF, hZmmMassTT, hZmmMassTF;
  std::vector<TH1F*> hZpt,hZeept,hZmmpt;

  std::vector<TH1F*> hMET, hMETee, hMETmm, hMETSig;

  std::vector<TH1F*> hWTransMass, hWenuTransMass, hWmnuTransMass;
  std::vector<TH1F*> hWpt;
  std::vector<TH1F*> hWQ;
  std::vector<TH1F*> hWTheta;

  std::vector<TH1F*> hNLElec, hNLMuon, hNLLeps;

  std::vector<TH1F*> hNTElec, hNTMuon, hNTLeps;

  std::vector<TH1F*> hNJets;
  std::vector<TH1F*> hNVtxs;

  std::vector<TH1F*> hWenuCombRelIso, hWmnuCombRelIso;

  TH1F* hVtxMatch;
  
  TTree* tWZCand;

};

#endif//#define _WZAnalyzer_h_
