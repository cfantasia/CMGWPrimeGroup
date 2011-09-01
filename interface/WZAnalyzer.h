#ifndef _WZAnalyzer_h_
#define _WZAnalyzer_h_

#include "UserCode/CMGWPrimeGroup/interface/AnalyzerBase.h"
#include "TTree.h"

class WZAnalyzer : public AnalyzerBase {
public:
  WZAnalyzer();                         // constructor; initialize the list to be empty
  WZAnalyzer(const edm::ParameterSet & cfg, WPrimeUtil * wprimeUtil);
  ~WZAnalyzer();

  //methods for stuff to be once per job
  void FillCutFns();

  //methods for stuff to be done for each sample
  void Declare_Histos(const TFileDirectory& dir);

  //methods for stuff to be done for each event
  void eventLoop(edm::EventBase const & event);
  bool PassCuts(const float& weight=1.);
  void ClearEvtVariables();
  void Fill_Histos(const int& index, const float& weight=1.);

  //methods for printers
  void PrintDebugEvent() const;
  void PrintEventDetails() const;
  void PrintEventLeptons() const;
  
//methods for utilities
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

  float WLepPt() const;
  float ZLepPt(int idx) const;
  
  float ElecPU(const heep::Ele & e) const;
  float MuonPU(const TeVMuon & m) const;

//methods for modifiers

//methods for the cuts
  bool PassTriggersCut() const;
  bool PassValidWElecCut() const;
  bool PassValidWMuonCut() const;
  bool PassLeadingLeptonPtCut() const;
  bool PassNumberOfZsCut() const;
  bool PassValidWZCut() const;
  bool PassHtCut() const;

  bool PassZLepPtCut() const;
  bool PassZLepTriggerMatchCut() const;
  bool PassZeePtCut() const;
  bool PassZmumuPtCut() const;

  bool PassWLepPtCut() const;
  bool PassWLepIsoCut() const;

  bool PassWLepTightCut() const;
  bool PassWFlavorElecCut() const;
  bool PassWFlavorMuonCut() const;
  bool PassFakeEvtCut() const;
  bool PassFakeLeptonTagCut() const;
  bool PassFakeLeptonProbeTightCut() const;

  bool PassTriggerMatch(const heep::Ele& e1, const heep::Ele& e2) const;
  bool PassTriggerMatch(const TeVMuon& m1, const TeVMuon& m2) const;
  bool PassTriggerMatch(const pat::Electron& p, const float cut, const vstring& triggers) const;
  bool PassTriggerMatch(const TeVMuon& p, const float cut, const vstring& triggers) const;
  bool PassTriggerEmulation(const heep::Ele& elec, const float minPt=0.) const;

//////////////////
/////Variables////
//////////////////

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

// +++++++++++++++++++General Cut values
  uint maxNumZs_;
  float minLeadPt_;

// +++++++++++++++++++Ht Cuts
  float minHt_;

// +++++++++++++++++++W Cuts
  float minDeltaR_;

  float minWlepPt_;
// +++++++++++++++++++Z Cuts
  float maxZMassDiff_;

  float minZeePt1_;
  float minZeePt2_;
  float minZmmPt1_;
  float minZmmPt2_;

  //Handles
  PatElectronVH patElectronsH_;
  PatMuonVH patMuonsH_;
  METVH metH_;
  PFCandidateVH pfCandidatesH_;

  edm::InputTag vertexLabel_;

//////Chosen Candidates
  DiBosonWLeptonic wzCand_;
  std::vector<reco::Vertex>  vertices_;

// +++++++++++++++++++ Histogram Definitions
  std::vector<TH1F*> hWZMass;
  std::vector<TH1F*> hWZ3e0muMass, hWZ2e1muMass, hWZ1e2muMass, hWZ0e3muMass;

  std::vector<TH1F*> hQ;
  std::vector<TH1F*> hWZTransMass;
  std::vector<TH1F*> hWZpt, hWZTheta;
  std::vector<TH1F*> hHt;
  std::vector<TH1F*> hTriLepMass;
  std::vector<TH1F*> hEvtType, hEvtTypeP, hEvtTypeM;
  std::vector<TH1F*> hLeadPt, hLeadPtZee, hLeadPtZmm;
  std::vector<TH1F*> hLeadElecPt, hLeadMuonPt;

  std::vector<TH1F*> hZMass, hZeeMass, hZmmMass ;
  std::vector<TH1F*> hZ3e0muMass, hZ2e1muMass, hZ1e2muMass, hZ0e3muMass;
  std::vector<TH1F*> hZeeMassTT, hZeeMassTF, hZmmMassTT, hZmmMassTF;
  std::vector<TH1F*> hZpt,hZeept,hZmmpt;

  std::vector<TH1F*> hMET, hMETee, hMETmm;
  std::vector<TH1F*> hMET3e0mu, hMET2e1mu, hMET1e2mu, hMET0e3mu;

  std::vector<TH1F*> hWTransMass, hWenuTransMass, hWmnuTransMass;
  std::vector<TH1F*> hW3e0muTransMass, hW2e1muTransMass, hW1e2muTransMass, hW0e3muTransMass;
  std::vector<TH1F*> hWpt,hWptZee,hWptZmm;
  std::vector<TH1F*> hWQ, hWenuQ, hWmnuQ;
  std::vector<TH1F*> hW3e0muQ, hW2e1muQ, hW1e2muQ, hW0e3muQ;
  std::vector<TH1F*> hWTheta;

  std::vector<TH1F*> hNLElec;
  std::vector<TH1F*> hNLMuon;
  std::vector<TH1F*> hNLLeps, hNLLepsZee, hNLLepsZmm;

  std::vector<TH1F*> hNTElec;
  std::vector<TH1F*> hNTMuon;
  std::vector<TH1F*> hNTLeps;

  std::vector<TH1F*> hNJets,hNJetsZee,hNJetsZmm;
  std::vector<TH1F*> hNVtxs, hNVtxsZee, hNVtxsZmm;

  std::vector<TH1F*> hWenuCombRelIso, hWmnuCombRelIso;

  TTree* tWZCand;

//Cuts 
  typedef bool (WZAnalyzer::*CutFnPtr)() const; 
#ifndef __CINT__
  std::map<std::string,CutFnPtr> mFnPtrs_;
  std::vector<CutFnPtr> CutFns_;
#endif

};

#endif//#define _WZAnalyzer_h_
