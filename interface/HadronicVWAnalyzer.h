#ifndef _HadronicVWAnalyzer_h_
#define _HadronicVWAnalyzer_h_

#include "UserCode/CMGWPrimeGroup/interface/AnalyzerBase.h"
#include "TTree.h"

class HadronicVWAnalyzer : public AnalyzerBase {
public:
  HadronicVWAnalyzer();                         // constructor; initialize the list to be empty
  HadronicVWAnalyzer(const edm::ParameterSet & cfg, WPrimeUtil * wprimeUtil);
  ~HadronicVWAnalyzer();

  //methods for stuff to be once per job
  void setupCutOrder();

  //methods for stuff to be done for each sample
  void defineHistos(const TFileDirectory& dir);

  //methods for stuff to be done for each event
  void eventLoop(edm::EventBase const & event);
  bool passCuts(const float& weight=1.);
  void clearEvtVariables();
  void fillHistos(const int& index, const float& weight=1.);

  //methods for printers
  void printDebugEvent() const;
  void printEventDetails() const;
  void printEventLeptons() const;
  
//methods for utilities
  int countZCands(ZCandV & Zs) const;
  void calcVVariables();
  void calcWVariables();
  void calcWElecVariables();
  void calcWMuonVariables();
  void calcVWVariables();
  void calcEventVariables();

  int   calcEvtType() const;
  float calcQ() const;
  float calcGenVWInvMass() const;
  bool inEE(const TeVMuon& mu) const;

  float WLepPt() const;
  float VLepPt(int idx) const;
  
  float ElecPU(const heep::Ele & e) const;
  float MuonPU(const TeVMuon & m) const;

//methods for modifiers

//methods for the cuts
  bool passTriggersCut() const;
  bool passValidVWCut() const;

  bool passTriggerMatch(const heep::Ele& e1, const heep::Ele& e2) const;
  bool passTriggerMatch(const TeVMuon& m1, const TeVMuon& m2) const;
  bool passTriggerMatch(const pat::Electron& p, const float cut, const vstring& triggers) const;
  bool passTriggerMatch(const TeVMuon& p, const float cut, const vstring& triggers) const;
  bool passTriggerEmulation(const heep::Ele& elec, const float minPt=0.) const;

//////////////////
/////Variables////
//////////////////

  double rhoFastJet_;
  std::vector<double> effectiveElecArea_;
  std::vector<double> effectiveMuonArea_;

///My calculated qualities//////////////////
  float VWMass_;
  float Vpt_;
  float Wpt_;
  float Q_;
  uint evtType_;

// +++++++++++++++++++General Cut values

// +++++++++++++++++++W Cuts

// +++++++++++++++++++V Cuts

  //Handles
  PatElectronVH patElectronsH_;
  PatMuonVH patMuonsH_;
  JetVH patJetsH_;
  METVH metH_;
  PFCandidateVH pfCandidatesH_;

//////Chosen Candidates
  XWLeptonic vwCand_;

// +++++++++++++++++++ Histogram Definitions
  std::vector<TH1F*> hVWMass;
  std::vector<TH1F*> hVW3e0muMass, hVW2e1muMass, hVW1e2muMass, hVW0e3muMass;

  std::vector<TH1F*> hQ;
  std::vector<TH1F*> hVWTransMass;
  std::vector<TH1F*> hVWpt;
  std::vector<TH1F*> hEvtType, hEvtTypeP, hEvtTypeM;

  std::vector<TH1F*> hVMass, hVeeMass, hVmmMass ;
  std::vector<TH1F*> hV3e0muMass, hV2e1muMass, hV1e2muMass, hV0e3muMass;
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

  TTree* tVWCand;

  //Cuts
  typedef bool (HadronicVWAnalyzer::*CutFnPtr)() const; 
  std::map<std::string,CutFnPtr> mFnPtrs_;
  std::vector<CutFnPtr> CutFns_;

};

#endif//#define _HadronicVWAnalyzer_h_
