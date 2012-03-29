#ifndef _HadronicVWAnalyzer_h_
#define _HadronicVWAnalyzer_h_

#include "UserCode/CMGWPrimeGroup/interface/AnalyzerBase.h"
#include "TTree.h"

class HadronicVWAnalyzer : public AnalyzerBase {
public:
  HadronicVWAnalyzer();                         // constructor; initialize the list to be empty
  HadronicVWAnalyzer(const edm::ParameterSet & cfg, int fileToRun);
  ~HadronicVWAnalyzer();

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
  void calcVVariables();
  void calcWVariables();
  void calcVWVariables();
  void calcEventVariables();

  int   calcEvtType() const;
  float calcQ() const;

//methods for modifiers

//methods for the cuts
  bool makeAndPassValidVCut();
  bool makeAndPassValidWCut();
  bool makeAndPassValidVWCut();

//////////////////
/////Variables////
//////////////////

///My calculated qualities//////////////////
  uint runNumber_, lumiNumber_, evtNumber_;
  float VWMass_;
  float VMass_, Vpt_;
  float WTransMass_, Wpt_;
  float Q_;
  float MET_, METSig_;
  uint NVtxs_;
  uint evtType_;

// +++++++++++++++++++General Cut values
  uint minNLeptons_;
  uint minNJets_;

  float minMET_;

// +++++++++++++++++++W Cuts
  float minWtransMass_;
  float minWpt_;

// +++++++++++++++++++V Cuts
  float minVmass_;
  float maxVmass_;
  float minVpt_;

  //Selectors
  ElectronSelector looseElectron_;
  MuonSelector looseMuon_;
  JetSelector looseJet_;

  //Handles
  PatElectronVH patElectronsH_;
  PatMuonVH patMuonsH_;
  JetVH patJetsH_;
  METVH metH_;
  PFCandidateVH pfCandidatesH_;
  edm::Handle<std::vector<reco::Vertex> > verticesH_;


//////Chosen Candidates
  ElectronV allElectrons_, looseElectrons_;
  MuonV allMuons_, looseMuons_;
  JetV  allJets_, looseJets_;
  pat::MET met_;

  NuAlgos nuAlgo_;
  ZCandidate vCand_;
  WCandidate wCand_;
  XWLeptonic vwCand_;

// +++++++++++++++++++ Histogram Definitions
  std::vector<TH1F*> hVWMass, hVWenuMass, hVWmnuMass;

  std::vector<TH1F*> hQ;
  std::vector<TH1F*> hVWTransMass;
  std::vector<TH1F*> hVWpt;
  std::vector<TH1F*> hEvtType, hEvtTypeP, hEvtTypeM;

  std::vector<TH1F*> hVMass;
  std::vector<TH1F*> hVpt;

  std::vector<TH1F*> hMET, hMETenu, hMETmnu, hMETSig;

  std::vector<TH1F*> hWTransMass, hWenuTransMass, hWmnuTransMass;
  std::vector<TH1F*> hWpt,hWenupt,hWmnupt;
  std::vector<TH1F*> hWQ, hWenuQ, hWmnuQ;

  std::vector<TH1F*> hNLElec;
  std::vector<TH1F*> hNLMuon;
  std::vector<TH1F*> hNLLeps;

  std::vector<TH1F*> hNJets;
  std::vector<TH1F*> hNVtxs;
  std::vector<TH1F*> hWeight;
  TTree* tVWCand;
};

#endif//#define _HadronicVWAnalyzer_h_
