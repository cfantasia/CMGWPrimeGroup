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
  void defineResolutionHistos(const TFileDirectory & dir, float Mass);

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
  float calcLt() const;
  float calcTriLepMass() const;
  float calcGenWZInvMass() const;
  bool inEE(const TeVMuon& mu) const;

  float ElecPU(const heep::Ele & e) const;
  float MuonPU(const TeVMuon & m) const;

//methods for modifiers

//methods for the cuts
  bool passLeadingLeptonPtCut() const;
  bool passNumberOfZsCut() const;
  bool passValidWCut(WCandidate& w);
  bool passValidZCut(ZCandidate& z);
  bool passValidWZCut(XWLeptonic& xw);
  bool passLtCut() const;

  bool passZLepPtCut() const;
  bool passZeePtCut() const;
  bool passZmumuPtCut() const;

  bool passWLepPtCut() const;
  bool passWLepIsoCut() const;

  bool passWLepTightCut();
  bool passZLepTightCut(bool firstDaughter);
  bool passZLepNOTTightCut(bool firstDaughter);
  bool passWFlavorElecCut() const;
  bool passWFlavorMuonCut() const;
  bool passFakeEvtCut() const;
  bool passMaxWtransMassCut(const WCandidate& w, const float& cut) const;
  bool passMaxDeltaWtransMassCut(const float& cut) const;
  bool passFakeLeptonProbeCut();

  bool passTriggerEmulation(const heep::Ele& elec, const float minPt=0.) const;

//////////////////
/////Variables////
//////////////////

  std::vector<double> effectiveElecArea_;
  std::vector<double> effectiveMuonArea_;
  
  bool doSystematics_, doMatrix_, removeTauEvents_, adjustMETPhi_;
  float elScaleFactor_, muScaleFactor_;

///My calculated qualities//////////////////
  uint runNumber_;
  uint lumiNumber_;
  uint evtNumber_;
  float WZMass_, WprimeGenMass_;
  float Zpt_, ZDr_;
  float ZMass_;
  float Wpt_;
  float WTransMass_;
  float WCharge_;
  float Lt_;
  float TriLepMass_;
  float Q_;
  float Discriminant_;
  float MET_, METPhi_, METSig_;
  uint evtType_;
  uint numZs_; 
  uint NVtxs_;
  float LeadPt_;
  float LeadElecPt_;
  float LeadMuonPt_;
  bool TT, TF, FT, FF;
  int ZTightCode_, WTightCode_;
  float ZLep1Pt_, ZLep1Eta_, ZLep1Phi_, ZLep2Pt_, ZLep2Eta_, ZLep2Phi_, WLepPt_, WLepEta_, WLepPhi_;
  float ZLep1PtGen_, ZLep1EtaGen_, ZLep1PhiGen_, ZLep2PtGen_, ZLep2EtaGen_, ZLep2PhiGen_;
  float WLepPtGen_, WLepEtaGen_, WLepPhiGen_, WNeuPtGen_, WNeuEtaGen_, WNeuPhiGen_;

// +++++++++++++++++++General Cut values
  uint minNLeptons_, maxNVLLeptons_, maxNJets_;
  uint minNTightLeptons_;
  uint maxNumZs_;
  float minLeadPt_;
  float minMET_;

// +++++++++++++++++++Lt Cuts
  float minLt_;

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
  ElectronSelector extraElectron_, looseZElectron_, tightZElectron_, looseWElectron_, tightWElectron_;
  MuonSelector     extraMuon_, looseZMuon_, tightZMuon_, looseWMuon_, tightWMuon_;
  JetSelector looseJet_;

  //Handles
  PatElectronVH patElectronsH_;
  PatMuonVH patMuonsH_;
  METVH metH_;
  PFCandidateVH pfCandidatesH_;
  edm::Handle<double> rhoFastJetH_;
  edm::Handle<std::vector<pat::Jet> > patJetsH_;
  edm::Handle<std::vector<reco::Vertex> > verticesH_;

  //Input Tags
  edm::InputTag rhoFastJetLabel_;

//////Chosen Candidates
  ElectronV allElectrons_, extraElectrons_, looseZElectrons_, tightZElectrons_, looseWElectrons_, tightWElectrons_;
  MuonV     allMuons_, extraMuons_, looseZMuons_, tightZMuons_, looseWMuons_, tightWMuons_;
  JetV looseJets_;
  pat::MET met_;

  ZCandidate zCand_;
  WCandidate wCand_;
  XWLeptonic wzCand_;
  NuAlgos wzAlgo_;

// +++++++++++++++++++ Histogram Definitions
  std::vector<TH1F*> hWZMass;
  std::vector<TH1F*> hWZ3e0mMass, hWZ2e1mMass, hWZ1e2mMass, hWZ0e3mMass;
  std::vector<TH1F*> hRes, hWprimeGenMass;

  std::vector<TH1F*> hQ;
  std::vector<TH1F*> hWZTransMass;
  std::vector<TH1F*> hWZpt;
  std::vector<TH1F*> hLt,hLtee,hLtmm,hLt3e0m,hLt2e1m,hLt1e2m,hLt0e3m;
  std::vector<TH1F*> hTriLepMass;
  std::vector<TH1F*> hEvtType, hEvtTypeP, hEvtTypeM;
  std::vector<TH1F*> hLeadPt, hLeadPtZee, hLeadPtZmm;
  std::vector<TH1F*> hLeadElecPt, hLeadMuonPt;

  std::vector<TH1F*> hZMass, hZeeMass, hZmmMass, hZ3e0mMass, hZ2e1mMass, hZ1e2mMass, hZ0e3mMass;
  std::vector<TH1F*> hZeeMassTT, hZeeMassTF, hZmmMassTT, hZmmMassTF;
  std::vector<TH1F*> hZpt,hZeept,hZmmpt, hZ3e0mpt, hZ2e1mpt, hZ1e2mpt, hZ0e3mpt;
  std::vector<TH1F*> hZeeDr,hZmmDr;


  std::vector<TH1F*> hMET, hMETee, hMETmm, hMET3e0m, hMET2e1m, hMET1e2m, hMET0e3m, hMETSig, hMETPhi;

  std::vector<TH1F*> hWTransMass, hWenuTransMass, hWmnuTransMass, hW3e0mTransMass, hW2e1mTransMass, hW1e2mTransMass, hW0e3mTransMass;
  std::vector<TH1F*> hWpt, hWenupt, hWmnupt, hW3e0mpt, hW2e1mpt, hW1e2mpt, hW0e3mpt;
  std::vector<TH1F*> hWQ, hWenuQ, hWmnuQ, hW3e0mQ, hW2e1mQ, hW1e2mQ, hW0e3mQ;
  std::vector<TH1F*> hWTheta;

  std::vector<TH1F*> hNLElec, hNLMuon, hNLLeps;

  std::vector<TH1F*> hNTElec, hNTMuon, hNTLeps;

  std::vector<TH1F*> hNJets;
  std::vector<TH1F*> hNVtxs, hNVtxsZee, hNVtxsZmm;
  std::vector<TH1F*> hWeight;
  std::vector<TH1F*> hL1FastJet;

  std::vector<TH1F*> hWenuCombRelIso, hWmnuCombRelIso;

  std::vector<TH1F*> hDeltaPhiWJet, hDeltaPhiLepJet, hDeltaRLepJet, hDeltaWMT;

  std::vector<TH2F*> hEtaVsPt, hEtaVsPtElec, hEtaVsPtMuon, hEtaVsPt3e0m, hEtaVsPt2e1m, hEtaVsPt1e2m, hEtaVsPt0e3m;
  std::vector<TH2F*> hPtVsZeeBarrelMassTT,hPtVsZeeBarrelMassTF,hPtVsZeeEndCapMassTT,hPtVsZeeEndCapMassTF,hPtVsZmmMassTT,hPtVsZmmMassTF;

  TH1F *hDiscriminant, *hDiscriminantFrac, *hDiscriminantAngle, *hDiscriminantReal, *hDiscriminantImag,;
  TH1F* hVtxMatch;

  std::vector<TTree*> tEvts;

};

#endif//#define _WZAnalyzer_h_
