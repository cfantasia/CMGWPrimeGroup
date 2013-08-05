#ifndef _TagAndProbe_h_
#define _TagAndProbe_h_

#include "UserCode/CMGWPrimeGroup/interface/AnalyzerBase.h"
#include "TTree.h"

class TagAndProbe : public AnalyzerBase {
public:
  TagAndProbe();                         // constructor; initialize the list to be empty
  TagAndProbe(const edm::ParameterSet & cfg, int fileToRun);
  ~TagAndProbe();

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
  void calcEventVariables();

  int   calcEvtType() const;
  float calcLeadPt(int type=0) const;
  bool inEE(const TeVMuon& mu) const;

  float ElecPU(const heep::Ele & e) const;
  float MuonPU(const TeVMuon & m) const;

//methods for modifiers

//methods for the cuts
  bool passNumberOfZsCut() const;
  bool passValidZCut(ZCandidate& z);

  bool passZLepPtCut() const;

  bool passZLepTightCut(bool firstDaughter);

  bool passTriggerEmulation(const heep::Ele& elec, const float minPt=0.) const;

//////////////////
/////Variables////
//////////////////

  bool removeTauEvents_, adjustMETPhi_;
  float elScaleFactor_, muScaleFactor_;

///My calculated qualities//////////////////
  uint runNumber_;
  uint lumiNumber_;
  uint evtNumber_;
  float ZMass_, Zpt_, ZDr_;
  uint numZs_, evtType_; 
  uint NVtxs_;
  float LeadPt_;
  float LeadElecPt_;
  float LeadMuonPt_;
  float MET_;
  int ZTightCode_;
  float ZLep1Pt_, ZLep1Eta_, ZLep1Phi_, ZLep2Pt_, ZLep2Eta_, ZLep2Phi_, WLepPt_, WLepEta_, WLepPhi_;
  float ZLep1PtGen_, ZLep1EtaGen_, ZLep1PhiGen_, ZLep2PtGen_, ZLep2EtaGen_, ZLep2PhiGen_;

// +++++++++++++++++++General Cut values
  uint minNLeptons_, maxNVLLeptons_, maxNJets_;
  uint minNTightLeptons_;
  uint maxNumZs_;
  float minLeadPt_;
  float minMET_;

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
  ElectronSelector countElectron_, TagElectron_, looseProbeElectron_, tightProbeElectron_;
  MuonSelector     countMuon_,     TagMuon_,     looseProbeMuon_,     tightProbeMuon_;
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
  ElectronV allElectrons_, countElectrons_, TagElectrons_, looseProbeElectrons_, tightProbeElectrons_;
  MuonV     allMuons_,     countMuons_,     TagMuons_,     looseProbeMuons_,     tightProbeMuons_;
  JetV looseJets_;
  pat::MET met_;

  ZCandidate zCand_;

// +++++++++++++++++++ Histogram Definitions

  std::vector<TTree*> tEvts;

};

#endif//#define _TagAndProbe_h_
