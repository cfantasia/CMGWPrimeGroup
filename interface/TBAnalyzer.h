#ifndef _TBAnalyzer_h_
#define _TBAnalyzer_h_

#include "UserCode/CMGWPrimeGroup/interface/AnalyzerBase.h"
#include "TTree.h"

class TBAnalyzer : public AnalyzerBase {
public:
  TBAnalyzer();                         // constructor; initialize the list to be empty
  TBAnalyzer(const edm::ParameterSet & cfg, int fileToRun);
  ~TBAnalyzer();

  //methods for stuff to be once per job
  void setupCutOrder();

  //methods for stuff to be done for each sample
  void defineHistos(const TFileDirectory& dir);

  //methods for stuff to be done for each event
  void eventLoop(edm::EventBase const & event);
  void clearEvtVariables();
  void fillHistos(const int& index, const float& weight=1.);

  //////Cuts//////////
  bool passValidBCut(pat::Jet & b);
  bool passValidTCut(XWLeptonic & t);
  bool passValidTBCut(XWLeptonic & tb);

  uint minNLeptons_;
  uint minNJets_;
  uint minNBJets_;

  float minMET_;

  float minWtransMass_;
  float minWpt_;
  NuAlgos nuAlgo_;

  float minTpt_;
  float minTMass_;
  float maxTMass_;

  float minBDisc_;
  std::string BDisc_;
  float maxBMass_;
  float minBpt_;

  const float MAX_MBL;

  //Selectors
  ElectronSelector looseElectron_;
  MuonSelector looseMuon_;
  JetSelector looseJet_;

  //Handles
  JetVH patJetsH_;


//////Chosen Candidates
  ElectronV allElectrons_, looseElectrons_, tightElectrons_;
  MuonV allMuons_, looseMuons_, tightMuons_;
  JetV looseJets_, looseBJets_, tightBJets_;
  pat::MET met_;

  pat::Jet bCand1_;
  pat::Jet bCand2_;
  WCandidate wCand_;
  XWLeptonic tCand_;
  XWLeptonic tbCand_;

  uint runNumber_, lumiNumber_, evtNumber_;
  uint evtType_;
  float TBMass_, TMass_, Tpt_;
  float BMass1_, BDisc1_, Bpt1_;
  float BMass2_, BDisc2_, Bpt2_;
  float WTransMass_, Wpt_;
  float Q_;
  float MET_, METSig_;
  uint NVtxs_;
  float Mbl_;

// +++++++++++++++++++ Histogram Definitions
  std::vector<TH1F*> hTBMass, hTBenuMass, hTBmnuMass;
  std::vector<TH1F*> hTMass, hTenuMass, hTmnuMass, hQ;
  std::vector<TH1F*> hEvtType;
  std::vector<TH1F*> hTpt, hTenupt, hTmnupt;

  std::vector<TH1F*> hMET, hMETSig;
  std::vector<TH1F*> hWTransMass, hWenuTransMass, hWmnuTransMass;
  std::vector<TH1F*> hWpt;
  std::vector<TH1F*> hB1Mass, hB1Disc, hB1pt;
  std::vector<TH1F*> hB2Mass, hB2Disc, hB2pt;
  std::vector<TH1F*> hMbl;
  std::vector<TH1F*> hNLElec, hNLMuon, hNLLeps;
  std::vector<TH1F*> hNLJets, hNLBJets, hNTBJets;
  std::vector<TH1F*> hWeight;

  TTree* tTBCand;

};

#endif//#define _TBAnalyzer_h_
