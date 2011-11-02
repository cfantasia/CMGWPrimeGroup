#ifndef _AnalyzerBase_h_
#define _AnalyzerBase_h_

#include "UserCode/CMGWPrimeGroup/interface/util.h"
#include "UserCode/CMGWPrimeGroup/interface/WPrimeUtil.h"
#include "UserCode/CMGWPrimeGroup/interface/BosonFinder.h"
#include "UserCode/CMGWPrimeGroup/interface/TeVMuon.h"

class AnalyzerBase {
public:

  AnalyzerBase();
  AnalyzerBase(const edm::ParameterSet & cfg, WPrimeUtil * wprimeUtil);
  virtual ~AnalyzerBase();

///////////////Utilities//////////////////

//Tabulate results after the cut has been passed
  virtual void tabulateEvent(const int& cut_index, const float& weight);
  virtual void tabulateFile(std::vector<wprime::InputFile>::const_iterator fi, wprime::EffV& results);
  virtual void printFileSummary(std::vector<wprime::InputFile>::const_iterator fi, ofstream& out);
 
  virtual void fillHistos(const int& index, const float& weight=1.) = 0;//Pure Virtual
  virtual void defineHistos(const TFileDirectory & dir);
  virtual void defineHistoset(const std::string& n, const std::string& t, const std::string& xtitle,
                               const int& nbins, const float& min, const float& max, const std::string& units,
                               std::vector<TH1F*>& h, const TFileDirectory& d);
  virtual void resetcounters();
  virtual void clearEvtVariables();

  //methods for printers
  virtual void printEventFull(edm::EventBase const & event) const;
  virtual void printPassingEvent(edm::EventBase const & event);
  virtual void printDebugEvent() const;
  virtual void printEventToFile(edm::EventBase const & event);
  virtual void printEventDetails() const;
  virtual void printEventLeptons() const;
  virtual void printLeptons() const;
  virtual void printElectrons() const;
  virtual void printMuons() const;
  virtual void printJets() const;
  virtual void printElectron(const pat::Electron& elec) const;
  virtual void printElectron(const heep::Ele& elec) const;
  virtual void printMuon(const TeVMuon& mu) const;
  virtual void printJet(const pat::Jet& jet) const;


////////////////////////////
//////////setters///////////
////////////////////////////
  void setCandEvtFile(const std::string & s);

////////////////////
//////Cuts//////////
////////////////////
  virtual bool passCuts(const float& weight);

  virtual bool passNoCut() const;
  virtual bool passTriggersCut() const;
  virtual bool passMinNLeptonsCut() const;
  virtual bool passMinNTightLeptonsCut() const;
  virtual bool passMaxNLeptonsCut() const;
  virtual bool passMinNJetsCut() const;

  virtual bool passMinMETCut() const;
  virtual bool passMinPtCut(const reco::Candidate& cand, const float& cut) const;

/////////Check Z Properties/////
  virtual bool passValidZCut() const;
  virtual bool passZMassCut() const;
  virtual bool passZptCut() const;

/////////Check W Properties/////
  virtual bool passValidWCut() const;
  virtual bool passWtransMassCut() const;
  virtual bool passWptCut() const;

///////Check Had. V Properties//
  virtual bool passValidVCut() const;
  virtual bool passVMassCut() const;
  virtual bool passVptCut() const;

//////////////////
//file stuff//////
//////////////////
  
  virtual void beginFile(std::vector<wprime::InputFile>::const_iterator fi);
  virtual void eventLoop(edm::EventBase const & event);

// operations to be done when closing input file 
  virtual void endFile(std::vector<wprime::InputFile>::const_iterator fi,
                       ofstream & out);
  virtual void endAnalysis(ofstream & out);

protected:

  bool debugme;//print stuff if active
  bool doPreselect_;

  std::string looseElectronType_, tightElectronType_;
  std::string looseMuonType_, tightMuonType_;
  std::string looseJetType_, tightJetType_;

  edm::InputTag electronsLabel_;
  edm::InputTag muonsLabel_;
  edm::InputTag jetsLabel_;
  edm::InputTag metLabel_;
  edm::InputTag pfCandsLabel_;
  edm::InputTag vertexLabel_;
  edm::InputTag hltEventLabel_;
  vstring triggersToUse_;

  WPrimeUtil * wprimeUtil_;
  ofstream outCandEvt_;

  uint muReconstructor_;
  bool useAdjustedMET_;

  TH1F * hNumEvts;

  float weight_;

  int NCuts_;
  std::vector<std::string> Cuts_;

  Pset eSelectorPset_;
  ElectronSelector looseElectron_;
  ElectronSelector tightElectron_;
  pat::strbitset electronLooseResult_;
  pat::strbitset electronTightResult_;

  Pset mSelectorPset_;
  MuonSelector looseMuon_;
  MuonSelector tightMuon_;
  pat::strbitset muonLooseResult_;
  pat::strbitset muonTightResult_;

  Pset jSelectorPset_;
  JetSelector looseJet_;
  JetSelector tightJet_;
  pat::strbitset jetLooseResult_;

  PatElectronVH patElectronsH_;
  PatMuonVH patMuonsH_;
  METVH metH_;
  PFCandidateVH pfCandidatesH_;
  std::vector<reco::Vertex>  vertices_;

//////Chosen Candidates
  ElectronV allElectrons_, looseElectrons_, tightElectrons_;
  MuonV allMuons_, looseMuons_, tightMuons_;
  JetV  allJets_, looseJets_, tightJets_;
  pat::MET met_;
  ZCandidate zCand_;
  WCandidate wCand_, vCand_;

  pat::TriggerEvent triggerEvent_; 
  wprime::EffV results_;

  float minNLeptons_;
  float maxNLeptons_;
  float minNTightLeptons_;
  float minNJets_;

  float minMET_;

  float minZmass_;
  float maxZmass_;
  float minVmass_;
  float maxVmass_;
  float minWtransMass_;

  float minZpt_;
  float minWpt_;
  float minVpt_;

  std::vector<unsigned> muon_reconstructors;

};

#endif//#define _AnalyzerBase_h_

