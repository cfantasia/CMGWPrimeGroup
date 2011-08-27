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
  virtual void Tabulate_Me(const int& cut_index, const float& weight);
  virtual void Fill_Histos(const int& index, const float& weight=1.);
  virtual void Declare_Histos(const TFileDirectory & dir);
  virtual void DeclareHistoSet(std::string n, std::string t, std::string xtitle,
                               int nbins, float min, float max, std::string units,
                               std::vector<TH1F*>& h, TFileDirectory& d);
  virtual void ResetCounters();

  //methods for printers
  virtual void PrintEventFull(edm::EventBase const & event) const;
  virtual void PrintPassingEvent(edm::EventBase const & event);
  virtual void PrintDebugEvent() const;
  virtual void PrintEventToFile(edm::EventBase const & event);
  virtual void PrintEventDetails() const;
  virtual void PrintEventLeptons() const;
  virtual void PrintLeptons() const;
  virtual void PrintElectrons() const;
  virtual void PrintMuons() const;
  virtual void PrintJets() const;
  virtual void PrintElectron(const pat::Electron& elec) const;
  virtual void PrintElectron(const heep::Ele& elec) const;
  virtual void PrintMuon(const TeVMuon& mu) const;
  virtual void PrintJet(const pat::Jet& jet) const;


////////////////////////////
//////////Setters///////////
////////////////////////////
  void SetCandEvtFile(const std::string & s);

////////////////////
//////Cuts//////////
////////////////////
  virtual bool PassNoCut() const;
  virtual bool PassTriggersCut() const;
  virtual bool PassMinNLeptonsCut() const;
  virtual bool PassMinNTightLeptonsCut() const;
  virtual bool PassMaxNLeptonsCut() const;
  virtual bool PassMinNJetsCut() const;

  virtual bool PassMinMETCut() const;

/////////Check Z Properties/////
  virtual bool PassValidZCut() const;
  virtual bool PassZMassCut() const;
  virtual bool PassZptCut() const;

/////////Check W Properties/////
  virtual bool PassValidWCut() const;
  virtual bool PassWtransMassCut() const;
  virtual bool PassWptCut() const;

///////Check Had. V Properties//
  virtual bool PassValidVCut() const;
  virtual bool PassVMassCut() const;
  virtual bool PassVptCut() const;

//////////////////
//file stuff//////
//////////////////
  virtual void beginFile(std::vector<wprime::InputFile>::const_iterator fi);

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

  edm::InputTag hltEventLabel_;
  edm::InputTag pileupLabel_;
  edm::InputTag vertexLabel_;
  vstring triggersToUse_;

  WPrimeUtil * wprimeUtil_;
  ofstream outCandEvt_;

  uint muonAlgo_;
  bool useAdjustedMET_;

  TH1F * hNumEvts;

  int NCuts_;
  std::vector<std::string> Cuts_;

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
  JetSelector tightJet_;
  pat::strbitset jetResult_;

  PatElectronVH patElectronsH_;
  PatMuonVH patMuonsH_;
  METVH metH_;
  PFCandidateVH pfCandidatesH_;

//////Chosen Candidates
  ElectronV electrons_, looseElectrons_, tightElectrons_;
  MuonV muons_, looseMuons_, tightMuons_;
  JetV  jets_, looseJets_, tightJets_;
  pat::MET met_;
  ZCandidate zCand_;
  WCandidate wCand_, vCand_;

  pat::TriggerEvent triggerEvent_; 
  std::vector< PileupSummaryInfo > PupInfo_; 
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
};

#endif//#define _AnalyzerBase_h_

