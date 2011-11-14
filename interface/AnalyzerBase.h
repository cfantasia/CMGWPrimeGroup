#ifndef _AnalyzerBase_h_
#define _AnalyzerBase_h_

#include "UserCode/CMGWPrimeGroup/interface/util.h"
#include "UserCode/CMGWPrimeGroup/interface/WPrimeUtil.h"
#include "UserCode/CMGWPrimeGroup/interface/BosonFinder.h"
#include "UserCode/CMGWPrimeGroup/interface/TeVMuon.h"

class AnalyzerBase{
public:

  AnalyzerBase();
  AnalyzerBase(const edm::ParameterSet & cfg, int fileToRun);
  virtual ~AnalyzerBase();

///////////////Utilities//////////////////

//Fill Vector of Cuts based on map
  void fillCuts();

//Tabulate results after the cut has been passed
  virtual void tabulateEvent(const int& cut_index, const float& weight);
  virtual void tabulateFile(std::vector<wprime::InputFile>::const_iterator fi, wprime::EffV& results);
  virtual void printFileSummary(std::vector<wprime::InputFile>::const_iterator fi, ofstream& out);
 
  virtual void fillHistos(const int& index, const float& weight=1.) = 0;//Pure Virtual
  virtual void defineHistos(const TFileDirectory & dir);
  virtual void defineResolutionHistos(const TFileDirectory & dir, float Mass){}
  virtual void defineHistoset(const std::string& n, const std::string& t, 
			      const std::string& xtitle, int nbins, 
			      float xmin, float xmax, 
			      const std::string& units,
                               std::vector<TH1F*>& h, const TFileDirectory& d);
  void defineOneHisto(const std::string & name, const std::string & title, 
		      const std::string & xtitle, int nbins, float xmin, float xmax,
		      const std::string & units,TH1F* & h,const TFileDirectory & d);

  // mass format expected in <x>.<y> TeV, e.g. "1.2", corresponding to 1.2 TeV
  // channel could be "e", "mu", "ee", "mumu", etc
  void createResolutionHist(const TFileDirectory & d, float Mass,
			    const std::string & channel, TH1F* & put_here);

  // mass format expected in <x>.<y> TeV, e.g. "1.2", corresponding to 1.2 TeV
  // channel could be "e", "mu", "ee", "mumu", etc
  void createGenMtHist(const TFileDirectory & d, float Mass,
		       const std::string & channel, TH1F* & put_here);


  virtual void resetCounters();
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

////////////////////
//////Cuts//////////
////////////////////
  virtual bool passCuts(const float& weight);

  virtual bool passNoCut() const;
  virtual bool passTriggersCut() const;
  virtual bool passMinNLeptonsCut(const ElectronV& electrons, const MuonV& muons, const float & cut) const;
  virtual bool passMaxNLeptonsCut(const ElectronV& electrons, const MuonV& muons, const float & cut) const;
  virtual bool passMinNJetsCut(const JetV& jets, const float & cut) const;

  virtual bool passMinMETCut(const pat::MET & met, const float& cut) const;
  virtual bool passMinPtCut(const reco::Candidate& cand, const float& cut) const;

/////////Check Z Properties/////
  virtual bool passValidZCut(const ZCandidate& z) const;
  virtual bool passZMassCut(const ZCandidate& z, const float& mincut, const float& maxcut) const;
  virtual bool passZptCut(const ZCandidate& z, const float& cut) const;

/////////Check W Properties/////
  virtual bool passValidWCut(const WCandidate& w) const;
  virtual bool passWtransMassCut(const WCandidate& w, const float& cut) const;
  virtual bool passVMassCut(const WCandidate& w, const float& mincut, const float& maxcut) const;
  virtual bool passWptCut(const WCandidate& w, const float& cut) const;

//////////////////
//file stuff//////
//////////////////
  
  void run();
  virtual void beginFile(std::vector<wprime::InputFile>::iterator fi);
  virtual void eventLoop(edm::EventBase const & event);
  
  void setEventWeight(edm::EventBase const & event);


// operations to be done when closing input file 
  virtual void endFile(std::vector<wprime::InputFile>::iterator fi,
                       ofstream & out);
  virtual void endAnalysis(ofstream & out);

protected:

  // print out event # 
  unsigned int reportAfter_;
  // maximum # of events to process (set to <0 for processing all events)
  int maxEvents_;
  // Should we use the json file
  bool useJSON_;
  // Should we count the number of gen evts in patTuple?
  bool countGenEvts_;

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

  ofstream outLogFile_;
  ofstream outCandEvtFile_;
  std::vector<wprime::InputFile> inputFiles_; 

  std::string outputFile_;
  std::string logFile_;
  std::string candEvtFile_;

  //Variables for counting # of gen events
  std::vector<uint> nEvents_;
  std::vector<std::string> ctrNames_;

  edm::InputTag genLabel_;
  edm::InputTag pfLabel_;
  edm::InputTag pileupLabel_;

  edm::Handle<std::vector< PileupSummaryInfo > > PupH_;

  float lumi_ipb; // in pb^-1, to be retrieved from samples_cross_sections.txt

  int fileToRun;
  fwlite::TFileService * fs;
  // directory containing all input samples
  std::string top_level_dir; 

  // file with samples & cross-sections
  std::string sample_cross_sections;

  std::vector< edm::LuminosityBlockRange > jsonVector;

  bool jsonContainsEvent (const std::vector<edm::LuminosityBlockRange>&jsonVec,
                          const edm::EventBase &event);

  std::string MCPUDistFile_;
  std::string MCPUDistHist_;
  std::string DataPUDistFile_;
  std::string DataPUDistHist_;

  WPrimeUtil* wprimeUtil_;

  bool doRecoilCorrectionForW_;

  uint muReconstructor_;
  bool useAdjustedMET_;

  TH1F * hNumEvts;

  float weight_;

  int NCuts_;
  std::vector<std::string> CutNames_;
  std::vector<std::string> CutDescs_;
  typedef boost::function<bool()> fnCut;
  std::map<std::string,fnCut > mFnPtrs_;
  std::vector<fnCut > CutFns_;

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

  int NumSigSamples;
  
  // return index of current signal sample (-1 if this is not a signal sample)
  int getSignalSampleIndex();

 private:
  int signalSample_index;
  void setNumSignalSamples();
};

#endif//#define _AnalyzerBase_h_

