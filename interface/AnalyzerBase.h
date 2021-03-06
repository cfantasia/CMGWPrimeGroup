#ifndef _AnalyzerBase_h_
#define _AnalyzerBase_h_

#include "UserCode/CMGWPrimeGroup/interface/util.h"
#include "UserCode/CMGWPrimeGroup/interface/selectors.h"
#include "UserCode/CMGWPrimeGroup/interface/WPrimeUtil.h"
#include "UserCode/CMGWPrimeGroup/interface/BosonFinder.h"
#include "UserCode/CMGWPrimeGroup/interface/TeVMuon.h"
#include "TH2.h"

class AnalyzerBase{
public:

  AnalyzerBase();
  AnalyzerBase(const edm::ParameterSet & cfg, int fileToRun);
  virtual ~AnalyzerBase();

///////////////Utilities//////////////////

//Fill Vector of Cuts based on map
  typedef boost::function<bool()> fnCut;
  void fillCuts(const std::map<std::string,fnCut >& mFnPtrs);

//Tabulate results after the cut has been passed
  virtual void tabulateEvent(const int& cut_index, const float& weight);
  virtual void tabulateFile(std::vector<wprime::InputFile>::iterator fi);
  virtual void printFileSummary(std::vector<wprime::InputFile>::const_iterator fi, ofstream& out);
 
  virtual void fillHistos(const int& index, const float& weight=1.) = 0;//Pure Virtual
  virtual void defineHistos(const TFileDirectory & dir) = 0;
  virtual void defineResolutionHistos(const TFileDirectory & dir, float Mass){}
  virtual void defineHistoSet(const std::string& n, const std::string& t, 
			      const std::string& xtitle, int nbins, 
			      float xmin, float xmax, 
			      const std::string& units,
                               std::vector<TH1F*>& h, const TFileDirectory& d);
  virtual void defineHistoSet(const std::string& n, const std::string& t, 
			      const std::string& xtitle, int nxbins, 
			      float xmin, float xmax, 
			      const std::string& ytitle, int nybins, 
			      float ymin, float ymax, 
                              std::vector<TH2F*>& h, const TFileDirectory& d);
  void defineOneHisto(const std::string & name, const std::string & title, 
		      const std::string & xtitle, int nbins, float xmin, float xmax,
		      const std::string & units,TH1F* & h,const TFileDirectory & d);
  void defineOneHisto(const std::string & name, const std::string & title, 
		      const std::string & xtitle, int nxbins, float xmin, float xmax,
          const std::string & ytitle, int nybins, float ymin, float ymax,
		      TH2F* & h,const TFileDirectory & d);

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
  virtual const reco::Vertex& findPV(const std::vector<reco::Vertex>& vtxs) const;
  virtual const reco::Vertex& findPV(const std::vector<reco::Vertex>& vtxs, const reco::Candidate& p) const;
  virtual bool sameVertex(const reco::Vertex& vtx, const reco::Candidate& p) const;

  //methods for printers
  virtual void printEventFull(edm::EventBase const & event) const;
  virtual void printPassingEvent(edm::EventBase const & event);
  virtual void printDebugEvent() const;
  virtual void printEventDetails() const;

  template<class T>
    void print(const std::vector<T> & particles, 
               const std::vector<bool> & mask) const{
    print(particles, "objects", mask);
    
  }

  template<class T>
    void print(const std::vector<T> & particles, const std::string & name="objects",
               const std::vector<bool> & mask=std::vector<bool>()) const{
    //Mask let's you ignore certain objects
    bool useMask = mask.size() == particles.size(); 
    std::cout<<"----There are "<<particles.size()<<" "<<name<<" ------\n";
    for(uint i=0; i<particles.size(); ++i){
      if(useMask && !mask[i]) continue;
      print(particles[i]);
      std::cout<<"   ------------   \n";
    }
  }
  
  virtual void print(const pat::Electron& elec) const;
  virtual void print(const heep::Ele& elec) const;
  virtual void print(const TeVMuon& mu) const;
  virtual void print(const pat::Jet& jet) const;


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
  virtual bool passMaxNJetsCut(const JetV& jets, const float & cut) const;

  virtual bool passMinMETCut(const pat::MET & met, const float& cut) const;
  virtual bool passMinPtCut(const reco::Candidate& cand, const float& cut) const;

/////////Check Z Properties/////
  virtual bool passValidZCut(const ZCandidate& z) const;
  virtual bool passZMassCut(const ZCandidate& z, const float& mincut, const float& maxcut) const;
  virtual bool passZptCut(const ZCandidate& z, const float& cut) const;

/////////Check W Properties/////
  virtual bool passValidWCut(const WCandidate& w) const;
  virtual bool passWtransMassCut(const WCandidate& w, const float& cut) const;
  virtual bool passWptCut(const WCandidate& w, const float& cut) const;

/////////Check XW Properties/////
  virtual bool passValidXWCut(const XWLeptonic& w) const;
  virtual bool passXWMassCut(const XWLeptonic& w, const NuAlgos & algo, const float& min, const float& max) const;
  virtual bool passXWptCut(const XWLeptonic& w, const float& cut) const;

//////////////////
//file stuff//////
//////////////////
  
  void run();
  virtual void beginFile(std::vector<wprime::InputFile>::iterator fi);
  virtual void eventLoop(edm::EventBase const & event) = 0;
  
// operations to be done when closing input file 
  virtual void endFile(std::vector<wprime::InputFile>::iterator fi,
                       ofstream & out);
  virtual void endAnalysis(ofstream & out);

protected:  //These are available to derived classes
  bool debug_;//print stuff if active
  bool doPreselect_;

  //Input Tags
  edm::InputTag electronsLabel_;
  edm::InputTag muonsLabel_;
  edm::InputTag jetsLabel_;
  edm::InputTag metLabel_;
  edm::InputTag pfCandsLabel_;
  edm::InputTag vertexLabel_;
  edm::InputTag hltEventLabel_;
  edm::InputTag genLabel_;
  edm::InputTag pileupLabel_;

  vstring triggersToUse_;

  ofstream outLogFile_;
  ofstream outCandEvtFile_;
  std::vector<wprime::InputFile> inputFiles_; 

  edm::Handle<std::vector< PileupSummaryInfo > > PupH_;
  edm::Handle<pat::TriggerEvent> triggerEventH_;

  fwlite::TFileService * fs;

  WPrimeUtil* wprimeUtil_;

  bool doRecoilCorrectionForW_;

  uint muReconstructor_;
  bool useAdjustedMET_;

  float weight_;

  //These will get merged into a struct soon
  int NCuts_;
  std::vector<std::string> CutNames_;
  std::vector<std::string> CutDescs_;
  std::vector<fnCut > CutFns_;
  TH1F * hNumEvts;

  //These should be done by the each analysis though
  PatElectronVH patElectronsH_;
  PatMuonVH patMuonsH_;
  METVH metH_;
  PFCandidateVH pfCandidatesH_;

  std::vector<unsigned> muon_reconstructors;

  int NumSigSamples;
  
  // return index of current signal sample (-1 if this is not a signal sample)
  int getSignalSampleIndex();

 private:  //These are not directly available to derived classes

  // print out event # (set to <0 for percentage of events)
  int reportAfter_;
  // maximum # of events to process (set to <0 for processing all events)
  int maxEvents_;
  // Should we use the json file
  bool useJSON_;

  int signalSample_index;
  void setNumSignalSamples();

  std::vector< edm::LuminosityBlockRange > jsonVector;

  bool jsonContainsEvent (const std::vector<edm::LuminosityBlockRange>&jsonVec,
                          const edm::EventBase &event);

  std::vector<edm::EventID> vEventsToDebug_;

  void setEventWeight(edm::EventBase const & event);
};

#endif//#define _AnalyzerBase_h_

