#ifndef _Ele_MET_Analyzer_h_
#define _Ele_MET_Analyzer_h_

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "UserCode/CMGWPrimeGroup/interface/AnalyzerBase.h"

#include "UserCode/CMGWPrimeGroup/interface/WPrimeUtil.h"
#include "UserCode/CMGWPrimeGroup/interface/elmet_histo_constants.h"

#include "UserCode/CMGWPrimeGroup/interface/WprimeTreeVariables.h"

// for Trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "SHarper/HEEPAnalyzer/interface/HEEPEleSelector.h"
#include "SHarper/HEEPAnalyzer/interface/HEEPEle.h"

#include "TLorentzVector.h"

class TH1F;

#define debugmeElMet 0

class EleMETAnalyzer;
// function signature: flag indicating whether particular electron satisfies the 
// given selection cut and needs to be histogrammed; returns false if rest 
// of selection cuts should be skipped based on some event property 
// (e.g. when the trigger has failed the event, or there are more than 
// one electrons in the event, etc)
typedef bool (EleMETAnalyzer::*funcPtrEl)(bool *, const heep::Ele &, edm::EventBase const &);

// key: cuts_desc_short[i], value: function pointer corresponding to selection cut
typedef std::map<std::string, funcPtrEl> selection_map_elmet;

class EleMETAnalyzer : public AnalyzerBase 
{
 public:
  explicit EleMETAnalyzer(const edm::ParameterSet& cfg, 
			 int fileToRun);
  ~EleMETAnalyzer();

  void eventLoop(edm::EventBase const & event);

 private:
  heep::EleSelector cuts_; //allows us to apply the heep selection

  // Handle to the electron collection
  edm::Handle<pat::ElectronCollection > electrons;
  // Handle to the (pf)MET collection
  edm::Handle<pat::METCollection > defMet;
  pat::MET met;

  // Trigger
  edm::InputTag triggerResults_;
  std::vector<std::string> HLTPathsByName_;
  std::vector<unsigned int> HLTPathsByIndex_;

  bool SameTrigger(std::string & A, std::string & B);

  bool highestEtElectronOnly_; // whether to only consider highest-pt electron in event
  bool dumpHighEtElectrons_; // whether to dump high-pt electrons for data
  float dumpHighEtElectronThreshold_;
  float dumpHighMtElectronThreshold_;

  bool mkTuple_;

  bool doEoP_;
  float electronEoverPthreshold_;

  void defineHistos(const TFileDirectory & dir);
  void defineTrees(const TFileDirectory & dir);
  void defineHistos_ElectronEt(const TFileDirectory & dir);
  void defineHistos_ElectronEta(const TFileDirectory & dir);
  void defineHistos_ElectronPhi(const TFileDirectory & dir);
  void defineHistos_TMass(const TFileDirectory & dir);
  void defineResolutionHistos(const TFileDirectory & dir, float Mass);

  void setupCutOrder();
  selection_map_elmet cuts;

  // get the hardest electron (based on HEEP Et) in event
  // (returns index in pat::ElectronCollection)
  int getTheHardestElectron();

  // fill histograms for electron if fill_entry=true; update book-keeping 
  // (via private member: stats); make sure stats gets updated maximum 
  // once per event
  void tabulateMe(int cut_index, bool accountMe[], 
		  edm::EventBase const & event, int theEle);
  
  // dump on screen info about high-Et electron
  void printHighEtElectron(edm::EventBase const & event);

  TLorentzVector el4D;
  // this is the return value of HEEPEleSelector::getCutCode; it is needed by
  // goodQualityElectron and isolatedElectron; cache here
  // so that we don't have to call method twice per electron
  int cutCode; 
  // make sure HEEP cuts are calculated once per electron
  void runHEEPcuts(const heep::Ele & el);

  static int ignoreIsolationMask;
  static int useOnlyIsolationMask;

  // set electron 4-d momentum (sets el4D)
  void setElectronMomentum(const heep::Ele & el);

  // whether HLT accepted the event
  bool passedHLT(bool *, const heep::Ele &, edm::EventBase const & event);

  // check if electron satisfies quality requirements
  // fill goodQual; always returns true
  bool goodQualityElectron(bool * goodQual, const heep::Ele & el, edm::EventBase const &);

  // true if only one electron with track pt > the threshold
  bool onlyOneHighEtElectron(bool *, const heep::Ele &, edm::EventBase const &);

  // returns # of (global) electrons with Et above <Et_thresh>
  unsigned nEleAboveThresh(float Et_thresh);

  // set bool flag to true if electron isolated
  // always returns true
  bool isolatedElectron(bool * goodQual, const heep::Ele & ele, edm::EventBase const &);

  // check if electron, MET pass kinematic cuts, updated goodQual
  // always returns true
  bool kinematicCuts(bool * goodQual, const heep::Ele &, edm::EventBase const & event);

  // print summary of efficiencies
  void printFileSummary(std::vector<wprime::InputFile>::const_iterator,
			ofstream & out);
  
  ///These are required functions when inheriting from AnalyzerBase
  void fillHistos(const int& index, const float& weight=1.);

  void FillNtuple(int & theEvent, int & iEleMin, edm::EventBase const & event);

  WCandidate Wcand;

  float electronPtThreshold_;
  float oneEleEtCut_;
 
  TH1F * hPT[Num_elmet_cuts];
  TH1F * hETA[Num_elmet_cuts];
  TH1F * hPHI[Num_elmet_cuts];
  TH1F * hTM[Num_elmet_cuts];

  // one El+MET resolution + one Mt histograms for each signal sample
  std::vector<TH1F*> elMetRes; 
  std::vector<TH1F*> elMetGenMt;


  std::string analysis;
  WprimeVariables vars;
  TTree * cloneTrees[Num_elmet_cuts];
  TTree * NTuple;

};


#endif //#define _Ele_MET_Analyzer_h_
