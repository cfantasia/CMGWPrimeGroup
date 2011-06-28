#ifndef _Ele_MET_Analyzer_h_
#define _Ele_MET_Analyzer_h_

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "UserCode/CMGWPrimeGroup/interface/WPrimeUtil.h"
#include "UserCode/CMGWPrimeGroup/interface/elmet_histo_constants.h"

#include "SHarper/HEEPAnalyzer/interface/HEEPEleSelector.h"
#include "SHarper/HEEPAnalyzer/interface/HEEPEle.h"

#include "TLorentzVector.h"

class TH1F;

#define debugmeElMet 0
#define dumpHighPtElectrons 0

class EleMETAnalyzer;
// function signature: flag indicating whether particular electron satisfies the 
// given selection cut and needs to be histogrammed; returns false if rest 
// of selection cuts should be skipped based on some event property 
// (e.g. when the trigger has failed the event, or there are more than 
// one electrons in the event, etc)
typedef bool (EleMETAnalyzer::*funcPtrEl)(bool *, const heep::Ele &, edm::EventBase const &);

// key: cuts_desc_short[i], value: function pointer corresponding to selection cut
typedef std::map<std::string, funcPtrEl> selection_map_elmet;

class EleMETAnalyzer
{
 public:
  explicit EleMETAnalyzer(const edm::ParameterSet& cfg, 
			 WPrimeUtil * wprimeUtil);
  ~EleMETAnalyzer();

  void eventLoop(edm::EventBase const & event);
  // operations to be done when changing input file (e.g. create new histograms)
  void beginFile(std::vector<wprime::InputFile>::const_iterator file);
  // operations to be done when closing input file 
  // (e.g. print summary)
  void endFile(std::vector<wprime::InputFile>::const_iterator it,
	       ofstream & out);
  // e.g. print summmary of expected events for all samples
  void endAnalysis(ofstream & out);

 private:
  WPrimeUtil * wprimeUtil_;

  heep::EleSelector cuts_; //allows us to apply the heep selection

  edm::InputTag electrons_;
  edm::InputTag met_;
  // Handle to the electron collection
  edm::Handle<pat::ElectronCollection > electrons;
  // Handle to the (pf)MET collection
  edm::Handle<pat::METCollection > met;
   // keeps track of selection efficiencies for all input samples & cuts
  wprime::SampleStat stats;

  bool highestEtElectronOnly_; // whether to only consider highest-pt electron in event
  bool dumpHighEtElectrons_; // whether to dump high-pt electrons for data
  float dumpHighEtElectronThreshold_;

  void defineHistos(TFileDirectory & dir);
  void defineHistos_ElectronEt(TFileDirectory & dir);
  void defineHistos_ElectronEta(TFileDirectory & dir);
  void defineHistos_ElectronPhi(TFileDirectory & dir);
  void defineHistos_TMass(TFileDirectory & dir);

  void setupCutOrder();
  selection_map_elmet cuts;

  // Get the hardest electron (based on HEEP Et) in event
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
  // run HEEP cuts, return HEEPEleSelector::getCutCode
  int runHEEPcuts(const heep::Ele & el);

  static int ignoreIsolationMask;
  static int useOnlyIsolationMask;

  // set electron 4-d momentum (sets el4D)
  void setElectronMomentum(const heep::Ele & el);

  // Get new MET: here (unlike the muon case) there is only one correction to be made:
  // the hadronic MET component (that needs to be corrected 
  // if applyCorrection=true) from Z data; this will be done according to hadronic 
  // activity from Z->mumu (not ee?) reconstructed events
  TVector2 getNewMET(edm::EventBase const & event, const TLorentzVector & el_p);

  // whether HLT accepted the event
  bool passedHLT(bool *, const heep::Ele &, edm::EventBase const &);

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
  
  float electronPtThreshold_;
  float oneEleEtCut_;

  TH1F * hPT[Num_elmet_cuts];
  TH1F * hETA[Num_elmet_cuts];
  TH1F * hPHI[Num_elmet_cuts];
  TH1F * hTM[Num_elmet_cuts];

};


#endif //#define _Ele_MET_Analyzer_h_
