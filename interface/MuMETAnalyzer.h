#ifndef _Mu_MET_Analyzer_h_
#define _Mu_MET_Analyzer_h_

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "UserCode/CMGWPrimeGroup/interface/AnalyzerBase.h"

#include "UserCode/CMGWPrimeGroup/interface/WPrimeUtil.h"
#include "UserCode/CMGWPrimeGroup/interface/mumet_histo_constants.h"
#include "UserCode/CMGWPrimeGroup/interface/TeVMuon.h"
#include "UserCode/CMGWPrimeGroup/interface/BosonFinder.h"

#include "TLorentzVector.h"

class TH1F;
class TH2F;
class TeVMuon;

#define debugmeMuMet 0

class MuMETAnalyzer;
// function signature: flag indicating whether particular muon satisfies the 
// given selection cut and needs to be histogrammed; returns false if rest 
// of selection cuts should be skipped based on some event property 
// (e.g. when the trigger has failed the event, or there are more than 
// one muons in the event, etc)
typedef bool (MuMETAnalyzer::*funcPtrMu)(bool *, const TeVMuon *, edm::EventBase const &);

// key: cuts_desc_short[i], value: function pointer corresponding to selection cut
typedef std::map<std::string, funcPtrMu> selection_map_mumet;

class MuMETAnalyzer : public AnalyzerBase 
{
 public:
  explicit MuMETAnalyzer(const edm::ParameterSet& cfg, 
			 WPrimeUtil * wprimeUtil);
  ~MuMETAnalyzer();

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
  // Handle to the muon collection
  edm::Handle<pat::MuonCollection > muons;
  // Handle to the (pf)MET collection
  edm::Handle<pat::METCollection > defMet;
  pat::MET met;
   // keeps track of selection efficiencies for all input samples & cuts
  wprime::SampleStat stats;

  MuonV vmuons;

  // true if TrackRef for chosen high-pt muon reconstructor is null;
  // to be reset at beginning of loop-over-muons
  bool isInvalidMuon_;

  // identifies muon reconstructor (see TeVMuon.h)
  bool highestPtMuonOnly_; // whether to only consider highest-pt muon in event
  bool dumpHighPtMuons_; // whether to dump high-pt muons for data
  float dumpHighPtMuonThreshold_;
  float dumpHighMtMuonThreshold_;

  void defineHistos(const TFileDirectory & dir);
  void defineHistos_MuonPt(const TFileDirectory & dir);
  void defineHistos_MuonEta(const TFileDirectory & dir);
  void defineHistos_MuonPhi(const TFileDirectory & dir);
  //  void defineHistos_MuonJetDPhi(const TFileDirectory & dir);
  void defineHistos_MuonIso(const TFileDirectory & dir);
  void defineHistos_TMass(const TFileDirectory & dir);
  void defineHistos_TMvPT(const TFileDirectory & dir);

  void setupCutOrder();
  selection_map_mumet cuts;

  // get the hardest muon (based on tracker-pt) in event
  // (returns index in pat::MuonCollection)
  int getTheHardestMuon();

  // fill histograms for muon if fill_entry=true; update book-keeping 
  // (via private member: stats); make sure stats gets updated maximum 
  // once per event
  void tabulateMe(int cut_index, bool accountMe[], 
		  edm::EventBase const & event, const TeVMuon * muon);
  
  // dump on screen info about high-pt muon
  void printHighPtMuon(edm::EventBase const & event, const pat::Muon & muon);

  TLorentzVector mu4D;

  // whether HLT accepted the event
  bool passedHLT(bool *, const TeVMuon * muon, edm::EventBase const &);

  // check if muon satisfies quality requirements
  // fill goodQual; always returns true
  bool goodQualityMuon(bool * goodQual, const TeVMuon * muon, edm::EventBase const &);

  // true if only one muon with track pt > the threshold
  bool onlyOneHighTrackPtMuon(bool *, const TeVMuon *, edm::EventBase const &);

  // returns # of (global) muons with tracker-pt above <tracker_muon_pt>
  unsigned nMuAboveThresh(float tracker_muon_pt);

  // set bool flag to true if muon isolated
  // always returns true
  bool isolatedMuon(bool * goodQual, const TeVMuon * muon, edm::EventBase const &);

  // check if muon, MET pass kinematic cuts, updated goodQual
  // always returns true
  bool kinematicCuts(bool * goodQual, const TeVMuon * muon, edm::EventBase const & event);

  // print summary of efficiencies
  void printFileSummary(std::vector<wprime::InputFile>::const_iterator,
			ofstream & out);

  ///These are required functions when inheriting from AnalyzerBase
  void fillHistos(const int& index, const float& weight=1.);

  WCandidate Wcand;
  
  float muonPtThreshold_;
  float chi2Cut_;
  float muonEtaCut_;
  float oneMuPtTrackCut_;
  float relIsoCut_;

  TH1F * hPTgen[Num_mumet_cuts];
  TH1F * hPT[Num_mumet_cuts];
  TH1F * hETA[Num_mumet_cuts];
  TH1F * hPHI[Num_mumet_cuts];
  //  TH1F * hMJDPHI[Num_mumet_cuts];
  TH1F * hISO[Num_mumet_cuts];
  TH1F * hTM[Num_mumet_cuts];
  TH2F * hTMvPT[Num_mumet_cuts];

  std::vector<unsigned> reconstructors;


};


#endif //#define _Mu_MET_Analyzer_h_
