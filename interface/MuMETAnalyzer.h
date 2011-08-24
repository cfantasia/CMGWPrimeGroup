#ifndef _Mu_MET_Analyzer_h_
#define _Mu_MET_Analyzer_h_

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "UserCode/CMGWPrimeGroup/interface/WPrimeUtil.h"
#include "UserCode/CMGWPrimeGroup/interface/mumet_histo_constants.h"
#include "UserCode/CMGWPrimeGroup/interface/TeVMuon.h"

#include "TLorentzVector.h"

class TH1F;
class TeVMuon;

#define debugmeMuMet 0
#define dumpHighPtMuons 0

class MuMETAnalyzer;
// function signature: flag indicating whether particular muon satisfies the 
// given selection cut and needs to be histogrammed; returns false if rest 
// of selection cuts should be skipped based on some event property 
// (e.g. when the trigger has failed the event, or there are more than 
// one muons in the event, etc)
typedef bool (MuMETAnalyzer::*funcPtrMu)(bool *, const TeVMuon *, edm::EventBase const &);

// key: cuts_desc_short[i], value: function pointer corresponding to selection cut
typedef std::map<std::string, funcPtrMu> selection_map_mumet;

class MuMETAnalyzer
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
  WPrimeUtil * wprimeUtil_;

  edm::InputTag muonsLabel_;
  edm::InputTag metLabel_;
  edm::InputTag pfLabel_;
  // Handle to the muon collection
  edm::Handle<pat::MuonCollection > muons;
  // Handle to the (pf)MET collection
  edm::Handle<pat::METCollection > defMet;
  pat::MET met;
   // keeps track of selection efficiencies for all input samples & cuts
  wprime::SampleStat stats;

  MuonV vmuons;

  bool useAdjustedMET_;

  // true if TrackRef for chosen high-pt muon reconstructor is null;
  // to be reset at beginning of loop-over-muons
  bool isInvalidMuon_;

  // identifies muon reconstructor (see TeVMuon.h)
  unsigned muReconstructor_; 
  bool highestPtMuonOnly_; // whether to only consider highest-pt muon in event
  bool dumpHighPtMuons_; // whether to dump high-pt muons for data
  float dumpHighPtMuonThreshold_;

  void defineHistos(TFileDirectory & dir);
  void defineHistos_MuonPt(TFileDirectory & dir);
  void defineHistos_MuonEta(TFileDirectory & dir);
  void defineHistos_MuonPhi(TFileDirectory & dir);
  //  void defineHistos_MuonJetDPhi(TFileDirectory & dir);
  void defineHistos_MuonIso(TFileDirectory & dir);
  void defineHistos_TMass(TFileDirectory & dir);

  void setupCutOrder();
  selection_map_mumet cuts;

  // Get the hardest muon (based on tracker-pt) in event
  // (returns index in pat::MuonCollection)
  int getTheHardestMuon();

  // fill histograms for muon if fill_entry=true; update book-keeping 
  // (via private member: stats); make sure stats gets updated maximum 
  // once per event
  void tabulateMe(int cut_index, bool accountMe[], 
		  edm::EventBase const & event, const TeVMuon * muon);
  
  // dump on screen info about high-pt muon
  void printHighPtMuon(edm::EventBase const & event, TeVMuon & muon);

  TLorentzVector mu4D;

  // Get new MET: there are two corrections to be made:
  // (a) the hadronic MET component (that needs to be corrected 
  // if applyCorrection=true) from Z data; this will be done according to hadronic 
  // activity from Z->mumu reconstructed events
  // (b) the muon-pt component that needs to be updated if we switch to one
  // of the dedicated high-pt muon reconstructors
  TVector2 getNewMET(edm::EventBase const & event, const TLorentzVector & mu_p);

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
  bool kinematicCuts(bool * goodQual, const TeVMuon *, edm::EventBase const & event);

  // print summary of efficiencies
  void printFileSummary(std::vector<wprime::InputFile>::const_iterator,
			ofstream & out);
  
  // get (PF) MET without the default-pt for the running muon in event (mu4D);
  // this is done so that we can adjust the muon-pt component of the MET by 
  // switching to one of the dedicated high-pt muon reconstructors
  TVector2 getPFMETwithoutMu(edm::EventBase const & event);
  bool pfMETwithoutMuCalculated_; // want to calculate this max. once for each muon
  TVector2 pfMETwithoutMuCached_; 


			
  float muonPtThreshold_;
  float chi2Cut_;
  float muonEtaCut_;
  float oneMuPtTrackCut_;
  float relIsoCut_;

  TH1F * hPT[Num_mumet_cuts];
  TH1F * hETA[Num_mumet_cuts];
  TH1F * hPHI[Num_mumet_cuts];
  //  TH1F * hMJDPHI[Num_mumet_cuts];
  TH1F * hISO[Num_mumet_cuts];
  TH1F * hTM[Num_mumet_cuts];

  std::vector<unsigned> reconstructors;


};


#endif //#define _Mu_MET_Analyzer_h_
