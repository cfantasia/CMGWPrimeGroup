#ifndef _W_gamma_Analyzer_h_
#define _W_gamma_Analyzer_h_

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "UserCode/CMGWPrimeGroup/interface/MuMETAnalyzer.h"
#include "UserCode/CMGWPrimeGroup/interface/WPrimeUtil.h"
#include "UserCode/CMGWPrimeGroup/interface/wgamma_histo_constants.h"

#include "TLorentzVector.h"
typedef math::XYZTLorentzVector LorentzVector;

#include <string>

class TH1F;


#define debugmepho 0
//#define dumpHighPtMuons 0

class WgammaAnalyzer;
// function signature: flag indicating whether particular muon satisfies the 
// given selection cut and needs to be histogrammed; returns false if rest 
// of selection cuts should be skipped based on some event property 
// (e.g. when the trigger has failed the event, or there are more than 
// one muons in the event, etc)
typedef bool (WgammaAnalyzer::*funcPtrPho)(bool *, int, edm::EventBase const &);

// key: cuts_desc_short[i], value: function pointer corresponding to selection cut
typedef std::map<std::string, funcPtrPho> selection_map_wgamma;

class WgammaAnalyzer
{
 public:
  explicit WgammaAnalyzer(const edm::ParameterSet& cfg, 
			 WPrimeUtil * wprimeUtil);
  ~WgammaAnalyzer();

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
  edm::InputTag muons_;
  edm::InputTag met_;
  edm::InputTag metLabel_;
  edm::InputTag particleFlow_;
  edm::InputTag photons_;


  class ParticleStruct {
  public:
    int   isValid, pdgId;
    float et, pt, eta, phi;
    ParticleStruct () : isValid(-999), pdgId(-999), et(-999.), 
			pt(-999.),     eta(-999.),  phi(-999.) {};
  };
  class MetStruct : public ParticleStruct {
  public:
    MetStruct () : ParticleStruct () {};
    MetStruct (pat::MET*);
  };
  MetStruct recMet;



  // Handle to the muon collection
  edm::Handle<pat::MuonCollection > muons;
  // Handle to the (pf)MET collection
  edm::Handle<pat::METCollection > met;
  // Handle to the pat::Photon collection
  edm::Handle<pat::PhotonCollection> photons;

   // keeps track of selection efficiencies for all input samples & cuts
  wprime::SampleStat stats;

  // true if TrackRef for chosen high-pt muon reconstructor is null;
  // to be reset at beginning of loop-over-muons
  bool isInvalidMuon_;
  bool isInvalidPhoton_;

  // identifies muon reconstructor (see mumet_histo_constants.h)
  unsigned muReconstructor_; 
  bool highestPtMuonOnly_; // whether to only consider highest-pt muon in event
  bool highestPtPhotonOnly_;
  bool dumpHighPtMuons_; // whether to dump high-pt muons for data
  bool dumpHighPtPhotons_;
  float dumpHighPtMuonThreshold_;
  float dumpHighPtPhotonThreshold_;

  float barJurECALIsoConst_;
  float barJurECALIsoSlope_;
  float barTowHCALIsoConst_;
  float barTowHCALIsoSlope_;
  float barMaxHOverE_;
  float barHConeTrkIsoConst_;
  float barHConeTrkIsoSlope_;
  float barMaxEtaWidth_; 
  float endJurECALIsoConst_;
  float endJurECALIsoSlope_;
  float endTowHCALIsoConst_;
  float endTowHCALIsoSlope_;
  float endMaxHOverE_;
  float endHConeTrkIsoConst_;
  float endHConeTrkIsoSlope_;
  float endMaxEtaWidth_; 
  float minPhotonPt_;
  float maxPhotonEta_;  
  bool applyTrackVeto_;


  void defineHistos(TFileDirectory & dir);
  void defineHistos_MuonPt(TFileDirectory & dir);
  void defineHistos_MuonEta(TFileDirectory & dir);
  void defineHistos_MuonPhi(TFileDirectory & dir);
  void defineHistos_MuonMETDPhi(TFileDirectory & dir);
  void defineHistos_MuonIso(TFileDirectory & dir);
  void defineHistos_WTMass(TFileDirectory & dir);
  void defineHistos_PhotonPt(TFileDirectory & dir);
  void defineHistos_PhotonEta(TFileDirectory & dir);
  void defineHistos_PhotonPhi(TFileDirectory & dir);
  void defineHistos_MWgamma(TFileDirectory & dir);

  void setupCutOrderMuons();
  void setupCutOrderPhotons();

  selection_map_wgamma cuts_mu;
  selection_map_wgamma cuts_pho;

  // Get the hardest muon (based on tracker-pt) in event
  // (returns index in pat::MuonCollection)
  int getTheHardestMuon();
  int getTheHardestPhoton();

  void setMuLorentzVector(TLorentzVector& P, const reco::TrackRef & trk);
  LorentzVector calculateNeutrinoP4(LorentzVector muonP4, TVector2 met);
 

  // fill histograms for muon if fill_entry=true; update book-keeping 
  // (via private member: stats); make sure stats gets updated maximum 
  // once per event
  void tabulateMu(int cut_index, bool accountMe[], 
		  edm::EventBase const & event, int theMu);
  void tabulatePho(int pho_cut_index, bool accountMe[Num_photon_cuts], 
                   edm::EventBase const & event, double & InvMass);


  
  // dump on screen info about high-pt muon
  void printHighPtMuon(edm::EventBase const & event);
  void printHighPtPhoton(edm::EventBase const & event, int thePho);

  TLorentzVector mu4D;
  // set muon 4-d momentum according to muonReconstructor_ value (sets mu4D)
  void setMuonMomentum(int theMu);
  TLorentzVector PhotonP4;
  void setPhotonMomentum(int thePhoton);

  // Get new MET: there are two corrections to be made:
  // (a) the hadronic MET component (that needs to be corrected 
  // if applyCorrection=true) from Z data; this will be done according to hadronic 
  // activity from Z->mumu reconstructed events
  // (b) the muon-pt component that needs to be updated if we switch to one
  // of the dedicated high-pt muon reconstructors
  TVector2 getNewMET(edm::EventBase const & event, const TLorentzVector & mu_p);

  //computes the combined rel isolation value
  float combRelIsolation(int theMu);

  // whether HLT accepted the event
  bool passedHLT(bool *, int, edm::EventBase const &);

  // check if muon has minimum pt, fill isThere accordingly
  // always returns true
  bool muonMinimumPt(bool * isThere, int, edm::EventBase const &);
    
  // check if muon satisfies quality requirements
  // fill goodQual; always returns true
  bool goodQualityMuon(bool * goodQual, int theMu, edm::EventBase const &);


  // true if only one muon with track pt > the threshold
  bool onlyOneHighTrackPtMuon(bool *, int, edm::EventBase const &);

  // returns # of (global) muons with tracker-pt above <tracker_muon_pt>
  unsigned nMuAboveThresh(float tracker_muon_pt);

  // set bool flag to true if muon isolated
  // always returns true
  bool isolatedMuon(bool * goodQual, int theMu, edm::EventBase const & event);

  // check if muon, MET pass kinematic cuts, updated goodQual
  // always returns true
  bool kinematicCuts(bool * goodQual, int, edm::EventBase const & event);

  // min Photon Pt, max Eta
  bool photonPt(bool * goodQual, int thePho, edm::EventBase const &);

  // Photon Isolation
  bool photonIsolation(bool * goodQual, int thePho, edm::EventBase const &);

  //Photon Hadronic over EM
  bool photonHOverE(bool * goodQual, int thePho, edm::EventBase const &);

  //Photon sigmaIetaIeta
  bool photonEtaWidth(bool * goodQual, int thePho, edm::EventBase const &);

  //Photon Track Veto
  bool photonTrackVeto(bool * goodQual, int thePho, edm::EventBase const &);


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
  float combRelCut_;

  TH1F * hPT[Num_mumet_cuts];
  TH1F * hETA[Num_mumet_cuts];
  TH1F * hPHI[Num_mumet_cuts];
  TH1F * hMUMETDPHI[Num_mumet_cuts];
  TH1F * hISO[Num_mumet_cuts];
  TH1F * hWTM[Num_mumet_cuts];


  TH1F * hPHOPT[Num_photon_cuts];
  TH1F * hPHOETA[Num_photon_cuts];
  TH1F * hPHOPHI[Num_photon_cuts];
  TH1F * hPHOISO[Num_photon_cuts];
  TH1F * hMWG[Num_photon_cuts];





};


#endif //#define _Mu_MET_Analyzer_h_
