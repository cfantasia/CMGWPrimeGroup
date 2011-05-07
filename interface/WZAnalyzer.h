#ifndef _WZAnalyzer_h_
#define _WZAnalyzer_h_

#include <vector>
#include <string>

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "UserCode/CMGWPrimeGroup/interface/WPrimeUtil.h"
#include "UserCode/CMGWPrimeGroup/interface/WZUtilities.h"

class WZAnalyzer  {
public:
  WZAnalyzer();                         // constructor; initialize the list to be empty
  WZAnalyzer(const edm::ParameterSet & cfg, WPrimeUtil * wprimeUtil);
  ~WZAnalyzer();

  void FillCutFns();
  void getEff(float & eff, float & deff, float Num, float Denom);
  double deltaEta(double eta1, double eta2);
  double deltaPhi(double phi1, double phi2);

  void ResetCounters();
  void Declare_Histos(TFileDirectory& dir);
  void DeclareHistoSet(std::string n, std::string t, std::string xtitle,
                       int nbins, float min, float max,
                       std::vector<TH1F*>& h, TFileDirectory& d);
  void DeclareHisto(std::string n, std::string t, std::string xtitle,
                    int nbins, float min, float max,
                    TH1F* h, TFileDirectory& d);

  double deltaR(double eta1, double phi1, double eta2, double phi2);

  void ScaleHistos();
  void Fill_Histos(int index, float weight=1.);
  void saveHistos(std::string dir);
  void deleteHistos();
  void printSummary(const std::string& dir);

  void eventLoop(edm::EventBase const & event);
  void Tabulate_Me(int& cut_index,const float& weight);
  void UseSample(std::string dir);
  void CalcEventVariables();
  bool PassCuts(const float& weight=1.);

  void PrintEventToFile(edm::EventBase const & event);
  void PrintEvent(edm::EventBase const & event);
  void PrintEventFull(edm::EventBase const & event);
  void PrintElectron(const pat::Electron* elec, int parent);
  void PrintMuon(const pat::Muon* mu, int parent);
  double CalcLeadPt(int type=0);
  double CalcQ();

//methods for utilities
  void CheckStream(ofstream& stream, std::string s);

//methods for modifiers
  void SetCandEvtFile(std::string s);
  void SetLogFile(std::string s);
  void SetOutputFile(std::string s);

//methods for the cuts
  bool PassNoCut();
  bool PassTriggersCut();
  bool PassValidWCut();
  bool PassValidZCut();
  bool PassValidWandZCut();
  bool PassValidWZCandCut();
  bool PassNumberOfZsCut();
  bool PassWptCut();
  bool PassZptCut();
  bool PassHtCut();
  bool PassHtMetCut();
  bool PassMETCut();

  bool PassZMassCut();
  bool PassZLepPtCut();
  bool PassZeePtCut();
  bool PassZmumuPtCut();

  bool PassWtransMassCut();
  bool PassWLepPtCut();
  bool PassWLepIsoCut();
  bool PassWenuIsoCut();
  bool PassWmunuIsoCut();

  bool PassElecLooseCut(const pat::Electron* elec);
  bool PassElecTightCut(const pat::Electron* elec);
  bool PassElecLooseEtCut(const pat::Electron* elec);
  bool PassElecTightEtCut(const pat::Electron* elec);
  bool PassElecLooseWPCut(const pat::Electron* elec);
  bool PassElecWPRelIsoCut(const pat::Electron* elec);

  bool PassMuonLooseCut(const pat::Muon* mu);
  bool PassMuonTightCut(const pat::Muon* mu);
  bool PassMuonLoosePtCut(const pat::Muon* mu);
  bool PassMuonTightPtCut(const pat::Muon* mu);
  bool PassMuonGlobalCut(const pat::Muon* mu);
  bool PassMuonDxyCut(const pat::Muon* mu);
  bool PassMuonNpixhitCut(const pat::Muon* mu);
  bool PassMuonNtrkhitCut(const pat::Muon* mu);
  bool PassMuonNormChi2Cut(const pat::Muon* mu);
  bool PassMuonHitsUsedCut(const pat::Muon* mu);
  bool PassMuonStationsCut(const pat::Muon* mu);
  bool PassMuonEtaCut(const pat::Muon* mu);
  bool PassMuonCombRelIsoCut(const pat::Muon* mu);

  int   Calc_EvtType();
  float Calc_Q();
  float Calc_Ht();
  float CalcElecSc(const pat::Electron* elec);
  float Calc_MuonNormChi2(const pat::Muon* mu);
  float Calc_MuonRelIso(const pat::Muon* mu);
  float Calc_GenWZInvMass();

//////////
  void reportProgress(int eventNum);
  void verbose(const char *string, ...);

  void beginFile(std::vector<wprime::InputFile>::const_iterator fi);
  void endFile(std::vector<wprime::InputFile>::const_iterator fi, ofstream & out);
  void endAnalysis(ofstream & out);

// +++++++++++++++++++useful constants
  bool debugme;//print stuff if active
  bool doPreselect_;

  int PDGMUON;
  int PDGELEC;
  int PDGW;
  int PDGZ;
  int PDGWPRIME;

  float PDGZMASS;
  float W_mass;

  double PI;
  double TWOPI;
  float NOCUT;

// +++++++++++++++++++
  std::vector<uint> nEvents;
  int eventNum;
  uint runNumber;
  uint lumiID;

  typedef pair< map<string, bool>, vector<string> > OptArgPair;
  std::map<std::string, int> intOptions_;
  std::map<std::string, std::string> stringOptions_;
  std::vector<std::string> filenames_;
  std::string datasetName;
  
// +++++++++++++++++++location of data files and samples info
  std::string top_level_dir;
  ofstream outCandEvt;
  ofstream outLogFile;

  WPrimeUtil * wprimeUtil_;

///My calculated qualities//////////////////
  float Ht;
  float Q;
  int   evtType;
  int numZs;
  float LeadPt;
  float LeadElecPt;
  float LeadMuonPt;
  bool TT, TF;

// +++++++++++++++++++General Cut values
  int maxNumZs;
  int minNumLeptons;
  float minMET;

// +++++++++++++++++++Ht Cuts
  float minHt;

// +++++++++++++++++++W Cuts
  float minWtransMass;
  float minWpt;

  float maxWmunuCombRelIso;
  int   cutWenuWPRelIsoMask;

// +++++++++++++++++++Z Cuts
  float minZpt;
  float minZmass;
  float maxZmass;

// +++++++++++++++++++Electron General Cuts
//VBTF Recommended Cuts
  float minElecLooseEt;
  float minElecTightEt;
  int cutElecWPLooseMask;
  std::vector<double> maxElecSigmaiEtaiEta;
  std::vector<double> maxElecDeltaPhiIn;
  std::vector<double> maxElecDeltaEtaIn;
  std::vector<double> maxElecHOverE    ;

// +++++++++++++++++++Muon General Cuts
  float maxMuonEta;
  float minMuonLoosePt;
  float minMuonTightPt;
//VBTF Recommended Cuts
  float maxMuonDxy;
  float maxMuonNormChi2;
  int minMuonNPixHit;
  int minMuonNTrkHit;
  int minMuonStations;
  int minMuonHitsUsed;

//////Chosen Candidates
  pat::MET met;
  ZCandidate zCand;
  WCandidate wCand;
  WZCandidate wzCand;
  
  std::vector<std::string>* TriggerHLTNames;
  std::vector<int>* TriggerHLTPrescales;
  std::vector<bool>* TriggerHLTDecisions;

  std::vector<float> Num_surv_cut_;

// +++++++++++++++++++ Histogram Definitions
  std::vector<TH1F*> listOfHists;
  TH1F * hEffRel;
  TH1F * hEffAbs;
  TH1F * hNumEvts;

  std::vector<TH1F*> hEvtType;

  std::vector<TH1F*> hWZInvMass     ;
  std::vector<TH1F*> hWZ3e0muInvMass;
  std::vector<TH1F*> hWZ2e1muInvMass;
  std::vector<TH1F*> hWZ1e2muInvMass;
  std::vector<TH1F*> hWZ0e3muInvMass;

  std::vector<TH1F*> hWZTransMass;
  std::vector<TH1F*> hHt;
  std::vector<TH1F*> hWpt;
  std::vector<TH1F*> hZpt;
  std::vector<TH1F*> hMET;

  std::vector<TH1F*> hZMass     ;
  std::vector<TH1F*> hZeeMass   ;
  std::vector<TH1F*> hZmumuMass ;
  std::vector<TH1F*> hZ3e0muMass;
  std::vector<TH1F*> hZ2e1muMass;
  std::vector<TH1F*> hZ1e2muMass;
  std::vector<TH1F*> hZ0e3muMass;
  std::vector<TH1F*> hZeeMassTT;
  std::vector<TH1F*> hZeeMassTF;
  std::vector<TH1F*> hZmumuMassTT ;
  std::vector<TH1F*> hZmumuMassTF ;

  std::vector<TH1F*> hWTransMass     ;
  std::vector<TH1F*> hWenuTransMass  ;
  std::vector<TH1F*> hWmunuTransMass ;
  std::vector<TH1F*> hW3e0muTransMass;
  std::vector<TH1F*> hW2e1muTransMass;
  std::vector<TH1F*> hW1e2muTransMass;
  std::vector<TH1F*> hW0e3muTransMass;

  std::vector<TH1F*> hQ;

  std::vector<TH1F*> hLeadPt;
  std::vector<TH1F*> hLeadElecPt;
  std::vector<TH1F*> hLeadMuonPt;

  std::vector<TH1F*> hElecPt;
  std::vector<TH1F*> hElecEt;
  std::vector<TH1F*> hElecdEta;
  std::vector<TH1F*> hElecdPhi;
  std::vector<TH1F*> hElecSigmann;
  std::vector<TH1F*> hElecEP;
  std::vector<TH1F*> hElecHE;
  std::vector<TH1F*> hElecTrkRelIso;
  std::vector<TH1F*> hElecECalRelIso;
  std::vector<TH1F*> hElecHCalRelIso;

  std::vector<TH1F*> hMuonPt;
  std::vector<TH1F*> hMuonDxy;
  std::vector<TH1F*> hMuonNormChi2;
  std::vector<TH1F*> hMuonNPix;
  std::vector<TH1F*> hMuonNTrk;
  std::vector<TH1F*> hMuonRelIso;
  std::vector<TH1F*> hMuonStation;
  std::vector<TH1F*> hMuonSip;

//Cuts
  
  int NCuts_;
  std::vector<std::string> Cuts_;
  int NLooseElecCuts_;
  std::vector<std::string> LooseElecCuts_;
  int NLooseMuonCuts_;
  std::vector<std::string> LooseMuonCuts_;
  int NTightElecCuts_;
  std::vector<std::string> TightElecCuts_;
  int NTightMuonCuts_;
  std::vector<std::string> TightMuonCuts_;
  typedef bool (WZAnalyzer::*    CutFnPtr)(); 
  typedef bool (WZAnalyzer::*ElecCutFnPtr)(const pat::Electron*); 
  typedef bool (WZAnalyzer::*MuonCutFnPtr)(const pat::Muon*); 
#ifndef __CINT__
  std::map<std::string,     CutFnPtr> mFnPtrs_;
  std::map<std::string, ElecCutFnPtr> mElecFnPtrs_;
  std::map<std::string, MuonCutFnPtr> mMuonFnPtrs_;
  std::vector<    CutFnPtr> CutFns_;
  std::vector<ElecCutFnPtr> LooseElecCutFns_;
  std::vector<MuonCutFnPtr> LooseMuonCutFns_;
  std::vector<ElecCutFnPtr> TightElecCutFns_;
  std::vector<MuonCutFnPtr> TightMuonCutFns_;
#endif
};

/// Return a pointer to product with productName
template <class T, class P>
const T * getPointerFWLite(const P & ev, string productName)
{
  fwlite::Handle<T> handle;
  handle.getByLabel(ev, productName.c_str());
  if (handle.isValid()) {
    return handle.ptr();
  }
  return 0;
}

#endif//#define _WZAnalyzer_h_
