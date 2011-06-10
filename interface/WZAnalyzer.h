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
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "UserCode/CMGWPrimeGroup/interface/WPrimeUtil.h"
#include "UserCode/CMGWPrimeGroup/interface/WZUtilities.h"
#include "UserCode/CMGWPrimeGroup/interface/BosonFinder.h"
#include "UserCode/CMGWPrimeGroup/interface/TeVMuon.h"

class WZAnalyzer  {
public:
  WZAnalyzer();                         // constructor; initialize the list to be empty
  WZAnalyzer(const edm::ParameterSet & cfg, WPrimeUtil * wprimeUtil);
  ~WZAnalyzer();

  void FillCutFns();
  void getEff(float & eff, float & deff, float Num, float Denom);

  void ResetCounters();
  void Declare_Histos(TFileDirectory& dir);
  void DeclareHistoSet(std::string n, std::string t, std::string xtitle,
                       int nbins, float min, float max,
                       std::vector<TH1F*>& h, TFileDirectory& d);
  void DeclareHisto(std::string n, std::string t, std::string xtitle,
                    int nbins, float min, float max,
                    TH1F* h, TFileDirectory& d);

  void ClearEvtVariables();
  void ClearAndResize(vector<TH1F*>& h, int& size, TH1F* ptr=NULL);

  void ScaleHistos();
  void Fill_Histos(int index, float weight=1.);
  void saveHistos(std::string dir);
  void deleteHistos();
  void printSummary(const std::string& dir, ofstream & out);

  void eventLoop(edm::EventBase const & event);
  void Tabulate_Me(int& cut_index,const float& weight);
  void CalcZVariables();
  void CalcWVariables();
  void CalcWZVariables();
  void CalcEventVariables();
  bool PassCuts(const float& weight=1.);

  void PrintEventToFile(edm::EventBase const & event);
  void PrintEvent(edm::EventBase const & event);
  void PrintEventFull(edm::EventBase const & event);
  void PrintTrigger();
  void PrintElectron(const heep::Ele* elec, int parent);
  void PrintMuon(const TeVMuon* mu, int parent);
  double CalcLeadPt(int type=0);
  double CalcQ();
  
  bool SameTrigger(string & A, string & B);

//methods for utilities

//methods for modifiers
  void SetCandEvtFile(std::string s);

//methods for the cuts
  bool PassNoCut();
  bool PassTriggersCut();
  bool PassNLeptonsCut();
  bool PassValidWCut();
  bool PassValidZCut();
  bool PassValidWandZCut();
  bool PassValidWZCandCut();
  bool PassLeadingLeptonPtCut();
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

  bool PassElecLooseCut(const heep::Ele* elec);
  bool PassElecTightCut(const heep::Ele* elec);
  bool PassElecLooseEtCut(const heep::Ele* elec);
  bool PassElecTightEtCut(const heep::Ele* elec);
  bool PassElecLooseWPCut(const heep::Ele* elec);
  bool PassElecWPRelIsoCut(const heep::Ele* elec);

  bool PassMuonLooseCut(const TeVMuon* mu);
  bool PassMuonTightCut(const TeVMuon* mu);
  bool PassMuonLoosePtCut(const TeVMuon* mu);
  bool PassMuonTightPtCut(const TeVMuon* mu);
  bool PassMuonGlobalCut(const TeVMuon* mu);
  bool PassMuonDxyCut(const TeVMuon* mu);
  bool PassMuonNpixhitCut(const TeVMuon* mu);
  bool PassMuonNtrkhitCut(const TeVMuon* mu);
  bool PassMuonNormChi2Cut(const TeVMuon* mu);
  bool PassMuonHitsUsedCut(const TeVMuon* mu);
  bool PassMuonStationsCut(const TeVMuon* mu);
  bool PassMuonEtaCut(const TeVMuon* mu);
  bool PassMuonCombRelIsoCut(const TeVMuon* mu);

  int   Calc_EvtType();
  float Calc_Q();
  float Calc_Ht();
  float CalcElecSc(const heep::Ele* elec);
  float Calc_MuonRelIso(const TeVMuon* mu);
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

  std::string electronsLabel_;
  std::string muonsLabel_;
  std::string metLabel_;
  std::string hltEventLabel_;
  std::string pileupLabel_;
  vstring triggersToUse_;

  int PDGMUON;
  int PDGELEC;
  int PDGW;
  int PDGZ;
  int PDGWPRIME;

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
  
// +++++++++++++++++++location of data files and samples info
  WPrimeUtil * wprimeUtil_;
  ofstream outCandEvt_;

  uint muonAlgo_;

///My calculated qualities//////////////////
  float Ht_;
  float Q_;
  uint evtType_;
  uint numZs_;
  float LeadPt_;
  float LeadElecPt_;
  float LeadMuonPt_;
  bool TT, TF;

// +++++++++++++++++++General Cut values
  uint maxNumZs_;
  uint minNLeptons_;
  float minLeadPt_;
  float minMET_;

// +++++++++++++++++++Ht Cuts
  float minHt_;

// +++++++++++++++++++W Cuts
  float minWtransMass_;
  float minWpt_;

  float maxWmunuCombRelIso_;
  int   cutWenuWPRelIsoMask_;
  string cutElecWPTightType_;
  float minDeltaR_;

// +++++++++++++++++++Z Cuts
  float minZpt_;
  float minZmass_;
  float maxZmass_;

// +++++++++++++++++++Electron General Cuts
//VBTF Recommended Cuts
  float minElecLooseEt_;
  float minElecTightEt_;
  int cutElecWPLooseMask_;
  string cutElecWPLooseType_;
  std::vector<double> maxElecSigmaiEtaiEta_;
  std::vector<double> maxElecDeltaPhiIn_;
  std::vector<double> maxElecDeltaEtaIn_;
  std::vector<double> maxElecHOverE_    ;

// +++++++++++++++++++Muon General Cuts
  float maxMuonEta_;
  float minMuonLoosePt_;
  float minMuonTightPt_;
//VBTF Recommended Cuts
  float maxMuonDxy_;
  float maxMuonNormChi2_;
  int minMuonNPixHit_;
  int minMuonNTrkHit_;
  int minMuonStations_;
  int minMuonHitsUsed_;

//////Chosen Candidates
  ElectronV electrons_, looseElectrons_, tightElectrons_;
  MuonV muons_, looseMuons_, tightMuons_;
  pat::MET met_;
  ZCandidate zCand_;
  WCandidate wCand_;
  WZCandidate wzCand_;
  pat::TriggerEvent triggerEvent_; 
  std::vector< PileupSummaryInfo > PupInfo_; 
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
  typedef bool (WZAnalyzer::*ElecCutFnPtr)(const heep::Ele*); 
  typedef bool (WZAnalyzer::*MuonCutFnPtr)(const TeVMuon*); 
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
