#include "UserCode/CMGWPrimeGroup/interface/WZAnalyzer.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
using namespace std;
WZAnalyzer::WZAnalyzer(){}
WZAnalyzer::WZAnalyzer(const edm::ParameterSet & cfg, WPrimeUtil * wprimeUtil){

  wprimeUtil_ = wprimeUtil;
  assert(wprimeUtil_);

   Cuts_          = cfg.getParameter<vstring>("Cuts");
  NCuts_          = Cuts_.size();

  FillCutFns();

  eSelectorPset_ = cfg.getParameter<PSet>("electronSelectors");
  looseElectronType_ = cfg.getParameter<string>("LooseElectronType");
  tightElectronType_ = cfg.getParameter<string>("TightElectronType");
  looseElectron_ = ElectronSelector(eSelectorPset_, looseElectronType_);
  tightElectron_ = ElectronSelector(eSelectorPset_, tightElectronType_);
  electronResult_ = looseElectron_.getBitTemplate();
  if(debugme) cout<<"Using "<<looseElectronType_<<" for Z electrons and "
                  <<tightElectronType_<<" for W electrons\n";

  mSelectorPset_ = cfg.getParameter<PSet>("muonSelectors");
  looseMuonType_ = cfg.getParameter<string>("LooseMuonType");
  tightMuonType_ = cfg.getParameter<string>("TightMuonType");
  looseMuon_ = MuonSelector(mSelectorPset_, looseMuonType_);
  tightMuon_ = MuonSelector(mSelectorPset_, tightMuonType_);
  muonResult_ = looseMuon_.getBitTemplate();
  if(debugme) cout<<"Using "<<looseElectronType_<<" for Z muons and "
                  <<tightElectronType_<<" for W muons\n";

  doPreselect_ = cfg.getParameter<bool>("preselect");

  debugme = cfg.getParameter<bool>("debugme");

  electronsLabel_ = cfg.getParameter<edm::InputTag>("electrons");
  muonsLabel_ = cfg.getParameter<edm::InputTag>("muons");
  pfCandsLabel_ = cfg.getParameter<edm::InputTag>("particleFlow");
  metLabel_ = cfg.getParameter<edm::InputTag>("met");

  muonAlgo_ = cfg.getParameter<uint>("muonReconstructor");if(debugme) cout<<"Using muon algo "<<algo_desc_long[muonAlgo_]<<endl;
  useAdjustedMET_ = cfg.getParameter<bool>("useAdjustedMET");

  hltEventLabel_ = cfg.getParameter<edm::InputTag>("hltEventTag");
  pileupLabel_ = cfg.getParameter<edm::InputTag>("pileupTag");
  vertexLabel_ = cfg.getParameter<edm::InputTag>("vertexTag");
  jetsLabel_ = cfg.getParameter<edm::InputTag>("jets");
  
  triggersToUse_ = cfg.getParameter<vstring>("triggersToUse");

  maxZMassDiff_ = cfg.getParameter<double>("maxZMassDiff");
  minDeltaR_ = cfg.getParameter<double>("minDeltaR");

  SetCandEvtFile(cfg.getParameter<string  >("candEvtFile" ));
  results_.assign(NCuts_,wprime::FilterEff());

  if(debugme) printf("Using %i cuts\n",NCuts_);

  effectiveElecArea_ = cfg.getParameter<vector<double> >("effectiveElecArea");
  effectiveMuonArea_ = cfg.getParameter<vector<double> >("effectiveMuonArea");
 
// +++++++++++++++++++General Cut values
  maxNumZs_ = cfg.getParameter<uint>("maxNumZs");
  minNLeptons_ = cfg.getParameter<uint>("minNLeptons");
  maxNLeptons_ = cfg.getParameter<uint>("maxNLeptons");
  minLeadPt_ = cfg.getParameter<double>("minLeadPt");
  minMET_ = cfg.getParameter<double>("minMET");

// +++++++++++++++++++Ht Cuts
  minHt_ = cfg.getParameter<double>("minHt");

// +++++++++++++++++++W Cuts
  minWtransMass_ = cfg.getParameter<double>("minWtransMass");
  minWpt_ = cfg.getParameter<double>("minWpt");
  minWlepPt_ = cfg.getParameter<double>("minWlepPt");

// +++++++++++++++++++Z Cuts
  minZpt_ = cfg.getParameter<double>("minZpt");
  minZmass_ = cfg.getParameter<double>("minZmass");
  maxZmass_ = cfg.getParameter<double>("maxZmass");
  minZeePt1_ = cfg.getParameter<double>("minZeePt1");
  minZeePt2_ = cfg.getParameter<double>("minZeePt2");
  minZmmPt1_ = cfg.getParameter<double>("minZmmPt1");
  minZmmPt2_ = cfg.getParameter<double>("minZmmPt2");

  ClearEvtVariables();
}

WZAnalyzer::~WZAnalyzer(){
  outCandEvt_.close();
}

void WZAnalyzer::FillCutFns(){
  mFnPtrs_["NoCuts"] = &WZAnalyzer::PassNoCut;
  mFnPtrs_["HLT"] = &WZAnalyzer::PassTriggersCut;
  mFnPtrs_["MinNLeptons"] = &WZAnalyzer::PassMinNLeptonsCut;
  mFnPtrs_["MaxNLeptons"] = &WZAnalyzer::PassMaxNLeptonsCut;
  mFnPtrs_["ValidW"] = &WZAnalyzer::PassValidWCut;
  mFnPtrs_["ValidWElec"] = &WZAnalyzer::PassValidWElecCut;
  mFnPtrs_["ValidWMuon"] = &WZAnalyzer::PassValidWMuonCut;
  mFnPtrs_["ValidZ"] = &WZAnalyzer::PassValidZCut;
  mFnPtrs_["ValidWandZ"] = &WZAnalyzer::PassValidWandZCut;
  mFnPtrs_["ValidWZCand"] = &WZAnalyzer::PassValidWZCandCut;
  mFnPtrs_["LeadLepPt"] = &WZAnalyzer::PassLeadingLeptonPtCut;
  mFnPtrs_["NumZs"] = &WZAnalyzer::PassNumberOfZsCut;
  mFnPtrs_["ZMass"] = &WZAnalyzer::PassZMassCut;
  mFnPtrs_["WTransMass"] = &WZAnalyzer::PassWtransMassCut;
  mFnPtrs_["MET"] = &WZAnalyzer::PassMETCut;
  mFnPtrs_["Ht"] = &WZAnalyzer::PassHtCut;
  mFnPtrs_["Zpt"] = &WZAnalyzer::PassZptCut;
  mFnPtrs_["Wpt"] = &WZAnalyzer::PassWptCut;
  mFnPtrs_["ZLepPt"] = &WZAnalyzer::PassZLepPtCut;
  mFnPtrs_["WLepPt"] = &WZAnalyzer::PassWLepPtCut;
  mFnPtrs_["ZLepTrigMatch"] = &WZAnalyzer::PassZLepTriggerMatchCut;
//  mFnPtrs_["WLepIso"] = &WZAnalyzer::PassWLepIsoCut;
  mFnPtrs_["AllCuts"] = &WZAnalyzer::PassNoCut;

  mFnPtrs_["WLepTight"] = &WZAnalyzer::PassWLepTightCut;
  mFnPtrs_["WFlavorElec"] = &WZAnalyzer::PassWFlavorElecCut;
  mFnPtrs_["WFlavorMuon"] = &WZAnalyzer::PassWFlavorMuonCut;
  mFnPtrs_["FakeEvt"] = &WZAnalyzer::PassFakeEvtCut;
  mFnPtrs_["FakeLepTag"] = &WZAnalyzer::PassFakeLeptonTagCut;
  mFnPtrs_["FakeLepProbeLoose"] = &WZAnalyzer::PassFakeLeptonProbeLooseCut;
  mFnPtrs_["FakeLepProbeTight"] = &WZAnalyzer::PassFakeLeptonProbeTightCut;

  CutFns_.resize(NCuts_);
  for(int i=0; i<NCuts_; ++i){
    if(mFnPtrs_.find(Cuts_[i]) == mFnPtrs_.end()){
      cout<<"Didn't find Cut named "<<Cuts_[i]<<endl;
      abort();
    }
    CutFns_[i] = mFnPtrs_.find(Cuts_[i])->second;
  }
  
}

inline void
WZAnalyzer::ResetCounters(){
  results_.assign(NCuts_,wprime::FilterEff());
}

//--------------------------------------------------------------
void WZAnalyzer::Declare_Histos(TFileDirectory & dir)
{
  if(debugme) printf("Declare histos\n");

  DeclareHistoSet("hWZMass", "Reconstructed WZ Invariant Mass",
                  "M_{WZ} (GeV)", 1200, 0, 1200, "GeV", hWZMass,dir);
  DeclareHistoSet("hWZ3e0muMass", "Reconstructed WZ(3e0#mu) Invariant Mass",
                  "M_{WZ}^{3e0#mu} (GeV)", 1200, 0, 1200, "GeV", hWZ3e0muMass,dir);
  DeclareHistoSet("hWZ2e1muMass", "Reconstructed WZ(2e1#mu) Invariant Mass",
                  "M_{WZ}^{2e1#mu} (GeV)", 1200, 0, 1200, "GeV", hWZ2e1muMass,dir);
  DeclareHistoSet("hWZ1e2muMass", "Reconstructed WZ(1e2#mu) Invariant Mass",
                  "M_{WZ}^{1e2#mu} (GeV)", 1200, 0, 1200, "GeV", hWZ1e2muMass,dir);
  DeclareHistoSet("hWZ0e3muMass", "Reconstructed WZ(0e3#mu) Invariant Mass",
                  "M_{WZ}^{0e3#mu} (GeV)", 1200, 0, 1200, "GeV", hWZ0e3muMass,dir);

//Q=M_{WZ} - M_W - M_Z
  DeclareHistoSet("hQ", "Q=M_{WZ} - M_{W} - M_{Z}",
                  "Q (GeV)", 50, 0, 500, "GeV", hQ,dir);
  DeclareHistoSet("hWZTransMass", "Reconstructed WZ Transverse Mass",
                  "M_{WZ}^{T} (GeV)", 100, 0, 1000, "GeV", hWZTransMass,dir);
//WZpt Histos
  DeclareHistoSet("hWZpt", "Reconstructed WZ Transverse Momentum",
                  "p_{WZ}^{T} (GeV)", 50, 0, 500, "GeV", hWZpt,dir);
//Ht Histos
  DeclareHistoSet("hHt", "H_{T}", 
                  "Lepton Pt Sum: H_{T} (GeV)", 80, 0, 800, "GeV", hHt,dir);
  DeclareHistoSet("hTriLepMass", "hTriLepMass",
                  "Trilepton Invariant Mass", 100, 0., 1000., "GeV", hTriLepMass, dir);
  DeclareHistoSet("hEvtType", "Event Type",
                  "N_{#mu}", 4, 0, 4, "NONE", hEvtType,dir);
  DeclareHistoSet("hEvtTypeP", "Event Type for Q=+1",
                  "N_{#mu},W^{+}", 4, 0, 4, "NONE", hEvtTypeP,dir);
  DeclareHistoSet("hEvtTypeM", "Event Type for Q=-1",
                  "N_{#mu},W^{-}", 4, 0, 4, "NONE", hEvtTypeM,dir);
//Lead Lepton Pt
  DeclareHistoSet("hLeadPt", "Leading Lepton Pt",
                  "p_{T}^{Max}", 40, 0, 400., "GeV", hLeadPt,dir);
  DeclareHistoSet("hLeadPtZee", "Leading Lepton Pt Zee",
                  "p_{T}^{Max, ee}", 40, 0, 400., "GeV", hLeadPtZee,dir);
  DeclareHistoSet("hLeadPtZmm", "Leading Lepton Pt Zmm",
                  "p_{T}^{Max #mu#mu}", 40, 0, 400., "GeV", hLeadPtZmm,dir);
  DeclareHistoSet("hLeadElecPt", "Leading Electron Pt",
                  "p_{T}^{Max e}", 40, 0, 400., "GeV", hLeadElecPt,dir);
  DeclareHistoSet("hLeadMuonPt", "Leading Muon Pt",
                  "p_{T}^{Max #mu}", 40, 0, 400., "GeV", hLeadMuonPt,dir);

///////////////////////////
//Z Mass Histos
  DeclareHistoSet("hZMass" , "Reconstructed Mass of Z",
                  "M_{Z} (GeV)", 30, 60, 120, "GeV", hZMass,dir);
  DeclareHistoSet("hZeeMass","Reconstructed Mass of Zee",
                  "M_{Z}^{ee} (GeV)", 30, 60, 120, "GeV", hZeeMass,dir);
  DeclareHistoSet("hZmmMass","Reconstructed Mass of Z#mu#mu",
                  "M_{Z}^{#mu#mu} (GeV)", 30, 60, 120, "GeV", hZmmMass,dir);
  DeclareHistoSet("hZ3e0muMass" , "Reconstructed Mass of Z(3e0#mu)",
                  "M_{Z}^{3e0#mu} (GeV)", 30, 60, 120, "GeV", hZ3e0muMass,dir);
  DeclareHistoSet("hZ2e1muMass" , "Reconstructed Mass of Z(2e1#mu)",
                  "M_{Z}^{2e1#mu} (GeV)", 30, 60, 120, "GeV", hZ2e1muMass,dir);
  DeclareHistoSet("hZ1e2muMass" , "Reconstructed Mass of Z(1e2#mu)",
                  "M_{Z}^{1e2#mu} (GeV)", 30, 60, 120, "GeV", hZ1e2muMass,dir);
  DeclareHistoSet("hZ0e3muMass" , "Reconstructed Mass of Z(0e3#mu)",
                  "M_{Z}^{0e3#mu} (GeV)", 30, 60, 120, "GeV", hZ0e3muMass,dir);
  DeclareHistoSet("hZeeMassTT","Reconstructed MassTT of ZeeTT",
                  "M_{Z}^{ee,TT} (GeV)", 30, 60, 120, "GeV", hZeeMassTT,dir);
  DeclareHistoSet("hZeeMassTF","Reconstructed Mass of ZeeTF",
                  "M_{Z}^{ee,TF} (GeV)", 30, 60, 120, "GeV", hZeeMassTF,dir);
  DeclareHistoSet("hZmmMassTT","Reconstructed Mass of Z#mu#muTT",
                  "M_{Z}^{#mu#mu,TT} (GeV)", 30, 60, 120, "GeV", hZmmMassTT,dir);
  DeclareHistoSet("hZmmMassTF","Reconstructed Mass of Z#mu#muTF",
                  "M_{Z}^{#mu#mu,TF} (GeV)", 30, 60, 120, "GeV", hZmmMassTF,dir);

//Zpt Histos
  DeclareHistoSet("hZpt", "p_{T}^{Z}", 
                  "p_{T}^{Z} (GeV)", 40, 0, 400, "GeV", hZpt,dir);
  DeclareHistoSet("hZeept", "p_{T}^{Z#rightarrowee}", 
                  "p_{T}^{Z#rightarrowee} (GeV)", 40, 0, 400, "GeV", hZeept,dir);
  DeclareHistoSet("hZmmpt", "p_{T}^{Z#rightarrow#mu#mu}", 
                  "p_{T}^{Z#rightarrow#mu#mu} (GeV)", 40, 0, 400, "GeV", hZmmpt,dir);
//MET Histos
  DeclareHistoSet("hMET", "MET",
                  "#slash{E}_{T} (GeV)", 30, 0, 300, "GeV", hMET,dir);
  DeclareHistoSet("hMETee", "MET",
                  "#slash{E}_{T}^{ee} (GeV)", 30, 0, 300, "GeV", hMETee,dir);
  DeclareHistoSet("hMETmm", "MET",
                  "#slash{E}_{T}^{#mu#mu} (GeV)", 30, 0, 300, "GeV", hMETmm,dir);
  DeclareHistoSet("hMET3e0mu", "MET",
                  "#slash{E}_{T}^{3e0#mu} (GeV)", 30, 0, 300, "GeV", hMET3e0mu,dir);
  DeclareHistoSet("hMET2e1mu", "MET",
                  "#slash{E}_{T}^{2e1#mu} (GeV)", 30, 0, 300, "GeV", hMET2e1mu,dir);
  DeclareHistoSet("hMET1e2mu", "MET",
                  "#slash{E}_{T}^{1e2#mu} (GeV)", 30, 0, 300, "GeV", hMET1e2mu,dir);
  DeclareHistoSet("hMET0e3mu", "MET",
                  "#slash{E}_{T}^{0e3#mu} (GeV)", 30, 0, 300, "GeV", hMET0e3mu,dir);

//W Trans Mass Histos
  DeclareHistoSet("hWTransMass", "Reconstructed Transverse Mass of W",
                  "M_{T} (GeV)", 20, 0, 100, "GeV", hWTransMass,dir);
  DeclareHistoSet("hWenuTransMass", "Reconstructed Transverse Mass of We\\nu",
                  "M_{T}^{e#nu} (GeV)", 20, 0, 100, "GeV", hWenuTransMass,dir);
  DeclareHistoSet("hWmnuTransMass", "Reconstructed TransverseMass of W#mu\\nu",
                  "M_{T}^{#mu#nu} (GeV)", 20, 0, 100, "GeV", hWmnuTransMass,dir);
  DeclareHistoSet("hW3e0muTransMass", "Reconstructed Transverse Mass of W(3e0#mu)",
                  "M_{T}^{3e0#mu} (GeV)", 20, 0, 100, "GeV", hW3e0muTransMass,dir);
  DeclareHistoSet("hW2e1muTransMass", "Reconstructed Transverse Mass of W(2e1#mu)",
                  "M_{T}^{2e1#mu} (GeV)", 20, 0, 100, "GeV", hW2e1muTransMass,dir);
  DeclareHistoSet("hW1e2muTransMass", "Reconstructed Transverse Mass of W(1e2#mu)",
                  "M_{T}^{1e2#mu} (GeV)", 20, 0, 100, "GeV", hW1e2muTransMass,dir);
  DeclareHistoSet("hW0e3muTransMass", "Reconstructed Transverse Mass of W(0e3#mu)",
                  "M_{T}^{0e3#mu} (GeV)", 20, 0, 100, "GeV", hW0e3muTransMass,dir);

//Wpt Histos
  DeclareHistoSet("hWpt", "p_{T}^{W}", 
                  "p_{T}^{W} (GeV)", 40, 0, 400, "GeV", hWpt,dir);
  DeclareHistoSet("hWptZee", "p_{T}^{W,Z#rightarrowee}", 
                  "p_{T}^{W,Z#rightarrowee} (GeV)", 40, 0, 400, "GeV", hWptZee,dir);
  DeclareHistoSet("hWptZmm", "p_{T}^{W,Z#rightarrow#mu#mu}", 
                  "p_{T}^{W,Z#rightarrow#mu#mu} (GeV)", 40, 0, 400, "GeV", hWptZmm,dir);

//W Charge Histos
  DeclareHistoSet("hWQ", "Reconstructed Charge of W",
                  "q_{W}", 3, -1, 1, "", hWQ,dir);
  DeclareHistoSet("hWenuQ", "Reconstructed Charge of We\\nu",
                  "q_{W}^{e#nu}", 3, -1.5, 1.5, "", hWenuQ,dir);
  DeclareHistoSet("hWmnuQ", "Reconstructed TransverseMass of W#mu\\nu",
                  "q_{W}^{#mu#nu}", 3, -1.5, 1.5, "", hWmnuQ,dir);
  DeclareHistoSet("hW3e0muQ", "Reconstructed Charge of W(3e0#mu)",
                  "q_{W}^{3e0#mu}", 3, -1.5, 1.5, "", hW3e0muQ,dir);
  DeclareHistoSet("hW2e1muQ", "Reconstructed Charge of W(2e1#mu)",
                  "q_{W}^{2e1#mu}", 3, -1.5, 1.5, "", hW2e1muQ,dir);
  DeclareHistoSet("hW1e2muQ", "Reconstructed Charge of W(1e2#mu)",
                  "q_{W}^{1e2#mu}", 3, -1.5, 1.5, "", hW1e2muQ,dir);
  DeclareHistoSet("hW0e3muQ", "Reconstructed Charge of W(0e3#mu)",
                  "q_{W}^{0e3#mu}", 3, -1.5, 1.5, "", hW0e3muQ,dir);

  DeclareHistoSet("hNLElec", "Number of Loose Electrons in Event",
                  "N_{e}^{Loose}", 10, 0, 10, "NONE", hNLElec,dir);
  DeclareHistoSet("hNLMuon", "Number of Loose Muons in Event",
                  "N_{#mu}^{Loose}", 10, 0, 10, "NONE", hNLMuon,dir);
  DeclareHistoSet("hNLLeps", "Number of Loose Leptons in Event",
                  "N_{l}^{Loose}", 10, 0, 10, "NONE", hNLLeps,dir);
  DeclareHistoSet("hNLLepsZee", "Number of Loose Leptons in Event, Z#rightarrowee",
                  "N_{l}^{Loose,Z#rightarrowee}", 10, 0, 10, "NONE", hNLLepsZee,dir);
  DeclareHistoSet("hNLLepsZmm", "Number of Loose Leptons in Event",
                  "N_{l}^{Loose,Z#rightarrow#mu#mu}", 10, 0, 10, "NONE", hNLLepsZmm,dir);

  DeclareHistoSet("hNTElec", "Number of Tight Electrons in Event",
                  "N_{e}", 10, 0, 10, "NONE", hNTElec,dir);
  DeclareHistoSet("hNTMuon", "Number of Tight Muons in Event",
                  "N_{#mu}", 10, 0, 10, "NONE", hNTMuon,dir);
  DeclareHistoSet("hNTLeps", "Number of Tight Leptons in Event",
                  "N_{l}", 10, 0, 10, "NONE", hNTLeps,dir);

  DeclareHistoSet("hNJets", "Number of Jets in Event",
                  "N_{Jets}", 10, 0, 10, "NONE", hNJets,dir);
  DeclareHistoSet("hNJetsZee", "Number of Jets in Event, Z#rightarrowee",
                  "N_{Jets}^{Z#rightarrowee}", 10, 0, 10, "NONE", hNJetsZee,dir);
  DeclareHistoSet("hNJetsZmm", "Number of Jets in Event, Z#rightarrow#mu#mu",
                  "N_{Jets}^{Z#rightarrow#mu#mu}", 10, 0, 10, "NONE", hNJetsZmm,dir);

  DeclareHistoSet("hNVtxs", "Number of Vertexs in Event",
                  "N_{Vtx}", 50, 0, 50, "NONE", hNVtxs,dir);
  DeclareHistoSet("hNVtxsZee", "Number of Vertexs in Event, Z#rightarrowee",
                  "N_{Vtx}^{Z#rightarrowee}", 50, 0, 50, "NONE", hNVtxsZee,dir);
  DeclareHistoSet("hNVtxsZmm", "Number of Vertexs in Event, Z#rightarrow#mu#mu",
                  "N_{Vtx}^{Z#rightarrow#mu#mu}", 50, 0, 50, "NONE", hNVtxsZmm,dir);

  DeclareHistoSet("hWenuCombRelIso", "Comb Rel Iso of W Electron",
                  "Electron Combined Relative Isolation", 20, 0, 0.2, "NONE", hWenuCombRelIso,dir);
  DeclareHistoSet("hWmnuCombRelIso", "Comb Rel Iso of W Muon",
                  "Muon Combined Relative Isolation", 20, 0, 0.2, "NONE", hWmnuCombRelIso,dir);
  

  tWZCand = dir.make<TTree>("tWZCand", "Analysis Variables after WZCand");//Only 1 for now;
  tWZCand->Branch("WZMass", &WZMass_);
  tWZCand->Branch("EvtType", &evtType_);
  tWZCand->Branch("Ht", &Ht_);
  tWZCand->Branch("Zpt", &Zpt_);
  tWZCand->Branch("Wpt", &Wpt_);
  tWZCand->Branch("weight", &weight_);

///Eff Plots///////
  string title = Form("Expected # of Events / %.0f pb^{-1}",  wprimeUtil_->getLumi_ipb());
  title = title + ";;" + title;
  hNumEvts = NULL; hNumEvts = dir.make<TH1F>("hNumEvts",title.c_str(),NCuts_,0,NCuts_);
  hEffRel  = NULL; hEffRel  = dir.make<TH1F>("hEffRel","Relative Efficiency",NCuts_,0,NCuts_);
  hEffAbs  = NULL; hEffAbs  = dir.make<TH1F>("hEffAbs","Absolute Efficiency",NCuts_,0,NCuts_);

//  DeclareHisto("hNumEvts",title.c_str(),       NCuts_,0,NCuts_,hNumEvts,dir);
//  DeclareHisto("hEffRel","Relative Efficiency",NCuts_,0,NCuts_,hEffRel,dir);
//  DeclareHisto("hEffAbs","Absolute Efficiency",NCuts_,0,NCuts_,hEffAbs,dir);

  for(int i=0; i<NCuts_; ++i) hNumEvts->GetXaxis()->SetBinLabel(i+1,Cuts_[i].c_str());
  for(int i=0; i<NCuts_; ++i) hEffRel ->GetXaxis()->SetBinLabel(i+1,Cuts_[i].c_str());
  for(int i=0; i<NCuts_; ++i) hEffAbs ->GetXaxis()->SetBinLabel(i+1,Cuts_[i].c_str());

}//Declare_Histos

//Fill Histograms
//-----------------------------------------------------------
void WZAnalyzer::Fill_Histos(int index, float weight)
{
//-----------------------------------------------------------
  if(debugme) printf("Filling Histos\n");
  if(wCand_ && zCand_){
    hWZMass[index]->Fill(wzCand_.mass("minPz"), weight);
    if     (evtType_ == 0) hWZ3e0muMass[index]->Fill(wzCand_.mass("minPz"), weight);
    else if(evtType_ == 1) hWZ2e1muMass[index]->Fill(wzCand_.mass("minPz"), weight);
    else if(evtType_ == 2) hWZ1e2muMass[index]->Fill(wzCand_.mass("minPz"), weight);
    else if(evtType_ == 3) hWZ0e3muMass[index]->Fill(wzCand_.mass("minPz"), weight);
    hQ[index]->Fill(Q_, weight); 
    hWZTransMass[index]->Fill(wzCand_.transMass(), weight);
    hWZpt[index]->Fill(wzCand_.pt(), weight);
    hHt[index]->Fill(Ht_, weight);
    hTriLepMass[index]->Fill(TriLepMass_, weight);
    hEvtType[index]->Fill(evtType_, weight);
    if     (wCand_.charge() > 0) hEvtTypeP[index]->Fill(evtType_, weight);
    else if(wCand_.charge() < 0) hEvtTypeM[index]->Fill(evtType_, weight);
    hLeadPt[index]->Fill(LeadPt_, weight);
    hLeadElecPt[index]->Fill(LeadElecPt_, weight);
    hLeadMuonPt[index]->Fill(LeadMuonPt_, weight);
    if     (zCand_.flavor() == PDGELEC){
      hLeadPtZee[index]->Fill(LeadPt_, weight);
      hWptZee[index]->Fill(wCand_.pt(), weight);
      hNLLepsZee[index]->Fill(looseElectrons_.size()+looseMuons_.size(), weight);
      hNJetsZee[index]->Fill(jets_.size(), weight);
      hNVtxsZee[index]->Fill(vertices_.size(), weight);
    }else if(zCand_.flavor() == PDGMUON){ 
      hLeadPtZmm[index]->Fill(LeadPt_, weight);
      hWptZmm[index]->Fill(wCand_.pt(), weight);
      hNLLepsZmm[index]->Fill(looseElectrons_.size()+looseMuons_.size(), weight);
      hNJetsZmm[index]->Fill(jets_.size(), weight);
      hNVtxsZmm[index]->Fill(vertices_.size(), weight);
    }
    if(Cuts_[index] == "ValidWZCand"){//All, Wpt, Zpt, Ht + 1 for starting @ 0
      tWZCand->Fill();
    }
  }
  if(zCand_){
    hZMass[index]->Fill(zCand_.mass(), weight);
    hZpt[index]->Fill(zCand_.pt(), weight);
    if      (zCand_.flavor() == PDGELEC){
      hZeeMass[index]->Fill(zCand_.mass(), weight);
      if(TT) hZeeMassTT[index]->Fill(zCand_.mass(), weight);
      if(TF) hZeeMassTF[index]->Fill(zCand_.mass(), weight);
      hZeept[index]->Fill(zCand_.pt(), weight);
      hMETee[index]->Fill(met_.et(), weight);
    }else if (zCand_.flavor() == PDGMUON){
      hZmmMass[index]->Fill(zCand_.mass(), weight);
      if(TT) hZmmMassTT[index]->Fill(zCand_.mass(), weight);
      if(TF) hZmmMassTF[index]->Fill(zCand_.mass(), weight);
      hMETmm[index]->Fill(met_.et(), weight);
      hZmmpt[index]->Fill(zCand_.pt(), weight);
    }
    if     (evtType_ == 0) hZ3e0muMass[index]->Fill(zCand_.mass(), weight);
    else if(evtType_ == 1) hZ2e1muMass[index]->Fill(zCand_.mass(), weight);
    else if(evtType_ == 2) hZ1e2muMass[index]->Fill(zCand_.mass(), weight);
    else if(evtType_ == 3) hZ0e3muMass[index]->Fill(zCand_.mass(), weight);

    if     (evtType_ == 0) hMET3e0mu[index]->Fill(met_.et(), weight);
    else if(evtType_ == 1) hMET2e1mu[index]->Fill(met_.et(), weight);
    else if(evtType_ == 2) hMET1e2mu[index]->Fill(met_.et(), weight);
    else if(evtType_ == 3) hMET0e3mu[index]->Fill(met_.et(), weight);

  }
  if(wCand_){
    hWTransMass[index]->Fill(wCand_.mt(), weight);
    hWpt[index]->Fill(wCand_.pt(), weight);
    hWQ[index]->Fill(wCand_.charge(), weight);
    if      (wCand_.flavor() == PDGELEC){
      hWenuTransMass[index]->Fill(wCand_.mt(), weight);
      hWenuQ[index]->Fill(wCand_.charge(), weight);
      const heep::Ele& e = *wCand_.elec();
      hWenuCombRelIso[index]->Fill(CalcCombRelIso(e.patEle(), ElecPU(e)), weight);
    }else if (wCand_.flavor() == PDGMUON){
      hWmnuTransMass[index]->Fill(wCand_.mt(), weight);
      hWmnuQ[index]->Fill(wCand_.charge(), weight);
      const TeVMuon& m = *wCand_.muon();
      hWmnuCombRelIso[index]->Fill(m.combRelIsolation03(MuonPU(m)), weight);
    }
    if(evtType_ == 0) hW3e0muTransMass[index]->Fill(wCand_.mt(), weight);
    if(evtType_ == 1) hW2e1muTransMass[index]->Fill(wCand_.mt(), weight);
    if(evtType_ == 2) hW1e2muTransMass[index]->Fill(wCand_.mt(), weight);
    if(evtType_ == 3) hW0e3muTransMass[index]->Fill(wCand_.mt(), weight);

    if(evtType_ == 0) hW3e0muQ[index]->Fill(wCand_.charge(), weight);
    if(evtType_ == 1) hW2e1muQ[index]->Fill(wCand_.charge(), weight);
    if(evtType_ == 2) hW1e2muQ[index]->Fill(wCand_.charge(), weight);
    if(evtType_ == 3) hW0e3muQ[index]->Fill(wCand_.charge(), weight);
  }  
  hMET[index]->Fill(met_.et(), weight);

  hNLElec[index]->Fill(looseElectrons_.size(), weight);
  hNLMuon[index]->Fill(looseMuons_    .size(), weight);
  hNLLeps[index]->Fill(looseElectrons_.size()+looseMuons_.size(), weight);

  hNTElec[index]->Fill(tightElectrons_.size(), weight);
  hNTMuon[index]->Fill(tightMuons_    .size(), weight);
  hNTLeps[index]->Fill(tightElectrons_.size()+tightMuons_.size(), weight);

  hNJets[index]->Fill(jets_.size(), weight);
  hNVtxs[index]->Fill(vertices_.size(), weight);

}//Fill_Histos

int
WZAnalyzer::CountZCands(ZCandV & Zs) const{
  int count =0;
  for(uint i=0; i<Zs.size(); ++i) 
    if( (Zs[i].mass() > minZmass_) && (Zs[i].mass() < maxZmass_)) 
      count++;
  return count;
}

void
WZAnalyzer::CalcZVariables(){
  if (debugme) cout<<"In Calc Z Variables\n";
  // Reconstruct the Z
  float matchptcut = 0.;
  bool minHighPt = false;
  
  matchptcut = 10.;
  ElectronV zElectrons;
  for (size_t i=0; i < looseElectrons_.size(); i++)
    if(PassTriggerEmulation(looseElectrons_[i], matchptcut))
      zElectrons.push_back(looseElectrons_[i]);
  matchptcut = 20.;
  minHighPt = false;
  for (size_t i=0; i < zElectrons.size(); i++){
    if(PassTriggerEmulation(looseElectrons_[i], matchptcut)){
      minHighPt= true; 
      break;  
    }
  }
  if(!minHighPt) zElectrons.clear();
  
  matchptcut = 8.;
  MuonV zMuons;
  for (size_t i=0; i < looseMuons_.size(); i++)
    if(PassTriggerMatch(looseMuons_[i], matchptcut, triggersToUse_))
      zMuons.push_back(looseMuons_[i]);
  
  matchptcut = 13.;
  minHighPt = false;
  for (size_t i=0; i < zMuons.size(); i++){
    if(PassTriggerMatch(zMuons[i], matchptcut, triggersToUse_)){
      minHighPt= true;
      break;
    }
  }
  if(!minHighPt) zMuons.clear();

  if(debugme) printf("    Contains: %i z electron(s), %i z muon(s)\n",
                     (int)zElectrons.size(), (int)zMuons.size());

  ZCandV zeeCands = getZCands(zElectrons, maxZMassDiff_, false);
  removeLowLepPtCands(zeeCands, minZeePt1_, minZeePt2_);
  ZCandV zmmCands = getZCands(zMuons    , maxZMassDiff_, false);
  removeLowLepPtCands(zmmCands, minZmmPt1_, minZmmPt2_);
  ZCandV zCands;
  zCands.insert(zCands.end(), zeeCands.begin(), zeeCands.end());
  zCands.insert(zCands.end(), zmmCands.begin(), zmmCands.end());
  removeWorstCands(zCands, minZmass_, maxZmass_);
  sort(zCands.begin(), zCands.end(), closestToZMass());
  removeOverlapping(zCands);
  zCand_ = zCands.size() ? zCands[0] : ZCandidate();

  ZCandV zCandsAll;
  if(zCand_){
    ZCandV zeeCandsAll = getZCands(looseElectrons_, maxZMassDiff_, false);
    //removeLowLepPtCands(zeeCandsAll, minZeePt1_, minZeePt2_);
    ZCandV zmmCandsAll = getZCands(looseMuons_    , maxZMassDiff_, false);
    //removeLowLepPtCands(zmmCandsAll, minZmmPt1_, minZmmPt2_);

    zCandsAll.insert(zCandsAll.end(), zeeCandsAll.begin(), zeeCandsAll.end());
    zCandsAll.insert(zCandsAll.end(), zmmCandsAll.begin(), zmmCandsAll.end());
    removeWorstCands(zCandsAll, minZmass_, maxZmass_);

    bool tight1=false, tight2=false;
    if(zCand_.flavor() == PDGELEC){
      for(uint i=0; i<tightElectrons_.size(); ++i){
        if(!tight1 && Match(tightElectrons_[i], *zCand_.daughter(0))) tight1 = true;
        if(!tight2 && Match(tightElectrons_[i], *zCand_.daughter(1))) tight2 = true;
      }
    }else if(zCand_.flavor() == PDGMUON){
      for(uint i=0; i<tightMuons_.size(); ++i){
        if(!tight1 && Match(tightMuons_[i], *zCand_.daughter(0))) tight1 = true;
        if(!tight2 && Match(tightMuons_[i], *zCand_.daughter(1))) tight2 = true;
      } 
    }
    TT = tight1 && tight2;
    TF = (tight1 && !tight2) || (!tight1 && tight2);
  }
  zCandsAll.insert(zCandsAll.begin(), zCand_);
  removeOverlapping(zCandsAll);
  numZs_ = CountZCands(zCandsAll); 
  Zpt_ = zCand_.pt();

  if(debugme){
    printf("    Contains: %i Z candidate(s)\n", (int)zCandsAll.size());
    PrintEventLeptons();
    PrintEventDetails();
  }
}

inline void
WZAnalyzer::CalcWVariables(){
  if (debugme) cout<<"In Calc W Variables\n";
  wCand_ = getWCand(tightElectrons_, tightMuons_, met_, zCand_, minDeltaR_);
  Wpt_ = wCand_.pt();
  if(debugme) printf("    Contains: %i tight W candidate(s)\n", (bool)wCand_);
  if(debugme){
    PrintEventLeptons(); 
    PrintEventDetails();
  }
}

inline void
WZAnalyzer::CalcWElecVariables(){
  if (debugme) cout<<"In Calc W Elec Variables\n";
  wCand_ = getWCand(tightElectrons_, met_);
  if(debugme) printf("    Contains: %i tight W candidate(s)\n", (bool)wCand_);
}

inline void
WZAnalyzer::CalcWMuonVariables(){
  if (debugme) cout<<"In Calc W Muon Variables\n";
  wCand_ = getWCand(tightMuons_, met_);
  if(debugme) printf("    Contains: %i tight W candidate(s)\n", (bool)wCand_);
}

inline void
WZAnalyzer::CalcWZVariables(){
  if (debugme) cout<<"In Calc WZ Variables\n";
  wzCand_ = (zCand_ && wCand_) ? WZCandidate(zCand_, wCand_) : WZCandidate();
  WZMass_ = wzCand_.mass("minPz");
  Q_ = (zCand_ && wCand_) ? Calc_Q() : -999.;
  if(debugme) PrintEventDetails();
}

void
WZAnalyzer::CalcEventVariables(){
  if (debugme) cout<<"In Calc Event Variables\n";
  evtType_ = (zCand_ && wCand_) ? Calc_EvtType() : -999;
  if(debugme) printf("evt Type: %i, Z Flav: %i, W Flav: %i\n", evtType_, (int)zCand_.flavor(), (int)wCand_.flavor());
  LeadPt_ = CalcLeadPt(); 
  LeadElecPt_ = CalcLeadPt(PDGELEC);
  LeadMuonPt_ = CalcLeadPt(PDGMUON);
  Ht_ = (zCand_ && wCand_) ? Calc_Ht() : -999.;
  TriLepMass_ = (zCand_ && wCand_) ? CalcTriLepMass() : -999.;
}

void
WZAnalyzer::DeclareHistoSet(string n, string t, string xtitle,
                            int nbins, float min, float max, string units,
                            vector<TH1F*>& h, TFileDirectory& d){
  h.assign(NCuts_,NULL);

  float binWidth = (max-min)/nbins;
  for(int i=0; i<NCuts_; ++i){
    
    string name = n + "_" + Cuts_[i];
    string title = t + " (After " + Cuts_[i] + " Cut);"; 
    title += xtitle + ";Events";
    if(units.compare("NONE"))
      title += Form(" / %.0f ",binWidth) + units;
    //title += Form(" / %.0f pb^{-1}", wprimeUtil_->getLumi_ipb());
    h[i] = d.make<TH1F>(name.c_str(),title.c_str(),nbins,min,max);
  }
}

//Tabulate results after the cut has been passed
//-----------------------------------------------------------
void WZAnalyzer::Tabulate_Me(int& cut_index, const float& weight)
{
//-----------------------------------------------------------
  if(debugme) cout<<"Tabulating results for cut_index = "
                  <<cut_index<<" = "<<Cuts_[cut_index]<<endl;

//increase the number of events passing the cuts
  hNumEvts->Fill(cut_index,weight);
  
  results_[cut_index].Nsurv_evt_cut_w += weight;
  results_[cut_index].Nsurv_evt_cut++;
//fill the histograms
  Fill_Histos(cut_index,weight);
    
}//Tabulate_Me

void 
WZAnalyzer::eventLoop(edm::EventBase const & event){
  ClearEvtVariables();
  if(debugme) WPrimeUtil::PrintEvent(event);

  // Preselection - skip events that don't look promising
  if (doPreselect_){
    if(debugme) cout<<"Testing Preselection...\n";
    if (getProduct<double>(event, 
                           "wzPreselectionProducer:ZMassDiff") > 30.0 ||
        getProduct<double>(event, 
                           "wzPreselectionProducer:highestLeptonPt") < 10 ||
        getProduct<vector<uint> >(event, 
                                  "wzPreselectionProducer:nLeptonsEid")[5] < 3)
      return;
  }
//  // Get information about flavorHistory
//  const uint flavorType = getProduct<uint>(event, "flavorHistoryFilter", 0);
//  if(debugme) printf("    FlavorType: %i", flavorType);

  // Get leptons
  event.getByLabel(electronsLabel_,patElectronsH_);
  //const vector<pat::Electron> patElectrons = getProduct<vector<pat::Electron> >(event, electronsLabel_);
  event.getByLabel(muonsLabel_,patMuonsH_);
  //const vector<pat::Muon    > patMuons     = getProduct<vector<pat::Muon    > >(event, muonsLabel_);
  WPrimeUtil::getLeptonsMET(patElectronsH_, electrons_,
                            patMuonsH_, muonAlgo_, muons_,
                            metH_, useAdjustedMET_, met_,
                            pfCandidatesH_);
  if(debugme) printf("    Contains: %i electron(s), %i muon(s)\n",
                          (int)electrons_.size(), (int)muons_.size());

  rhoFastJet_ = getProduct<double>(event,"kt6PFJets:rho");

  // Make vectors of leptons passing various criteria
  for (size_t i = 0; i < electrons_.size(); i++) {
    if(EMuOverlap(electrons_[i].patEle(), muons_)) continue;
    const float pu = ElecPU(electrons_[i]);
    if (looseElectron_(electrons_[i].patEle(), electronResult_, pu))
      looseElectrons_.push_back(electrons_[i]);

    if (tightElectron_(electrons_[i].patEle(), electronResult_, pu))
      tightElectrons_.push_back(electrons_[i]);
  }

  for (size_t i = 0; i < muons_.size(); i++) {
    const float pu = MuonPU(muons_[i]);
    if (looseMuon_(muons_[i], muonResult_,pu))
      looseMuons_.push_back(muons_[i]);

    if (tightMuon_(muons_[i], muonResult_,pu))
      tightMuons_.push_back(muons_[i]);
  }
  if(looseElectrons_.size() + looseMuons_.size() == 0) return;

  if(debugme){
    PrintLeptons();
    printf("    Contains: %i loose electron(s), %i loose muon(s)\n",
           (int)looseElectrons_.size(), (int)looseMuons_.size());
    printf("    Contains: %i tight electron(s), %i tightmuon(s)\n",
           (int)tightElectrons_.size(), (int)tightMuons_.size());
  }

  //Get Jets
  jets_  = getProduct< JetV>(event,jetsLabel_);

  //Get Trigger 
  triggerEvent_ = getProduct<pat::TriggerEvent>(event,hltEventLabel_); 

  //Get Vertex
  vertices_ = getProduct<vector<reco::Vertex> >(event,vertexLabel_);

  float PU_Weight = 1.;
  if(!wprimeUtil_->runningOnData()){//Don't do this for data
    if(debugme){
      GenParticleV genParticles = getProduct<GenParticleV>(event, "genParticles");
      const reco::Candidate * genZ = 0;
      const reco::Candidate * genW = 0;
      for (size_t i = 0; i < genParticles.size(); i++){
        if (abs(genParticles[i].pdgId()) == 23){
          genZ = & genParticles[i];
          cout<<"Mass of gen Z is "<<genZ->mass()<<endl;
        }else if (abs(genParticles[i].pdgId()) == 24){ 
          genW = & genParticles[i];
        cout<<"Mass of gen W is "<<genW->mass()<<endl;
        }
      }
    }
    PupInfo_ = getProduct<std::vector< PileupSummaryInfo > >(event, pileupLabel_);   
    PU_Weight = wprimeUtil_->getPUWeight3BX(PupInfo_);

    if(debugme) 
      cout<<" PU Weight: "<<PU_Weight<<endl;   

  }//MC Only If

  if(wprimeUtil_->DebugEvent(event)){
    cout<<"This is a debug event\n";
    PrintPassingEvent(event);
    PrintDebugEvent();
  }

  weight_ = wprimeUtil_->getWeight()*PU_Weight;
  if(!PassCuts(weight_)) return;
  if(wprimeUtil_->runningOnData()){
    cout<<" The following data events passed All Cuts!!!\n";
    PrintPassingEvent(event);
    if(debugme) PrintEventLeptons();
    cout<<" ------------------\n";
  }
  if(debugme) PrintEventLeptons();
}

void
WZAnalyzer::PrintEventFull(edm::EventBase const & event) const{
  WPrimeUtil::PrintEvent(event);
  WPrimeUtil::PrintPassingTriggers(triggerEvent_,triggersToUse_);
  PrintEventDetails();
  PrintEventLeptons();
}

void WZAnalyzer::PrintPassingEvent(edm::EventBase const & event){
  PrintEventToFile(event);
  WPrimeUtil::PrintEvent(event);
  PrintEventDetails();
}

void WZAnalyzer::PrintDebugEvent() const{
  WPrimeUtil::PrintPassingTriggers(triggerEvent_,triggersToUse_);
  PrintEventDetails();
  PrintEventLeptons();
  PrintLeptons();
}

void WZAnalyzer::PrintEventToFile(edm::EventBase const & event){
  outCandEvt_<<event.id().run()<<":"
             <<event.id().luminosityBlock()<<":"
             <<event.id().event()<<endl;
}

void WZAnalyzer::PrintEventDetails() const{
  if(zCand_){
    cout<<" Z Flavor: "<<zCand_.flavor()
        <<" Z Mass: "<<zCand_.mass()
        <<" Z lep1 pt "<<ZLepPt(0)
        <<" Z lep2 pt "<<ZLepPt(1)
        <<endl;
  }
  if(wCand_){
    cout<<" W Flavor: "<<wCand_.flavor()
        <<" W MT: "<<wCand_.mt()
        <<" W lep pt "<<WLepPt()
        <<" pfMet et: "<<met_.et()
        <<" pfMet phi: "<<met_.phi()
        <<endl;
  }
  if(zCand_ && wCand_ && wzCand_.mass("minPz")>0.){
    cout<<" WZ Mass: "<<wzCand_.mass("minPz")
        <<" Neu Pz: "<<wzCand_.neutrinoPz("minPz")
        <<" Ht: "<<Ht_
        <<" Zpt: "<<zCand_.pt()
        <<" Wpt: "<<wCand_.pt()
        <<endl;
  }
  return;
}

void
WZAnalyzer::PrintEventLeptons() const{
  if     (zCand_.flavor() == PDGELEC){
//    PrintElectron(*zCand_.elec1(), PDGZ);
//    PrintElectron(*zCand_.elec2(), PDGZ);
    PrintElectron(Find(*zCand_.daughter(0), electrons_), PDGZ);
    PrintElectron(Find(*zCand_.daughter(1), electrons_), PDGZ);
  }else if(zCand_.flavor() == PDGMUON){
    PrintMuon(Find(*zCand_.daughter(0), muons_), PDGZ);
    PrintMuon(Find(*zCand_.daughter(1), muons_), PDGZ);
//    PrintMuon(*zCand_.muon1(), PDGZ);
//    PrintMuon(*zCand_.muon2(), PDGZ);
  }

  if     (wCand_.flavor() == PDGELEC){   
    PrintElectron(Find(*wCand_.daughter(0), electrons_), PDGW);
//    PrintElectron(*wCand_.elec(), PDGW);
  }else if(wCand_.flavor() == PDGMUON){
    PrintMuon(Find(*wCand_.daughter(0), muons_), PDGW);
//    PrintMuon    (*wCand_.muon(), PDGW);
  }
}

void
WZAnalyzer::PrintLeptons() const{
  cout<<"----All Electrons:"<<electrons_.size()<<"------\n";
  for(uint i=0; i<electrons_.size(); ++i) PrintElectron(electrons_[i]);
  cout<<"----All Muons:"<<muons_.size()<<"------\n";
  for(uint i=0; i<muons_    .size(); ++i) PrintMuon    (muons_[i]);
  cout<<"----Loose Electrons:"<<looseElectrons_.size()<<"------\n";
  for(uint i=0; i<looseElectrons_.size(); ++i) PrintElectron(looseElectrons_[i]);
  cout<<"----Loose Muons:"<<looseMuons_.size()<<"------\n";
  for(uint i=0; i<looseMuons_    .size(); ++i) PrintMuon    (looseMuons_[i]);
  cout<<"----Tight Electrons:"<<tightElectrons_.size()<<"------\n";
  for(uint i=0; i<tightElectrons_.size(); ++i) PrintElectron(tightElectrons_[i]);
  cout<<"----Tight Muons:"<<tightMuons_.size()<<"------\n";
  for(uint i=0; i<tightMuons_    .size(); ++i) PrintMuon    (tightMuons_[i]);
  cout<<"----------------------\n";
  cout<<"----------------------\n";
}

void
WZAnalyzer::PrintElectron(const heep::Ele& elec, int parent) const{
  cout << setiosflags(ios::fixed) << setprecision(3);
  if     (parent == PDGZ) cout<<"-----Electron from Z-------------------------"<<endl;
  else if(parent == PDGW) cout<<"-----Electron from W-------------------------"<<endl;
  else                    cout<<"-----Electron from ?-------------------------"<<endl;
  cout<<" Elec ScEt: "<<elec.et()<<endl; //ScEt
  if(!elec.isPatEle()){
    cout<<"Not a pat electron, whye???\n";
    return;
  }
  PrintElectron(elec.patEle(), parent);
}

void
WZAnalyzer::PrintElectron(const pat::Electron& elec, int parent) const{
  cout<<" Elec Pt: "<<elec.pt()<<endl
      <<" Elec P4Pt: "<<elec.p4().Pt()<<endl
      <<" Elec energy: "<<elec.energy()<<endl
      <<" Elec Charge: "<<elec.charge()<<endl
      <<" Elec Eta: "<<elec.eta()<<", isEB="<<elec.isEB()<<endl //Eta
      <<" Elec Phi: "<<elec.phi()<<endl
      <<" Elec NMiss: "<<elec.gsfTrack().get()->trackerExpectedHitsInner().numberOfHits()<<endl
      <<" Elec Dist: "<<elec.convDist()<<endl
      <<" Elec dCotTheta: "<<elec.convDcot()<<endl
      <<" Elec SigmaNN: "<<elec.sigmaIetaIeta()<<endl //sigmaNN
      <<" Elec dPhi: "<<elec.deltaPhiSuperClusterTrackAtVtx()<<endl //DeltaPhi
      <<" Elec dEta: "<<elec.deltaEtaSuperClusterTrackAtVtx()<<endl //DeltaEta
      <<" Elec HoverE: "<<elec.hadronicOverEm()<<endl// H/E
      <<" Elec EoverP: "<<elec.eSuperClusterOverP()<<endl;// E/P
  if(1 || parent == PDGW){
    cout<<" rhoFastJet: "<<rhoFastJet_<<endl;   
    cout<<" adj rhoFastJet: "<<rhoFastJet_*effectiveElecArea_[elec.isEE()]<<endl;   
    
    cout<<" Elec TrkIso: "<<CalcTrkIso(elec)<<endl;
    cout<<" Elec ECALIso: "<<CalcECalIso(elec)<<endl;
    cout<<" Elec HCALIso: "<<CalcHCalIso(elec)<<endl;
    cout<<" Elec CombRelIso: "<<CalcCombRelIso(elec, ElecPU(elec))<<endl;
    
//  cout<<" Elec WP95: "<<elec.electronID("simpleEleId95relIso")<<endl
//      <<" Elec WP90: "<<elec.electronID("simpleEleId90relIso")<<endl
//      <<" Elec WP85: "<<elec.electronID("simpleEleId85relIso")<<endl
//      <<" Elec WP80: "<<elec.electronID("simpleEleId80relIso")<<endl;
  }
}

void
WZAnalyzer::PrintMuon(const TeVMuon& mu, int parent) const{
  cout << setiosflags(ios::fixed) << setprecision(3);
  if     (parent == PDGZ) cout<<"-----Muon from Z-------------------------"<<endl;
  else if(parent == PDGW) cout<<"-----Muon from W-------------------------"<<endl;
  else                    cout<<"-----Muon from ?-------------------------"<<endl;
  reco::TrackRef gm = mu.globalTrack();
  cout<<" Muon Pt: "  <<mu.pt()<<endl
      <<" Muon Global Pt: "  <<mu.pt(kGLOBAL)<<endl
      <<" Muon Inner  Pt: "  <<mu.pt(kINNER)<<endl
      <<" Muon TPFMS  Pt: "  <<mu.pt(kTPFMS)<<endl
      <<" Muon COCKTA Pt: "  <<mu.pt(kCOCKTAIL)<<endl
      <<" Muon Picky  Pt: "  <<mu.pt(kPICKY)<<endl
      <<" Muon TeV    Pt: "  <<mu.pt(kTEV)<<endl
      <<" Muon DYT    Pt: "  <<mu.pt(kDYT)<<endl
      <<" Muon Pat    Pt: "  <<mu.pt(kPAT)<<endl
      <<" Muon Charge: "<<mu.charge()<<endl
      <<" Muon Eta: " <<mu.eta()<<endl
      <<" Muon Phi: " <<mu.phi()<<endl
      <<" Muon Dxy: " <<mu.userFloat("d0")<<endl //Dxy
      <<" Muon NormX2: "<<gm->normalizedChi2()<<endl //NormX2
      <<" Muon NPix: "  <<gm->hitPattern().numberOfValidPixelHits()<<endl //Npixhit
      <<" Muon NTrk: "  <<gm->hitPattern().numberOfValidTrackerHits()<<endl //Ntrk hit
      <<" Muon NMatches: "<<mu.numberOfMatches()<<endl //MuonStations
      <<" Muon Hits: "  <<gm->hitPattern().numberOfValidMuonHits()<<endl; //Muon Hits

  if(1 || parent == PDGW){
    cout<<" Muon EcalIso: "<<mu.isolationR03().emEt<<endl
        <<" Muon HcalIso: "<<mu.isolationR03().hadEt<<endl
        <<" Muon TrkIso: "<<mu.isolationR03().sumPt<<endl
        <<" Muon PU Corr: "<<MuonPU(mu)<<endl
        <<" Muon RelIso: "<<mu.combRelIsolation03(MuonPU(mu))<<endl;// CombRelIso
  }
}

float
WZAnalyzer::CalcLeadPt(int type) const{
  if(type){
    float leadpt = -999.;
    if(wCand_ && wCand_.flavor() == type) 
      leadpt = TMath::Max(leadpt, WLepPt());
    if(zCand_ && zCand_.flavor() == type){
      leadpt = TMath::Max(leadpt, ZLepPt(0));
      leadpt = TMath::Max(leadpt, ZLepPt(1));
    }
    return leadpt;
  }
  return TMath::Max(CalcLeadPt(PDGELEC), CalcLeadPt(PDGMUON));
}

/////////////////Accessors///////////////////////

/////////////////Modifiers///////////////////////

inline void WZAnalyzer::SetCandEvtFile(const string& s){
  outCandEvt_.open(s.c_str());
  WPrimeUtil::CheckStream(outCandEvt_, s);
}

/////////////////Cuts///////////////////////
bool
WZAnalyzer::PassCuts(const float& weight){
  if (debugme) cout<<"In Pass Cuts\n";
  
  for(int i=0; i<NCuts_; ++i){
    if(!(this->*CutFns_[i])()) return false;
    Tabulate_Me(i,weight); 
  }
  return true;
}

inline bool WZAnalyzer::PassNoCut(){ 
  return true;
}

inline bool WZAnalyzer::PassWLepTightCut(){
  if(wCand_.flavor() == PDGELEC){
    const heep::Ele & e = *wCand_.elec();
    return tightElectron_(e.patEle(),electronResult_,ElecPU(e));
  }else if(wCand_.flavor() == PDGMUON){
    const TeVMuon & m = *wCand_.muon();
    return tightMuon_(m, muonResult_,MuonPU(m));
}
  return false;
}
/*
inline bool WZAnalyzer::PassWLepIsoCut() const{
  if(wCand_.flavor() == PDGELEC) return PassElecTightCombRelIsoCut(Find(*wCand_.daughter(0))); 
  if(wCand_.flavor() == PDGMUON) return Find(*wCand_.daughter(0)).combRelIsolation03();
  return false;
}
*/
//Trigger requirements
//-----------------------------------------------------------
bool WZAnalyzer::PassTriggersCut(){
  if(debugme) cout<<"Trigger requirements"<<endl;
  //Apply the trigger if running on data or MC 
  //If MC, apply if no Z or if Z exists, zCand == PDGMuon)
  if(wprimeUtil_->runningOnData() || !zCand_ || zCand_.flavor() == PDGMUON){
    return WPrimeUtil::PassTriggersCut(triggerEvent_,triggersToUse_);
  }else{
    return true;//Cory: This is not good, but will pass HLT in the meantime.
  }
  return false;
}//--- PassTriggersCut()

bool
WZAnalyzer::EMuOverlap(const pat::Electron & e, 
                       const MuonV & ms) const{
  //Eliminate electrons that fall within a cone of dR=0.01 around a muon
  for (size_t i = 0; i < ms.size(); i++) {
    if(debugme) cout<<" with muon "<<i<<": "<<reco::deltaR(e, ms[i])<<endl;
    if(reco::deltaR(e, ms[i]) < 0.01) return true;
  }
  return false;
}

inline bool
WZAnalyzer::PassMinNLeptonsCut(){
  return (looseElectrons_.size() + looseMuons_.size()) >= minNLeptons_;
}

inline bool
WZAnalyzer::PassMaxNLeptonsCut(){
  return (looseElectrons_.size() + looseMuons_.size()) <= maxNLeptons_;
}

inline bool
WZAnalyzer::PassValidWandZCut(){
  return PassValidZCut() && PassValidWCut();
}

inline bool
WZAnalyzer::PassValidWCut(){
  CalcWVariables();
  CalcEventVariables();
  return wCand_ && wCand_.mt()>0.;
}

inline bool
WZAnalyzer::PassValidWElecCut(){
  CalcWElecVariables();
  return wCand_ && wCand_.mt()>0.;
}

inline bool
WZAnalyzer::PassValidWMuonCut(){
  CalcWMuonVariables();
  return wCand_ && wCand_.mt()>0.;
}

inline bool
WZAnalyzer::PassValidZCut(){
  CalcZVariables();
  return zCand_ && zCand_.mass()>0.;
}

inline bool
WZAnalyzer::PassValidWZCandCut(){
  CalcWZVariables();
  return wzCand_.mass("minPz")>0.;
}

inline bool
WZAnalyzer::PassNumberOfZsCut(){
  return numZs_ <= maxNumZs_;
}

inline bool
WZAnalyzer::PassLeadingLeptonPtCut(){
  return LeadPt_ > minLeadPt_;
}

inline bool
WZAnalyzer::PassMETCut(){
  return met_.et() > minMET_;
}

////////////////////////////////
/////////Check Z Properties/////
////////////////////////////////
inline bool
WZAnalyzer::PassZMassCut(){
  return (zCand_.mass() > minZmass_) && (zCand_.mass() < maxZmass_);  
}

inline bool
WZAnalyzer::PassZptCut(){
  return zCand_.pt() > minZpt_;
}

bool
WZAnalyzer::PassZLepPtCut(){
  if     (zCand_.flavor() == PDGELEC){
    return ( max(ZLepPt(0),ZLepPt(1)) > minZeePt1_ && 
             min(ZLepPt(0),ZLepPt(1)) > minZeePt2_ );
  }else if(zCand_.flavor() == PDGMUON){

    return ( max(ZLepPt(0),ZLepPt(1)) > minZmmPt1_ && 
             min(ZLepPt(0),ZLepPt(1)) > minZmmPt2_ );
  }
  return false;
}

bool
WZAnalyzer::PassZLepTriggerMatchCut(){
  if     (zCand_.flavor() == PDGELEC){ 
    const heep::Ele& e1 = *zCand_.elec1();
    const heep::Ele& e2 = *zCand_.elec2();
    return (PassTriggerEmulation(e1) && PassTriggerEmulation(e2));
  }else if(zCand_.flavor() == PDGMUON){
    return true;//Trigger matching now before making zs
    //TeVMuon& m1 = Find(*zCand_.daughter(0));
    //TeVMuon& m2 = Find(*zCand_.daughter(1));
    //return PassTriggerMatch(m1, m2); 
  }
  return false;
}

bool 
WZAnalyzer::PassTriggerMatch(const pat::Electron & p, const float cut, const vstring& triggers) const{
  for (size_t i=0; i < triggers.size(); ++i){
    if (p.triggerObjectMatchesByPath(triggers[i], true, false).size() > 0){
      const pat::TriggerObjectStandAlone * trigRef = p.triggerObjectMatchByPath(triggers[i], true, false);
      if(trigRef->et() > cut) return true;
    }
  }
  return false;
}

bool 
WZAnalyzer::PassTriggerMatch(const TeVMuon & p, const float cut, const vstring& triggers) const{
  for(uint i=0; i<p.triggerObjectMatches().size(); ++i){
    vector<string> names = p.triggerObjectMatches()[i].pathNames(true, false);
    for(uint j=0; j<names.size(); ++j){
      for (size_t k=0; k < triggers.size(); ++k){
        if(WPrimeUtil::SameTrigger(names[j], triggers[k])){
          if (p.triggerObjectMatchesByPath(names[j], true, false).size() > 0){
            if(p.triggerObjectMatchByPath(names[j], true, false)->pt() > cut) return true;
          }
        }
      }
    }
  }
  return false;
}


inline bool
WZAnalyzer::PassTriggerMatch(const heep::Ele& e1, const heep::Ele& e2) const{
  return (e1.patEle().triggerObjectMatches().size() > 0 &&
          e2.patEle().triggerObjectMatches().size() > 0 &&
          max(e1.patEle().triggerObjectMatches()[0].et(), e2.patEle().triggerObjectMatches()[0].et()) > 17. &&
          min(e1.patEle().triggerObjectMatches()[0].et(), e2.patEle().triggerObjectMatches()[0].et()) > 8.);
}

inline bool
WZAnalyzer::PassTriggerMatch(const TeVMuon& m1, const TeVMuon& m2) const{
  return (m1.triggerObjectMatches().size() > 0 &&
          m2.triggerObjectMatches().size() > 0 &&
          max(m1.triggerObjectMatches()[0].pt(), m2.triggerObjectMatches()[0].pt()) > 13. &&
          min(m1.triggerObjectMatches()[0].pt(), m2.triggerObjectMatches()[0].pt()) > 8.);
}


////////////////////////////////
/////////Check W Properties/////
////////////////////////////////

//Check W Transverse Mass
//-------------------s----------------------------------------
inline bool WZAnalyzer::PassWtransMassCut(){
  return wCand_.mt() > minWtransMass_;
}//--- PassWtransMassCut

inline bool
WZAnalyzer::PassWptCut(){
  return wCand_.pt() > minWpt_;
}

inline bool
WZAnalyzer::PassWLepPtCut(){
  return WLepPt() > minWlepPt_;
}

////////////////////////////////
/////////Check Electron Properties/////
////////////////////////////////
bool WZAnalyzer::PassTriggerEmulation(const heep::Ele& elec, const float minPt) const{
  if(!elec.patEle().ecalDrivenSeed()) return false;
  if(elec.patEle().pt() < minPt) return false;
  float e = elec.patEle().energy() * 
    fabs(sin(elec.patEle().superCluster()->position().theta()));
  if(elec.isEB()){
      return elec.patEle().sigmaIetaIeta() < 0.014 &&
        elec.patEle().hadronicOverEm() < 0.15 && 
      elec.patEle().dr03EcalRecHitSumEt() / e < 0.2  &&
        elec.patEle().dr03HcalTowerSumEt() / e < 0.19;
  }else if(elec.patEle().isEE()){
      return elec.patEle().sigmaIetaIeta() < 0.035 &&
        elec.patEle().hadronicOverEm() <0.10 && 
        elec.patEle().dr03EcalRecHitSumEt() / e < 0.2 &&
        elec.patEle().dr03HcalTowerSumEt() / e < 0.19;
  }
  return false;
}

//Check Ht Properties
//-----------------------------------------------------------
inline bool WZAnalyzer::PassHtCut(){
//-----------------------------------------------------------
  if(debugme) cout<<"Check Ht Cuts"<<endl;
  return Ht_ > minHt_;   
}//--- PassHtCut

///////////////////////////////////
//Fake Rate Cuts
inline bool WZAnalyzer::PassWFlavorElecCut(){
  return wCand_ && wCand_.flavor() == PDGELEC;
}

inline bool WZAnalyzer::PassWFlavorMuonCut(){
  return wCand_ && wCand_.flavor() == PDGMUON;
}

bool WZAnalyzer::PassFakeEvtCut(){
  if(looseElectrons_.size() != 1) return false;
  if(looseMuons_    .size() != 1) return false;
  if(looseMuons_[0].charge() != looseElectrons_[0].charge()) return false;
  return true;
}

//Pass Fake Lepton Tag Cut
//-----------------------------------------------------------
bool WZAnalyzer::PassFakeLeptonTagCut(){
  if(wCand_.flavor() == PDGELEC){
    return tightElectron_(looseElectrons_[0].patEle(), electronResult_,ElecPU(looseElectrons_[0]));
  }else if(wCand_.flavor() == PDGMUON){
    return tightMuon_(looseMuons_[0],muonResult_,MuonPU(looseMuons_[0]));
  }
  return true;
}//--- Tag Cut

//Pass Fake Lepton Probe Cut
//-----------------------------------------------------------
bool WZAnalyzer::PassFakeLeptonProbeLooseCut(){
  if(wCand_.flavor() == PDGELEC){
    return looseMuon_(looseMuons_[0], muonResult_, MuonPU(looseMuons_[0])); //Check the other lepton
  }else if(wCand_.flavor() == PDGMUON){
    return looseElectron_(looseElectrons_[0].patEle(), electronResult_,ElecPU(looseElectrons_[0]));
  }
  return true;
}//--- Probe Cut

bool WZAnalyzer::PassFakeLeptonProbeTightCut(){
  if(wCand_.flavor() == PDGELEC){
    return tightMuon_(looseMuons_[0], muonResult_, MuonPU(looseMuons_[0])); //Check the other lepton
  }else if(wCand_.flavor() == PDGMUON){
      return tightElectron_(looseElectrons_[0].patEle(), electronResult_, ElecPU(looseElectrons_[0]));
  }
  return true;
}//--- Probe Cut

///////////////////////////////////

//Calc Ht
//-----------------------------------------------------------
inline float WZAnalyzer::Calc_Ht() const{
  return WLepPt() + ZLepPt(0) + ZLepPt(1);
}//--- CalcHt

inline float WZAnalyzer::CalcTriLepMass() const{
  return (zCand_.daughter(0)->p4() +
          zCand_.daughter(1)->p4() + 
          wCand_.daughter(0)->p4()).M();
}//--- CalcTriLepMass

inline float WZAnalyzer::Calc_Q() const{
  return wzCand_.mass("minPz") - zCand_.mass() - WMASS;
}

inline int WZAnalyzer::Calc_EvtType() const{
  return (zCand_ && wCand_) ?  2 * (zCand_.flavor() != 11) + (wCand_.flavor() != 11) : -999;
}

inline bool WZAnalyzer::inEE(const TeVMuon& mu) const{
  return fabs(mu.eta()) >= 1.05;
}

///////////////Utilities//////////////////
//--------------------------------------------------------------

inline void
WZAnalyzer::ClearEvtVariables(){
  jets_.clear();
  electrons_.clear();
  looseElectrons_.clear();
  tightElectrons_.clear();
  muons_.clear();
  looseMuons_.clear();
  tightMuons_.clear();
  met_ = pat::MET();
  zCand_ = ZCandidate();
  wCand_ = WCandidate();
  wzCand_ = WZCandidate();
  evtType_ = -999;
  LeadPt_ = -999;
  LeadElecPt_ = -999;
  LeadMuonPt_ = -999;
  WZMass_ = -999;
  Ht_= -999;
  TriLepMass_ = -999;
  Zpt_ = -999;
  Wpt_ = -999;
  Q_ = -999;
  TT = TF = false;
  weight_ = 0;
}

void WZAnalyzer::beginFile(std::vector<wprime::InputFile>::const_iterator fi){
  TFileDirectory dir = wprimeUtil_->getFileService()->mkdir(fi->samplename); 
  Declare_Histos(dir);
  ResetCounters();
}

// operations to be done when closing input file 
// (e.g. print summary)
void WZAnalyzer::endFile(std::vector<wprime::InputFile>::const_iterator fi,
                         ofstream & out){
  WPrimeUtil::tabulateSummary(results_);
  WPrimeUtil::printSummary(fi->samplename, fi->description, Cuts_, results_, out);  
}

void WZAnalyzer::endAnalysis(ofstream & out){
}

float
WZAnalyzer::WLepPt() const{
  if(wCand_.flavor() == PDGELEC){
    return Find(*wCand_.daughter(0), electrons_).patEle().pt();
  }else if(wCand_.flavor() == PDGMUON){
    return Find(*wCand_.daughter(0), muons_).pt();
  }
  return -999.;
}

inline float
WZAnalyzer::ZLepPt(int idx) const{
  if(zCand_.flavor() == PDGELEC)
    return Find(*zCand_.daughter(idx), electrons_).patEle().pt();
  else if(zCand_.flavor() == PDGMUON)
    return Find(*zCand_.daughter(idx), muons_).pt();
  return -999.;
}

inline float
  WZAnalyzer::ElecPU(const heep::Ele & e) const{
  return rhoFastJet_*effectiveElecArea_[e.patEle().isEE()];
}

inline float
  WZAnalyzer::MuonPU(const TeVMuon & m) const{
  return rhoFastJet_*effectiveMuonArea_[inEE(m)];
}
