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

  electronsLabel_ = cfg.getParameter<string>("electrons");
  muonsLabel_ = cfg.getParameter<string>("muons");
  metLabel_ = cfg.getParameter<string>("met");

  muonAlgo_ = cfg.getParameter<uint>("muonAlgo");if(debugme) cout<<"Using muon algo "<<muonAlgo_<<endl;
  useAdjustedMET_ = cfg.getParameter<bool>("useAdjustedMET");

  hltEventLabel_ = cfg.getParameter<string>("hltEventTag");
  pileupLabel_ = cfg.getParameter<string>("pileupTag");

  triggersToUse_          = cfg.getParameter<vstring>("triggersToUse");

  maxZMassDiff_ = cfg.getParameter<double>("maxZMassDiff");
  minDeltaR_ = cfg.getParameter<double>("minDeltaR");

  SetCandEvtFile(cfg.getParameter<string  >("candEvtFile" ));
  Num_surv_cut_.assign(NCuts_,0.);

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
}

WZAnalyzer::~WZAnalyzer(){
  outCandEvt_.close();
}

void WZAnalyzer::FillCutFns(){
  mFnPtrs_["NoCuts"] = &WZAnalyzer::PassNoCut;
  mFnPtrs_["HLT"] = &WZAnalyzer::PassTriggersCut;
  mFnPtrs_["MinNLeptons"] = &WZAnalyzer::PassMinNLeptonsCut;
  mFnPtrs_["MaxNLeptons"] = &WZAnalyzer::PassMaxNLeptonsCut;
  mFnPtrs_["EvtSetup"] = &WZAnalyzer::PassEvtSetupCut;
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
  if(debugme) printf("Reset Counters\n");
  for(uint i=0; i<Num_surv_cut_.size(); ++i)
    Num_surv_cut_[i] = 0.;
}

//--------------------------------------------------------------
void WZAnalyzer::Declare_Histos(TFileDirectory & dir)
{
  if(debugme) printf("Declare histos\n");

  DeclareHistoSet("hEvtType", "Event Type",
                  "N_{#mu}", 4, 0, 4, "NONE", hEvtType,dir);

  DeclareHistoSet("hWZInvMass", "Reconstructed WZ Invariant Mass",
                  "m_{WZ} (GeV)", 100, 0, 1000, "GeV", hWZInvMass,dir);
  DeclareHistoSet("hWZ3e0muInvMass", "Reconstructed WZ(3e0\\mu) Invariant Mass",
                  "m_{WZ}^(3e0\\mu) (GeV)", 100, 0, 1000, "GeV", hWZ3e0muInvMass,dir);
  DeclareHistoSet("hWZ2e1muInvMass", "Reconstructed WZ(2e1\\mu) Invariant Mass",
                  "m_{WZ}^(2e1\\mu) (GeV)", 100, 0, 1000, "GeV", hWZ2e1muInvMass,dir);
  DeclareHistoSet("hWZ1e2muInvMass", "Reconstructed WZ(1e2\\mu) Invariant Mass",
                  "m_{WZ}^(1e2\\mu) (GeV)", 100, 0, 1000, "GeV", hWZ1e2muInvMass,dir);
  DeclareHistoSet("hWZ0e3muInvMass", "Reconstructed WZ(0e3\\mu) Invariant Mass",
                  "m_{WZ}^(0e3\\mu) (GeV)", 100, 0, 1000, "GeV", hWZ0e3muInvMass,dir);

  DeclareHistoSet("hWZTransMass", "Reconstructed WZ Transverse Mass",
                  "m_{WZ} (GeV)", 100, 0, 1000, "GeV", hWZTransMass,dir);
//Ht Histos
  DeclareHistoSet("hHt", "H_{T}", 
                  "Lepton Pt Sum: H_{T} (GeV)", 80, 0, 800, "GeV", hHt,dir);
//Wpt Histos
  DeclareHistoSet("hWpt", "p_{T}^{W}", 
                  "p_{T}^{W} (GeV)", 40, 0, 400, "GeV", hWpt,dir);
//Zpt Histos
  DeclareHistoSet("hZpt", "p_{T}^{Z}", 
                  "p_{T}^{Z} (GeV)", 40, 0, 400, "GeV", hZpt,dir);
//MET Histos
  DeclareHistoSet("hMET", "MET",
                  "MET (GeV)", 30, 0, 300, "GeV", hMET,dir);
  DeclareHistoSet("hMETee", "MET",
                  "MET^{ee} (GeV)", 30, 0, 300, "GeV", hMETee,dir);
  DeclareHistoSet("hMETmumu", "MET",
                  "MET^{#mu#mu} (GeV)", 30, 0, 300, "GeV", hMETmumu,dir);
  DeclareHistoSet("hMET3e0mu", "MET",
                  "MET^{3e0#mu} (GeV)", 30, 0, 300, "GeV", hMET3e0mu,dir);
  DeclareHistoSet("hMET2e1mu", "MET",
                  "MET^{2e1#mu} (GeV)", 30, 0, 300, "GeV", hMET2e1mu,dir);
  DeclareHistoSet("hMET1e2mu", "MET",
                  "MET^{1e2#mu} (GeV)", 30, 0, 300, "GeV", hMET1e2mu,dir);
  DeclareHistoSet("hMET0e3mu", "MET",
                  "MET^{0e3#mu} (GeV)", 30, 0, 300, "GeV", hMET0e3mu,dir);

//Z Mass Histos
  DeclareHistoSet("hZMass" , "Reconstructed Mass of Z",
                  "m_{Z} (GeV)", 30, 60, 120, "GeV", hZMass,dir);
  DeclareHistoSet("hZeeMass","Reconstructed Mass of Zee",
                  "m_{Z}^{ee} (GeV)", 30, 60, 120, "GeV", hZeeMass,dir);
  DeclareHistoSet("hZmumuMass","Reconstructed Mass of Z\\mu\\mu",
                  "m_{Z}^{#mu#mu} (GeV)", 30, 60, 120, "GeV", hZmumuMass,dir);
  DeclareHistoSet("hZ3e0muMass" , "Reconstructed Mass of Z(3e0\\mu)",
                  "m_{Z}^{3e0#mu} (GeV)", 30, 60, 120, "GeV", hZ3e0muMass,dir);
  DeclareHistoSet("hZ2e1muMass" , "Reconstructed Mass of Z(2e1\\mu)",
                  "m_{Z}^{2e1#mu} (GeV)", 30, 60, 120, "GeV", hZ2e1muMass,dir);
  DeclareHistoSet("hZ1e2muMass" , "Reconstructed Mass of Z(1e2\\mu)",
                  "m_{Z}^{1e2#mu} (GeV)", 30, 60, 120, "GeV", hZ1e2muMass,dir);
  DeclareHistoSet("hZ0e3muMass" , "Reconstructed Mass of Z(0e3\\mu)",
                  "m_{Z}^{0e3#mu} (GeV)", 30, 60, 120, "GeV", hZ0e3muMass,dir);
  DeclareHistoSet("hZeeMassTT","Reconstructed MassTT of ZeeTT",
                  "m_{Z}^{ee,TT} (GeV)", 30, 60, 120, "GeV", hZeeMassTT,dir);
  DeclareHistoSet("hZeeMassTF","Reconstructed Mass of ZeeTF",
                  "m_{Z}^{ee,TF} (GeV)", 30, 60, 120, "GeV", hZeeMassTF,dir);
  DeclareHistoSet("hZmumuMassTT","Reconstructed Mass of Z\\mu\\muTT",
                  "m_{Z}^{#mu#mu,TT} (GeV)", 30, 60, 120, "GeV", hZmumuMassTT,dir);
  DeclareHistoSet("hZmumuMassTF","Reconstructed Mass of Z\\mu\\muTF",
                  "m_{Z}^{#mu#mu,TF} (GeV)", 30, 60, 120, "GeV", hZmumuMassTF,dir);

//W Trans Mass Histos
  DeclareHistoSet("hWTransMass", "Reconstructed Transverse Mass of W",
                  "m_{T} (GeV)", 20, 0, 100, "GeV", hWTransMass,dir);
  DeclareHistoSet("hWenuTransMass", "Reconstructed Transverse Mass of We\\nu",
                  "m_{T}^{e#nu} (GeV)", 20, 0, 100, "GeV", hWenuTransMass,dir);
  DeclareHistoSet("hWmunuTransMass", "Reconstructed TransverseMass of W\\mu\\nu",
                  "m_{T}^{#mu#nu} (GeV)", 20, 0, 100, "GeV", hWmunuTransMass,dir);
  DeclareHistoSet("hW3e0muTransMass", "Reconstructed Transverse Mass of W(3e0\\mu)",
                  "m_{T}^{3e0#mu} (GeV)", 20, 0, 100, "GeV", hW3e0muTransMass,dir);
  DeclareHistoSet("hW2e1muTransMass", "Reconstructed Transverse Mass of W(2e1\\mu)",
                  "m_{T}^{2e1#mu} (GeV)", 20, 0, 100, "GeV", hW2e1muTransMass,dir);
  DeclareHistoSet("hW1e2muTransMass", "Reconstructed Transverse Mass of W(1e2\\mu)",
                  "m_{T}^{1e2#mu} (GeV)", 20, 0, 100, "GeV", hW1e2muTransMass,dir);
  DeclareHistoSet("hW0e3muTransMass", "Reconstructed Transverse Mass of W(0e3\\mu)",
                  "m_{T}^{0e3#mu} (GeV)", 20, 0, 100, "GeV", hW0e3muTransMass,dir);

//Q=M_{WZ} - M_W - M_Z
  DeclareHistoSet("hQ", "Q=M_{WZ} - M_{W} - M_{Z}",
                  "Q (GeV)", 50, 0, 500, "GeV", hQ,dir);

  DeclareHistoSet("hNLElec", "Number of Loose Electrons in Event",
                  "N_{e}", 10, 0, 10, "NONE", hNLElec,dir);
  DeclareHistoSet("hNLMuon", "Number of Loose Muons in Event",
                  "N_{#mu}", 10, 0, 10, "NONE", hNLMuon,dir);
  DeclareHistoSet("hNLLeps", "Number of Loose Leptons in Event",
                  "N_{l}", 10, 0, 10, "NONE", hNLLeps,dir);

  DeclareHistoSet("hNTElec", "Number of Tight Electrons in Event",
                  "N_{e}", 10, 0, 10, "NONE", hNTElec,dir);
  DeclareHistoSet("hNTMuon", "Number of Tight Muons in Event",
                  "N_{#mu}", 10, 0, 10, "NONE", hNTMuon,dir);
  DeclareHistoSet("hNTLeps", "Number of Tight Leptons in Event",
                  "N_{l}", 10, 0, 10, "NONE", hNTLeps,dir);

  DeclareHistoSet("hNJets", "Number of Jets in Event",
                  "N_{j}", 20, 0, 20, "NONE", hNJets,dir);
  DeclareHistoSet("hNVtxs", "Number of Vertexs in Event",
                  "N_{vtx}", 20, 0, 20, "NONE", hNVtxs,dir);

  DeclareHistoSet("hWenuCombRelIso", "Comb Rel Iso of W Electron",
                  "Electron Combined Relative Isolation", 20, 0, 0.2, "NONE", hWenuCombRelIso,dir);
  DeclareHistoSet("hWmunuCombRelIso", "Comb Rel Iso of W Muon",
                  "Muon Combined Relative Isolation", 20, 0, 0.2, "NONE", hWmunuCombRelIso,dir);
  
  DeclareHistoSet("hLeadPt", "Leading Lepton Pt",
                  "Lead Lepton Pt", 20, 0, 200., "GeV", hLeadPt,dir);
  DeclareHistoSet("hLeadElecPt", "Leading Electron Pt",
                  "Lead Electron Pt", 20, 0, 200., "GeV", hLeadElecPt,dir);
  DeclareHistoSet("hLeadMuonPt", "Leading Muon Pt",
                  "Lead Muon Pt", 20, 0, 200., "GeV", hLeadMuonPt,dir);

  DeclareHistoSet("hMuonTightCombIso", "hMuonTightCombIso",
                  "Muon Combined Isolation", 20, 0., 1., "NONE", hMuonTightCombIso, dir);

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
    hWZInvMass[index]->Fill(wzCand_.mass("minPz"), weight);
    if     (evtType_ == 0) hWZ3e0muInvMass[index]->Fill(wzCand_.mass("minPz"), weight);
    else if(evtType_ == 1) hWZ2e1muInvMass[index]->Fill(wzCand_.mass("minPz"), weight);
    else if(evtType_ == 2) hWZ1e2muInvMass[index]->Fill(wzCand_.mass("minPz"), weight);
    else if(evtType_ == 3) hWZ0e3muInvMass[index]->Fill(wzCand_.mass("minPz"), weight);
    hEvtType[index]->Fill(evtType_, weight);
    hQ[index]->Fill(Q_, weight); 
    hWZTransMass[index]->Fill(wzCand_.transMass(), weight);
    hHt[index]->Fill(Ht_, weight);
    hLeadPt[index]->Fill(LeadPt_, weight);
    hLeadElecPt[index]->Fill(LeadElecPt_, weight);
    hLeadMuonPt[index]->Fill(LeadMuonPt_, weight);
  }
  if(zCand_){
    hZpt[index]->Fill(zCand_.pt(), weight);
    hZMass[index]->Fill(zCand_.mass(), weight);
    if      (zCand_.flavor() == PDGELEC){
      hZeeMass[index]->Fill(zCand_.mass(), weight);
      if(TT) hZeeMassTT[index]->Fill(zCand_.mass(), weight);
      if(TF) hZeeMassTF[index]->Fill(zCand_.mass(), weight);
      hMETee[index]->Fill(met_.et(), weight);
    }else if (zCand_.flavor() == PDGMUON){
      hZmumuMass[index]->Fill(zCand_.mass(), weight);
      if(TT) hZmumuMassTT[index]->Fill(zCand_.mass(), weight);
      if(TF) hZmumuMassTF[index]->Fill(zCand_.mass(), weight);
      hMETmumu[index]->Fill(met_.et(), weight);
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
    hWpt[index]->Fill(wCand_.pt(), weight);
    hWTransMass[index]->Fill(wCand_.mt(), weight);
    if      (wCand_.flavor() == PDGELEC){
      hWenuTransMass[index]->Fill(wCand_.mt(), weight);
      const heep::Ele & e = FindElectron(*wCand_.daughter(0));
      hWenuCombRelIso[index]->Fill(CalcCombRelIso(e.patEle(), ElecPU(e)), weight);
    }else if (wCand_.flavor() == PDGMUON){
      hWmunuTransMass[index]->Fill(wCand_.mt(), weight);
      const TeVMuon m = FindMuon(*wCand_.daughter(0));
      hWmunuCombRelIso[index]->Fill(m.combRelIsolation03(MuonPU(m)), weight);
    }
    if(evtType_ == 0) hW3e0muTransMass[index]->Fill(wCand_.mt(), weight);
    if(evtType_ == 1) hW2e1muTransMass[index]->Fill(wCand_.mt(), weight);
    if(evtType_ == 2) hW1e2muTransMass[index]->Fill(wCand_.mt(), weight);
    if(evtType_ == 3) hW0e3muTransMass[index]->Fill(wCand_.mt(), weight);
  }  
  hMET[index]->Fill(met_.et(), weight);

  hNLElec[index]->Fill(looseElectrons_.size(), weight);
  hNLMuon[index]->Fill(looseMuons_    .size(), weight);
  hNLLeps[index]->Fill(looseElectrons_.size()+looseMuons_.size(), weight);
  hNTElec[index]->Fill(tightElectrons_.size(), weight);
  hNTMuon[index]->Fill(tightMuons_    .size(), weight);
  hNTLeps[index]->Fill(tightElectrons_.size()+tightMuons_.size(), weight);

  hNJets[index]->Fill(jets_.size(), weight);
  hNVtxs[index]->Fill(PU_NumInteractions_, weight);

  for(uint i=0; i<tightMuons_.size(); ++i){
    hMuonTightCombIso[index]->Fill(tightMuons_[i].combRelIsolation03(MuonPU(tightMuons_[i])), weight);
  }
}//Fill_Histos
int
WZAnalyzer::CountZCands(ZCandV & Zs){
  int count =0;
  for(uint i=0; i<Zs.size(); ++i) 
    if( (Zs[i].mass() > minZmass_) && (Zs[i].mass() < maxZmass_)) 
      count++;
  return count;
}

inline void
WZAnalyzer::CalcZVariables(){
  if (debugme) cout<<"In Calc Z Variables\n";
  // Reconstruct the Z
  ZCandV looseZCands = getZCands(looseElectrons_, looseMuons_, maxZMassDiff_);
  zCand_ = looseZCands.size() ? looseZCands[0] : ZCandidate();
  if(debugme) printf("    Contains: %i loose Z candidate(s)\n", (int)looseZCands.size());
  numZs_ = CountZCands(looseZCands); 
  if(debugme) PrintEventLeptons();
}

inline void
WZAnalyzer::CalcWVariables(){
  if (debugme) cout<<"In Calc W Variables\n";
  wCand_ = getWCand(tightElectrons_, tightMuons_, met_, zCand_, minDeltaR_);
  if(debugme) printf("    Contains: %i tight W candidate(s)\n", (bool)wCand_);
  if(debugme) PrintEventLeptons(); 
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
  Q_ = (zCand_ && wCand_) ? Calc_Q() : -999.;
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
  TT = TF = false;
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
  Num_surv_cut_[cut_index] += weight;
//fill the histograms
  Fill_Histos(cut_index,weight);
    
}//Tabulate_Me

//Writing results to a txt file
//--------------------------------------------------------------------------
void WZAnalyzer::printSummary(const string& dir, ofstream & out) 
{ 
//------------------------------------------------------------------------
  if(debugme) cout<<"Writing results to a txt file"<<endl;

  out << setiosflags(ios::fixed) << setprecision(2);
  out<<"$$$$$$$$$$$$$$$$$$$$$$$ Type of sample: "<<dir<<endl;
//  out << " xsec*lumi expected events = " << Nthe_evt << endl;
//  out << " # of evt passing preselection = " << Nexp_evt << " per "<<Form("%.0f", wprimeUtil_->getLumi_ipb())<<" inv pb"<<endl;
  
  for(int i = 0; i < NCuts_; ++i){
    out<<right<<"Cut " << setw(2) << i << "("
               <<left<< setw(15) << Cuts_[i]
               <<right << "): " <<"expected evts = " << setw(10) << Num_surv_cut_[i];
    hNumEvts->SetBinContent(i+1,Num_surv_cut_[i]);
    
//calculate efficiencies
    float eff, deff;
    if(i == 0){
      WPrimeUtil::getEff(eff, deff, Num_surv_cut_[i], Num_surv_cut_[0]);
    }else{
      WPrimeUtil::getEff(eff, deff, Num_surv_cut_[i], Num_surv_cut_[i-1]);
    }
    hEffRel->SetBinContent(i+1,eff*100);
    out << setw(15) <<"\tRelative eff = "<<setw(6)<<eff*100 << " +/- " << setw(6)<<deff*100 << "%";
    WPrimeUtil::getEff(eff, deff, Num_surv_cut_[i], Num_surv_cut_[0]);
    hEffAbs->SetBinContent(i+1,eff*100);
    out << setw(15) <<"\tAbsolute eff = "<<setw(6)<<eff*100 << " +/- " << setw(6)<<deff*100 << "%"
        << endl;
        
  } // loop over different cuts
}//printSummary


void 
WZAnalyzer::eventLoop(edm::EventBase const & event){
  ClearEvtVariables();

  if(debugme) PrintEvent(event);

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

  rhoFastJet_ = getProduct<double>(event,"kt6PFJets:rho");

  // Get leptons
  const vector<pat::Electron> patElectrons = getProduct<vector<pat::Electron> >(event, electronsLabel_);
  const vector<pat::Muon    > patMuons     = getProduct<vector<pat::Muon    > >(event, muonsLabel_);
  jets_  = getProduct< JetV>(event, "selectedPatJets");
  if(debugme) printf("    Contains: %i electron(s), %i muon(s)\n",
                     (int)patElectrons.size(), (int)patMuons.size());

  // Make vectors of leptons passing various criteria
  for (size_t i = 0; i < patElectrons.size(); i++) {
    electrons_.push_back(heep::Ele(patElectrons[i]));   
    if(EMuOverlap(patElectrons[i], patMuons)) continue;
    const float pu = ElecPU(electrons_[i]);
    if (looseElectron_(electrons_[i].patEle(), electronResult_, pu))
      looseElectrons_.push_back(electrons_[i]);

    if (tightElectron_(electrons_[i].patEle(), electronResult_, pu))
      tightElectrons_.push_back(electrons_[i]);
  }
  for (size_t i = 0; i < patMuons.size(); i++) {
    muons_.push_back(TeVMuon(patMuons[i],muonAlgo_));   
    const float pu = MuonPU(muons_[i]);
    if (looseMuon_(muons_[i], muonResult_,pu))
      looseMuons_.push_back(muons_[i]);

    if (tightMuon_(muons_[i], muonResult_,pu))
      tightMuons_.push_back(muons_[i]);

  }

  if(debugme) PrintLeptons();
  if(debugme) printf("    Contains: %i loose electron(s), %i loose muon(s)\n",
                     (int)looseElectrons_.size(), (int)looseMuons_.size());
  if(debugme) printf("    Contains: %i tight electron(s), %i tightmuon(s)\n",
                     (int)tightElectrons_.size(), (int)tightMuons_.size());

  // Get MET
  met_ = getProduct<METV>(event, metLabel_)[0];
  if(useAdjustedMET_){
    if(debugme) cout<<"Before met et: "<<met_.et()<<" met phi: "<<met_.phi()<<endl;
    met_ = AdjustedMET(looseElectrons_,looseMuons_, met_);
    if(debugme) cout<<"After  met et: "<<met_.et()<<" met phi: "<<met_.phi()<<endl;
  }

  triggerEvent_ = getProduct<pat::TriggerEvent>(event,hltEventLabel_); 

  float PU_Weight = 1.;
  if(!wprimeUtil_->runningOnData()){//Don't do this for data
/*    
    GenParticleV genParticles = getUntrackedProduct<GenParticleV>(event, "genParticles");
    const Candidate * genZ = 0;
    const Candidate * genW = 0;
    for (size_t i = 0; i < genParticles.size(); i++){
      if (abs(genParticles[i].pdgId()) == 23) genZ = & genParticles[i];
      else if (abs(genParticles[i].pdgId()) == 24) genW = & genParticles[i];
    }
*/
    PupInfo_ = getProduct<std::vector< PileupSummaryInfo > >(event, pileupLabel_);   
    PU_NumInteractions_ = wprimeUtil_->GetPU1BX(PupInfo_);
    PU_Weight = wprimeUtil_->getLumiWeight(PU_NumInteractions_);
    if(debugme) 
      cout<<" PU_NumInteractions: "<<PU_NumInteractions_
          <<" PU Weight: "<<PU_Weight
          <<endl;   
  }//MC Only If

  if((event.id().run() == 166033 && event.id().luminosityBlock() == 801 && event.id().event() == 1091994700) ||
     (event.id().run() == 166699 && event.id().luminosityBlock() == 664 && event.id().event() == 714495816)){
    cout<<"This is a missing event\n";
    PrintEvent(event);
    PrintLeptons();
  }

/////////////////////
  if(!PassCuts(wprimeUtil_->getWeight()*PU_Weight)) return;
  if(wprimeUtil_->runningOnData()){
    cout<<" The following data events passed All Cuts!!!\n";
    PrintPassingEvent(event);
    if(debugme) PrintEventLeptons();
    cout<<" ------------------\n";
  }
}

void
WZAnalyzer::PrintEventFull(edm::EventBase const & event){
  PrintEvent(event);
  PrintTrigger();
  PrintEventDetails();
  PrintEventLeptons();
}

void WZAnalyzer::PrintPassingEvent(edm::EventBase const & event){
  PrintEventToFile(event);
  PrintEvent(event);
  PrintEventDetails();
}

void WZAnalyzer::PrintDebugEvent(){
  PrintTrigger();
  PrintEventDetails();
  PrintEventLeptons();
}

void WZAnalyzer::PrintEventToFile(edm::EventBase const & event){
  outCandEvt_<<event.id().run()<<":"
             <<event.id().luminosityBlock()<<":"
             <<event.id().event()<<endl;
}

void WZAnalyzer::PrintEvent(edm::EventBase const & event){
  cout<<"run #: "<<event.id().run()
      <<" lumi: "<<event.id().luminosityBlock()
      <<" eventID: "<<event.id().event()<<endl;
}

void WZAnalyzer::PrintEventDetails(){
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
        <<endl;
  }
  if(zCand_ && wCand_ && wzCand_.mass("minPz")>0.){
    cout<<" WZ Mass: "<<wzCand_.mass("minPz")
        <<" Ht: "<<Ht_
        <<" Zpt: "<<zCand_.pt()
        <<" Wpt: "<<wCand_.pt()
        <<endl;
  }
  return;
}

void
WZAnalyzer::PrintEventLeptons(){
  if     (zCand_.flavor() == PDGELEC){
    PrintElectron(FindElectron(*zCand_.daughter(0)), PDGZ);
    PrintElectron(FindElectron(*zCand_.daughter(1)), PDGZ);
  }else if(zCand_.flavor() == PDGMUON){
    PrintMuon(FindMuon(*zCand_.daughter(0)), PDGZ);
    PrintMuon(FindMuon(*zCand_.daughter(1)), PDGZ);
  }

  if     (wCand_.flavor() == PDGELEC){   
    PrintElectron(FindElectron(*wCand_.daughter(0)), PDGW);
  }else if(wCand_.flavor() == PDGMUON) 
    PrintMuon    (FindMuon(*wCand_.daughter(0)), PDGW);
}

void
WZAnalyzer::PrintTrigger(){
  const pat::TriggerPathRefVector acceptedPaths = triggerEvent_.acceptedPaths();
  for (size_t i = 0; i < acceptedPaths.size(); i++){
    string A = acceptedPaths[i]->name();
    for (size_t j = 0; j < triggersToUse_.size(); j++){
      if(SameTrigger(A, triggersToUse_[j])){
        if(acceptedPaths[i]->prescale() == 1  && acceptedPaths[i]->wasAccept()) cout<<"Passed path: "<<A<<endl;
      }
    }
  }
}

void
WZAnalyzer::PrintLeptons(){
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
WZAnalyzer::PrintElectron(const heep::Ele& elec, int parent){
  cout << setiosflags(ios::fixed) << setprecision(3);
  if     (parent == PDGZ) cout<<"-----Electron from Z-------------------------"<<endl;
  else if(parent == PDGW) cout<<"-----Electron from W-------------------------"<<endl;
  else                    cout<<"-----Electron from ?-------------------------"<<endl;
  cout<<" Elec ScEt: "<<elec.et()<<endl; //ScEt
  if(!elec.isPatEle()) return;
  cout<<" Elec Pt: "<<elec.patEle().pt()<<endl
      <<" Elec P4Pt: "<<elec.patEle().p4().Pt()<<endl
      <<" Elec energy: "<<elec.patEle().energy()<<endl
      <<" Elec Charge: "<<elec.patEle().charge()<<endl
      <<" Elec Eta: "<<elec.patEle().eta()<<", isEB="<<elec.patEle().isEB()<<endl //Eta
      <<" Elec NMiss: "<<elec.patEle().gsfTrack().get()->trackerExpectedHitsInner().numberOfHits()<<endl
      <<" Elec Dist: "<<elec.patEle().convDist()<<endl
      <<" Elec dCotTheta: "<<elec.patEle().convDcot()<<endl
      <<" Elec SigmaNN: "<<elec.patEle().sigmaIetaIeta()<<endl //sigmaNN
      <<" Elec dPhi: "<<elec.patEle().deltaPhiSuperClusterTrackAtVtx()<<endl //DeltaPhi
      <<" Elec dEta: "<<elec.patEle().deltaEtaSuperClusterTrackAtVtx()<<endl //DeltaEta
      <<" Elec HoverE: "<<elec.patEle().hadronicOverEm()<<endl// H/E
      <<" Elec EoverP: "<<elec.patEle().eSuperClusterOverP()<<endl;// E/P
  if(1 || parent == PDGW){
    cout<<" rhoFastJet: "<<rhoFastJet_<<endl;   
    cout<<" adj rhoFastJet: "<<rhoFastJet_*effectiveElecArea_[elec.patEle().isEE()]<<endl;   
    
    cout<<" Elec TrkIso: "<<CalcTrkIso(elec.patEle())<<endl;
    cout<<" Elec ECALIso: "<<CalcECalIso(elec.patEle())<<endl;
    cout<<" Elec HCALIso: "<<CalcHCalIso(elec.patEle())<<endl;
    cout<<" Elec CombRelIso: "<<CalcCombRelIso(elec.patEle(), ElecPU(elec))<<endl;
    
//  cout<<" Elec WP95: "<<elec.patEle().electronID("simpleEleId95relIso")<<endl
//      <<" Elec WP90: "<<elec.patEle().electronID("simpleEleId90relIso")<<endl
//      <<" Elec WP85: "<<elec.patEle().electronID("simpleEleId85relIso")<<endl
//      <<" Elec WP80: "<<elec.patEle().electronID("simpleEleId80relIso")<<endl;
  }
}

void
WZAnalyzer::PrintMuon(const TeVMuon& mu, int parent){
  cout << setiosflags(ios::fixed) << setprecision(3);
  if     (parent == PDGZ) cout<<"-----Muon from Z-------------------------"<<endl;
  else if(parent == PDGW) cout<<"-----Muon from W-------------------------"<<endl;
  else                    cout<<"-----Muon from ?-------------------------"<<endl;
  reco::TrackRef gm = mu.globalTrack();
  cout<<" Muon Pt: "  <<mu.pt()<<endl
      <<" Muon Charge: "<<mu.charge()<<endl
      <<" Muon Eta: " <<mu.eta()<<endl
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
WZAnalyzer::CalcLeadPt(int type){
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

inline bool
WZAnalyzer::SameTrigger(string & A, string & B){
  return (B.find("*") == string::npos) ? !A.compare(B) : !A.compare(0, A.size()-1, B, 0, B.size()-1);
}

/////////////////Accessors///////////////////////

/////////////////Modifiers///////////////////////

inline void WZAnalyzer::SetCandEvtFile(string s){
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

inline bool WZAnalyzer::PassEvtSetupCut(){ 
  CalcEventVariables();
  if(debugme) PrintDebugEvent();
  return true;
}

inline bool WZAnalyzer::PassWLepTightCut(){
  if(wCand_.flavor() == PDGELEC){
    const heep::Ele & e = FindElectron(*wCand_.daughter(0));
    return tightElectron_(e.patEle(),electronResult_,ElecPU(e));
  }else if(wCand_.flavor() == PDGMUON){
    const TeVMuon & m = FindMuon(*wCand_.daughter(0));
    return tightMuon_(m, muonResult_,MuonPU(m));
}
  return false;
}
/*
inline bool WZAnalyzer::PassWLepIsoCut(){
  if(wCand_.flavor() == PDGELEC) return PassElecTightCombRelIsoCut(FindElectron(*wCand_.daughter(0))); 
  if(wCand_.flavor() == PDGMUON) return FindMuon(*wCand_.daughter(0)).combRelIsolation03();
  return false;
}
*/
//Trigger requirements
//-----------------------------------------------------------
bool WZAnalyzer::PassTriggersCut()
{
//-----------------------------------------------------------
  if(debugme) cout<<"Trigger requirements"<<endl;
  
  if(wprimeUtil_->runningOnData() || (zCand_ && zCand_.flavor() == PDGMUON)){
    const pat::TriggerPathRefVector acceptedPaths = triggerEvent_.acceptedPaths();
    if(debugme) cout<<"Using "<<acceptedPaths.size()<<" accepted paths from HLT"<<endl;
    for (size_t i = 0; i < acceptedPaths.size(); i++){
      string A = acceptedPaths[i]->name();
      for (size_t j = 0; j < triggersToUse_.size(); j++){
        if(SameTrigger(A, triggersToUse_[j])){
          if(debugme) cout<<"Match A: "<<acceptedPaths[i]->name()<<" B: "<<triggersToUse_[j]<<endl;
          if(acceptedPaths[i]->prescale() == 1  && acceptedPaths[i]->wasAccept()) return true;
          break;
        }
      }
    }
  }else{
    /*
    const pat::TriggerPathCollection* allPaths = triggerEvent_.paths();
    //const pat::TriggerPathRefVector allPaths = triggerEvent_.pathRefs();
    //for (size_t i = 0; i < allPaths->size(); i++){
    for ( pat::TriggerPathCollection::const_iterator iPath = allPaths->begin(); iPath != allPaths->end(); ++iPath ) {
//      string A = allPaths[i]->name();
      cout<<"Name: "<<iPath->name()<<endl;
      
      //if(A.find("HLT_Mu"))
      
    } 
    */   
    return true;//Cory: This is not good, but will pass HLT in the meantime.
  }
  return false;
}//--- PassTriggersCut()

bool
WZAnalyzer::EMuOverlap(const pat::Electron & e, 
                       const vector<pat::Muon    > & ms){
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
    heep::Ele& e1 = FindElectron(*zCand_.daughter(0));
    heep::Ele& e2 = FindElectron(*zCand_.daughter(1));
    return (PassTriggerEmulation(e1) && PassTriggerEmulation(e2));
  }else if(zCand_.flavor() == PDGMUON){
    TeVMuon& m1 = FindMuon(*zCand_.daughter(0));
    TeVMuon& m2 = FindMuon(*zCand_.daughter(1));
    return PassTriggerMatch(m1, m2); 
  }
  return false;
}

inline bool
WZAnalyzer::PassTriggerMatch(const heep::Ele& e1, const heep::Ele& e2){
  return (e1.patEle().triggerObjectMatches().size() > 0 &&
          e2.patEle().triggerObjectMatches().size() > 0 &&
          max(e1.patEle().triggerObjectMatches()[0].et(), e2.patEle().triggerObjectMatches()[0].et()) > 17. &&
          min(e1.patEle().triggerObjectMatches()[0].et(), e2.patEle().triggerObjectMatches()[0].et()) > 8.);
}

inline bool
WZAnalyzer::PassTriggerMatch(const TeVMuon& m1, const TeVMuon& m2){
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
bool WZAnalyzer::PassTriggerEmulation(const heep::Ele& elec){
  if(!elec.isEcalDriven()) return false;
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
inline float WZAnalyzer::Calc_Ht(){
  return WLepPt() + ZLepPt(0) + ZLepPt(1);
}//--- CalcHt

inline float WZAnalyzer::Calc_Q(){
  return wzCand_.mass("minPz") - zCand_.mass() - WMASS;
}

inline int WZAnalyzer::Calc_EvtType(){
  return (zCand_ && wCand_) ?  2 * (zCand_.flavor() != 11) + (wCand_.flavor() != 11) : -999;
}

inline bool WZAnalyzer::inEE(const TeVMuon& mu){
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
  Ht_= -999;
  Q_ = -999;
  TT = TF = false;
  PU_NumInteractions_ = -999.;
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
  printSummary(fi->samplename, out);  
}

void WZAnalyzer::endAnalysis(ofstream & out){
}

heep::Ele &
WZAnalyzer::FindElectron(reco::Candidate & p){
  for(uint i=0; i<electrons_.size(); ++i){
    if(Match(electrons_[i], p)) return electrons_[i];
  }
  cout<<"Didn't find match for electron, returning random one!!!\n";
  return electrons_[0];
}

TeVMuon &
WZAnalyzer::FindMuon(reco::Candidate & p){
  for(uint i=0; i<muons_.size(); ++i){
    if(Match(muons_[i], p)) return muons_[i];
  }
  cout<<"Didn't find match for muon!!!, returning random one\n";
  return muons_[0];
}

bool
WZAnalyzer::Match(heep::Ele & p1, reco::Candidate & p2){
  float tolerance = 0.0001;
  if (p1.patEle().pdgId() == p2.pdgId() &&
      fabs(p1.eta() - p2.eta()) < tolerance &&
      fabs(p1.phi() - p2.phi()) < tolerance
    )
    return true;
  return false;
}

bool
WZAnalyzer::Match(TeVMuon & p1, reco::Candidate & p2){
  float tolerance = 0.0001;
  if (p1.pdgId() == p2.pdgId() &&
      fabs(p1.eta() - p2.eta()) < tolerance &&
      fabs(p1.phi() - p2.phi()) < tolerance
    )
    return true;
  return false;
}

float
WZAnalyzer::WLepPt(){
  if(wCand_.flavor() == PDGELEC){
    return FindElectron(*wCand_.daughter(0)).patEle().pt();
  }else if(wCand_.flavor() == PDGMUON){
    return FindMuon(*wCand_.daughter(0)).pt();
  }
  return -999.;
}

inline float
WZAnalyzer::ZLepPt(int idx){
  if(zCand_.flavor() == PDGELEC)
    return FindElectron(*zCand_.daughter(idx)).patEle().pt();
  else if(zCand_.flavor() == PDGMUON)
    return FindMuon(*zCand_.daughter(idx)).pt();
  return -999.;
}

inline float
  WZAnalyzer::ElecPU(const heep::Ele & e){
  return rhoFastJet_*effectiveElecArea_[e.patEle().isEE()];
}

inline float
  WZAnalyzer::MuonPU(const TeVMuon & m){
  return rhoFastJet_*effectiveMuonArea_[inEE(m)];
}
