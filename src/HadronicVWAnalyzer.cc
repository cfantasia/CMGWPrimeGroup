#include "UserCode/CMGWPrimeGroup/interface/HadronicVWAnalyzer.h"

using namespace std;

HadronicVWAnalyzer::HadronicVWAnalyzer(){}
HadronicVWAnalyzer::HadronicVWAnalyzer(const edm::ParameterSet & cfg, WPrimeUtil * wprimeUtil) :
  AnalyzerBase(cfg, wprimeUtil){
  FillCutFns();
  if(debugme) printf("Using %i cuts\n",NCuts_);

  vertexLabel_ = cfg.getParameter<edm::InputTag>("vertexTag");

  effectiveElecArea_ = cfg.getParameter<vector<double> >("effectiveElecArea");
  effectiveMuonArea_ = cfg.getParameter<vector<double> >("effectiveMuonArea");
 
// +++++++++++++++++++General Cut values

// +++++++++++++++++++W Cuts

// +++++++++++++++++++V Cuts

  ClearEvtVariables();
}

HadronicVWAnalyzer::~HadronicVWAnalyzer(){
  outCandEvt_.close();
}

void HadronicVWAnalyzer::FillCutFns(){
  mFnPtrs_["NoCuts"] = &HadronicVWAnalyzer::PassNoCut;
  mFnPtrs_["HLT"] = &HadronicVWAnalyzer::PassTriggersCut;
  mFnPtrs_["MinNLeptons"] = &HadronicVWAnalyzer::PassMinNLeptonsCut;
  mFnPtrs_["MaxNLeptons"] = &HadronicVWAnalyzer::PassMaxNLeptonsCut;
  mFnPtrs_["MinNJets"] = &HadronicVWAnalyzer::PassMinNJetsCut;
  mFnPtrs_["ValidW"] = &HadronicVWAnalyzer::PassValidWCut;
  mFnPtrs_["ValidV"] = &HadronicVWAnalyzer::PassValidVCut;
  mFnPtrs_["ValidWVCand"] = &HadronicVWAnalyzer::PassValidWVCut;
  mFnPtrs_["VMass"] = &HadronicVWAnalyzer::PassVMassCut;
  mFnPtrs_["WTransMass"] = &HadronicVWAnalyzer::PassWtransMassCut;
  mFnPtrs_["MET"] = &HadronicVWAnalyzer::PassMinMETCut;
  mFnPtrs_["Vpt"] = &HadronicVWAnalyzer::PassVptCut;
  mFnPtrs_["Wpt"] = &HadronicVWAnalyzer::PassWptCut;
  mFnPtrs_["AllCuts"] = &HadronicVWAnalyzer::PassNoCut;

  CutFns_.resize(NCuts_);
  for(int i=0; i<NCuts_; ++i){
    if(mFnPtrs_.find(Cuts_[i]) == mFnPtrs_.end()){
      cout<<"Didn't find Cut named "<<Cuts_[i]<<endl;
      abort();
    }
    CutFns_[i] = mFnPtrs_.find(Cuts_[i])->second;
  }
  
}

//--------------------------------------------------------------
void HadronicVWAnalyzer::Declare_Histos(const TFileDirectory & dir)
{
  if(debugme) printf("Declare histos\n");

  DeclareHistoSet("hWVMass", "Reconstructed WV Invariant Mass",
                  "M_{WV} (GeV)", 1200, 0, 1200, "GeV", hWVMass,dir);
  DeclareHistoSet("hWV3e0muMass", "Reconstructed WV(3e0#mu) Invariant Mass",
                  "M_{WV}^{3e0#mu} (GeV)", 1200, 0, 1200, "GeV", hWV3e0muMass,dir);
  DeclareHistoSet("hWV2e1muMass", "Reconstructed WV(2e1#mu) Invariant Mass",
                  "M_{WV}^{2e1#mu} (GeV)", 1200, 0, 1200, "GeV", hWV2e1muMass,dir);
  DeclareHistoSet("hWV1e2muMass", "Reconstructed WV(1e2#mu) Invariant Mass",
                  "M_{WV}^{1e2#mu} (GeV)", 1200, 0, 1200, "GeV", hWV1e2muMass,dir);
  DeclareHistoSet("hWV0e3muMass", "Reconstructed WV(0e3#mu) Invariant Mass",
                  "M_{WV}^{0e3#mu} (GeV)", 1200, 0, 1200, "GeV", hWV0e3muMass,dir);

//Q=M_{WV} - M_W - M_V
  DeclareHistoSet("hQ", "Q=M_{WV} - M_{W} - M_{V}",
                  "Q (GeV)", 50, 0, 500, "GeV", hQ,dir);
  DeclareHistoSet("hWVTransMass", "Reconstructed WV Transverse Mass",
                  "M_{WV}^{T} (GeV)", 100, 0, 1000, "GeV", hWVTransMass,dir);
//WVpt Histos
  DeclareHistoSet("hWVpt", "Reconstructed WV Transverse Momentum",
                  "p_{WV}^{T} (GeV)", 50, 0, 500, "GeV", hWVpt,dir);

  DeclareHistoSet("hEvtType", "Event Type",
                  "N_{#mu}", 4, 0, 4, "NONE", hEvtType,dir);
  DeclareHistoSet("hEvtTypeP", "Event Type for Q=+1",
                  "N_{#mu},W^{+}", 4, 0, 4, "NONE", hEvtTypeP,dir);
  DeclareHistoSet("hEvtTypeM", "Event Type for Q=-1",
                  "N_{#mu},W^{-}", 4, 0, 4, "NONE", hEvtTypeM,dir);

///////////////////////////
//V Mass Histos
  DeclareHistoSet("hVMass" , "Reconstructed Mass of V",
                  "M_{V} (GeV)", 30, 60, 120, "GeV", hVMass,dir);
  DeclareHistoSet("hVeeMass","Reconstructed Mass of Vee",
                  "M_{V}^{ee} (GeV)", 30, 60, 120, "GeV", hVeeMass,dir);
  DeclareHistoSet("hVmmMass","Reconstructed Mass of V#mu#mu",
                  "M_{V}^{#mu#mu} (GeV)", 30, 60, 120, "GeV", hVmmMass,dir);
  DeclareHistoSet("hV3e0muMass" , "Reconstructed Mass of V(3e0#mu)",
                  "M_{V}^{3e0#mu} (GeV)", 30, 60, 120, "GeV", hV3e0muMass,dir);
  DeclareHistoSet("hV2e1muMass" , "Reconstructed Mass of V(2e1#mu)",
                  "M_{V}^{2e1#mu} (GeV)", 30, 60, 120, "GeV", hV2e1muMass,dir);
  DeclareHistoSet("hV1e2muMass" , "Reconstructed Mass of V(1e2#mu)",
                  "M_{V}^{1e2#mu} (GeV)", 30, 60, 120, "GeV", hV1e2muMass,dir);
  DeclareHistoSet("hV0e3muMass" , "Reconstructed Mass of V(0e3#mu)",
                  "M_{V}^{0e3#mu} (GeV)", 30, 60, 120, "GeV", hV0e3muMass,dir);

//Vpt Histos
  DeclareHistoSet("hVpt", "p_{T}^{V}", 
                  "p_{T}^{V} (GeV)", 40, 0, 400, "GeV", hVpt,dir);
  DeclareHistoSet("hVeept", "p_{T}^{V#rightarrowee}", 
                  "p_{T}^{V#rightarrowee} (GeV)", 40, 0, 400, "GeV", hVeept,dir);
  DeclareHistoSet("hVmmpt", "p_{T}^{V#rightarrow#mu#mu}", 
                  "p_{T}^{V#rightarrow#mu#mu} (GeV)", 40, 0, 400, "GeV", hVmmpt,dir);
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
  DeclareHistoSet("hWptVee", "p_{T}^{W,V#rightarrowee}", 
                  "p_{T}^{W,V#rightarrowee} (GeV)", 40, 0, 400, "GeV", hWptVee,dir);
  DeclareHistoSet("hWptVmm", "p_{T}^{W,V#rightarrow#mu#mu}", 
                  "p_{T}^{W,V#rightarrow#mu#mu} (GeV)", 40, 0, 400, "GeV", hWptVmm,dir);

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
  DeclareHistoSet("hNLLepsVee", "Number of Loose Leptons in Event, V#rightarrowee",
                  "N_{l}^{Loose,V#rightarrowee}", 10, 0, 10, "NONE", hNLLepsVee,dir);
  DeclareHistoSet("hNLLepsVmm", "Number of Loose Leptons in Event",
                  "N_{l}^{Loose,V#rightarrow#mu#mu}", 10, 0, 10, "NONE", hNLLepsVmm,dir);

  DeclareHistoSet("hNTElec", "Number of Tight Electrons in Event",
                  "N_{e}", 10, 0, 10, "NONE", hNTElec,dir);
  DeclareHistoSet("hNTMuon", "Number of Tight Muons in Event",
                  "N_{#mu}", 10, 0, 10, "NONE", hNTMuon,dir);
  DeclareHistoSet("hNTLeps", "Number of Tight Leptons in Event",
                  "N_{l}", 10, 0, 10, "NONE", hNTLeps,dir);

  DeclareHistoSet("hNJets", "Number of Jets in Event",
                  "N_{Jets}", 10, 0, 10, "NONE", hNJets,dir);
  DeclareHistoSet("hNJetsVee", "Number of Jets in Event, V#rightarrowee",
                  "N_{Jets}^{V#rightarrowee}", 10, 0, 10, "NONE", hNJetsVee,dir);
  DeclareHistoSet("hNJetsVmm", "Number of Jets in Event, V#rightarrow#mu#mu",
                  "N_{Jets}^{V#rightarrow#mu#mu}", 10, 0, 10, "NONE", hNJetsVmm,dir);

  DeclareHistoSet("hNVtxs", "Number of Vertexs in Event",
                  "N_{Vtx}", 50, 0, 50, "NONE", hNVtxs,dir);
  DeclareHistoSet("hNVtxsVee", "Number of Vertexs in Event, V#rightarrowee",
                  "N_{Vtx}^{V#rightarrowee}", 50, 0, 50, "NONE", hNVtxsVee,dir);
  DeclareHistoSet("hNVtxsVmm", "Number of Vertexs in Event, V#rightarrow#mu#mu",
                  "N_{Vtx}^{V#rightarrow#mu#mu}", 50, 0, 50, "NONE", hNVtxsVmm,dir);

  DeclareHistoSet("hWenuCombRelIso", "Comb Rel Iso of W Electron",
                  "Electron Combined Relative Isolation", 20, 0, 0.2, "NONE", hWenuCombRelIso,dir);
  DeclareHistoSet("hWmnuCombRelIso", "Comb Rel Iso of W Muon",
                  "Muon Combined Relative Isolation", 20, 0, 0.2, "NONE", hWmnuCombRelIso,dir);
  

  tWVCand = dir.make<TTree>("tWVCand", "Analysis Variables after WVCand");//Only 1 for now;
  tWVCand->Branch("WVMass", &WVMass_);
  tWVCand->Branch("EvtType", &evtType_);
  tWVCand->Branch("Vpt", &Vpt_);
  tWVCand->Branch("Wpt", &Wpt_);
  tWVCand->Branch("weight", &weight_);

///Eff Plots///////
  string title = Form("Expected # of Events / %.0f pb^{-1}",  wprimeUtil_->getLumi_ipb());
  title = title + ";;" + title;
  hNumEvts = NULL; hNumEvts = dir.make<TH1F>("hNumEvts",title.c_str(),NCuts_,0,NCuts_);
  for(int i=0; i<NCuts_; ++i) hNumEvts->GetXaxis()->SetBinLabel(i+1,Cuts_[i].c_str());

}//Declare_Histos

//Fill Histograms
//-----------------------------------------------------------
void HadronicVWAnalyzer::Fill_Histos(const int& index, const float& weight)
{
//-----------------------------------------------------------
  if(debugme) printf("Filling Histos\n");
  if(wCand_ && vCand_){
    hWVMass[index]->Fill(wvCand_.mass("minPz"), weight);
    if     (evtType_ == 0) hWV3e0muMass[index]->Fill(wvCand_.mass("minPz"), weight);
    else if(evtType_ == 1) hWV2e1muMass[index]->Fill(wvCand_.mass("minPz"), weight);
    else if(evtType_ == 2) hWV1e2muMass[index]->Fill(wvCand_.mass("minPz"), weight);
    else if(evtType_ == 3) hWV0e3muMass[index]->Fill(wvCand_.mass("minPz"), weight);
    hQ[index]->Fill(Q_, weight); 
    hWVTransMass[index]->Fill(wvCand_.transMass(), weight);
    hWVpt[index]->Fill(wvCand_.pt(), weight);
    hEvtType[index]->Fill(evtType_, weight);
    if     (wCand_.charge() > 0) hEvtTypeP[index]->Fill(evtType_, weight);
    else if(wCand_.charge() < 0) hEvtTypeM[index]->Fill(evtType_, weight);
    if     (vCand_.flavor() == PDGELEC){
      hWptVee[index]->Fill(wCand_.pt(), weight);
      hNLLepsVee[index]->Fill(looseElectrons_.size()+looseMuons_.size(), weight);
      hNJetsVee[index]->Fill(looseJets_.size(), weight);
      hNVtxsVee[index]->Fill(vertices_.size(), weight);
    }else if(vCand_.flavor() == PDGMUON){ 
      hWptVmm[index]->Fill(wCand_.pt(), weight);
      hNLLepsVmm[index]->Fill(looseElectrons_.size()+looseMuons_.size(), weight);
      hNJetsVmm[index]->Fill(looseJets_.size(), weight);
      hNVtxsVmm[index]->Fill(vertices_.size(), weight);
    }
    if(Cuts_[index] == "ValidWVCand"){
      tWVCand->Fill();
    }
  }
  if(vCand_){
    hVMass[index]->Fill(vCand_.mass(), weight);
    hVpt[index]->Fill(vCand_.pt(), weight);
    if      (vCand_.flavor() == PDGELEC){
      hVeeMass[index]->Fill(vCand_.mass(), weight);
      hVeept[index]->Fill(vCand_.pt(), weight);
      hMETee[index]->Fill(met_.et(), weight);
    }else if (vCand_.flavor() == PDGMUON){
      hVmmMass[index]->Fill(vCand_.mass(), weight);
      hMETmm[index]->Fill(met_.et(), weight);
      hVmmpt[index]->Fill(vCand_.pt(), weight);
    }
    if     (evtType_ == 0) hV3e0muMass[index]->Fill(vCand_.mass(), weight);
    else if(evtType_ == 1) hV2e1muMass[index]->Fill(vCand_.mass(), weight);
    else if(evtType_ == 2) hV1e2muMass[index]->Fill(vCand_.mass(), weight);
    else if(evtType_ == 3) hV0e3muMass[index]->Fill(vCand_.mass(), weight);

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

  hNJets[index]->Fill(looseJets_.size(), weight);
  hNVtxs[index]->Fill(vertices_.size(), weight);

}//Fill_Histos

inline void
HadronicVWAnalyzer::CalcVVariables(){
  if (debugme) cout<<"In Calc V Variables\n";
  vCand_ = getWCand(looseJets_);
  Vpt_ = vCand_.pt();
  if(debugme) printf("    Contains: %i tight V candidate(s)\n", (bool)vCand_);
  if(debugme){
    PrintEventLeptons(); 
    PrintEventDetails();
  }
}

inline void
HadronicVWAnalyzer::CalcWVariables(){
  if (debugme) cout<<"In Calc W Variables\n";
  wCand_ = getWCand(tightElectrons_, tightMuons_, met_);
  Wpt_ = wCand_.pt();
  if(debugme) printf("    Contains: %i tight W candidate(s)\n", (bool)wCand_);
  if(debugme){
    PrintEventLeptons(); 
    PrintEventDetails();
  }
}

inline void
HadronicVWAnalyzer::CalcWElecVariables(){
  if (debugme) cout<<"In Calc W Elec Variables\n";
  wCand_ = getWCand(tightElectrons_, met_);
  if(debugme) printf("    Contains: %i tight W candidate(s)\n", (bool)wCand_);
}

inline void
HadronicVWAnalyzer::CalcWMuonVariables(){
  if (debugme) cout<<"In Calc W Muon Variables\n";
  wCand_ = getWCand(tightMuons_, met_);
  if(debugme) printf("    Contains: %i tight W candidate(s)\n", (bool)wCand_);
}

inline void
HadronicVWAnalyzer::CalcWVVariables(){
  if (debugme) cout<<"In Calc WV Variables\n";
  wvCand_ = (vCand_ && wCand_) ? WVCandidate(vCand_, wCand_) : WVCandidate();
  WVMass_ = wvCand_.mass("minPz");
  Q_ = (vCand_ && wCand_) ? Calc_Q() : -999.;
  if(debugme) PrintEventDetails();
}

void
HadronicVWAnalyzer::CalcEventVariables(){
  if (debugme) cout<<"In Calc Event Variables\n";
  evtType_ = (vCand_ && wCand_) ? Calc_EvtType() : -999;
  if(debugme) printf("evt Type: %i, V Flav: %i, W Flav: %i\n", evtType_, (int)vCand_.flavor(), (int)wCand_.flavor());
}

void 
HadronicVWAnalyzer::eventLoop(edm::EventBase const & event){
  ClearEvtVariables();
  if(debugme) WPrimeUtil::PrintEvent(event);

  // Preselection - skip events that don't look promising
  if (doPreselect_){
    if(debugme) cout<<"Testing Preselection...\n";
  }
//  // Get information about flavorHistory
//  const uint flavorType = getProduct<uint>(event, "flavorHistoryFilter", 0);
//  if(debugme) printf("    FlavorType: %i", flavorType);

  //Get Jets
  event.getByLabel(jetsLabel_,patJetsH_);
  if(debugme) printf("    Contains: %i pat jets(s)\n",
                     (int)patJetsH_->size());
  for (size_t i = 0; i < patJetsH_->size(); i++) {
    const pat::Jet & jet = (*patJetsH_.product())[i];
    if (looseJet_(jet, jetResult_))
      looseJets_.push_back(jet);
  }
  if(looseJets_.size() == 0) return;


  // Get leptons
  ////////Cory:Not using electrons yet
  //event.getByLabel(electronsLabel_,patElectronsH_);
  event.getByLabel(muonsLabel_,patMuonsH_);
  //Cory: Didn't save pfMET
  //event.getByLabel(metLabel_, metH_);
  if(useAdjustedMET_) event.getByLabel(pfCandsLabel_, pfCandidatesH_);
  WPrimeUtil::getMuonsMET(patMuonsH_, muonAlgo_, muons_,
                          metH_, useAdjustedMET_, met_,
                          pfCandidatesH_);
  
/*
  WPrimeUtil::getLeptonsMET(patElectronsH_, electrons_,
                            patMuonsH_, muonAlgo_, muons_,
                            metH_, useAdjustedMET_, met_,
                            pfCandidatesH_);
  if(debugme) printf("    Contains: %i electron(s), %i muon(s)\n",
                          (int)electrons_.size(), (int)muons_.size());

  rhoFastJet_ = getProduct<double>(event,"kt6PFJets:rho");


  // Make vectors of leptons passing various criteria
  for (size_t i = 0; i < electrons_.size(); i++) {
    if(Overlap(electrons_[i].patEle(), muons_)) continue;
    const float pu = ElecPU(electrons_[i]);
    if (looseElectron_(electrons_[i].patEle(), electronResult_, pu))
      looseElectrons_.push_back(electrons_[i]);

    if (tightElectron_(electrons_[i].patEle(), electronResult_, pu))
      tightElectrons_.push_back(electrons_[i]);
  }
*/
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
    printf("    Contains: %i loose electron(s), %i loose muon(s), %i loose jet(s)\n",
           (int)looseElectrons_.size(), (int)looseMuons_.size(), (int)looseJets_.size());
    printf("    Contains: %i tight electron(s), %i tightmuon(s)\n",
           (int)tightElectrons_.size(), (int)tightMuons_.size());
  }

  ///////////////////
  //Cory: Remove overlap btw jets and leptons! 
  //(maybe I don't have to for pfJets)
  ///////////////////


  //Get Trigger 
  triggerEvent_ = getProduct<pat::TriggerEvent>(event,hltEventLabel_); 

  //Get Vertex
  vertices_ = getProduct<vector<reco::Vertex> >(event,vertexLabel_);

  float PU_Weight = 1.;
  if(!wprimeUtil_->runningOnData()){//Don't do this for data
    if(debugme){
      //Cory: Update this for VW!
      GenParticleV genParticles = getProduct<GenParticleV>(event, "genParticles");
      const reco::Candidate * genV = 0;
      const reco::Candidate * genW = 0;
      for (size_t i = 0; i < genParticles.size(); i++){
        if (abs(genParticles[i].pdgId()) == 23){
          genV = & genParticles[i];
          cout<<"Mass of gen V is "<<genV->mass()<<endl;
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

void HadronicVWAnalyzer::PrintDebugEvent() const{
  WPrimeUtil::PrintPassingTriggers(triggerEvent_,triggersToUse_);
  PrintEventDetails();
  PrintEventLeptons();
  PrintLeptons();
}

void HadronicVWAnalyzer::PrintEventDetails() const{
  if(vCand_){
    cout<<" V Flavor: "<<vCand_.flavor()
        <<" V Mass: "<<vCand_.mass()
//        <<" V lep1 pt "<<VLepPt(0)
//        <<" V lep2 pt "<<VLepPt(1)
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
  if(vCand_ && wCand_ && wvCand_.mass("minPz")>0.){
    cout<<" WV Mass: "<<wvCand_.mass("minPz")
        <<" Neu Pz: "<<wvCand_.neutrinoPz("minPz")
        <<" Vpt: "<<vCand_.pt()
        <<" Wpt: "<<wCand_.pt()
        <<endl;
  }
  return;
}

void
HadronicVWAnalyzer::PrintEventLeptons() const{
  if     (vCand_.flavor() == PDGELEC){
//    PrintElectron(*vCand_.elec1(), PDGZ);
//    PrintElectron(*vCand_.elec2(), PDGZ);
    PrintElectron(WPrimeUtil::Find(*vCand_.daughter(0), electrons_));
    PrintElectron(WPrimeUtil::Find(*vCand_.daughter(1), electrons_));
  }else if(vCand_.flavor() == PDGMUON){
    PrintMuon(WPrimeUtil::Find(*vCand_.daughter(0), muons_));
    PrintMuon(WPrimeUtil::Find(*vCand_.daughter(1), muons_));
//    PrintMuon(*vCand_.muon1(), PDGZ);
//    PrintMuon(*vCand_.muon2(), PDGZ);
  }

  if     (wCand_.flavor() == PDGELEC){   
    PrintElectron(WPrimeUtil::Find(*wCand_.daughter(0), electrons_));
//    PrintElectron(*wCand_.elec(), PDGW);
  }else if(wCand_.flavor() == PDGMUON){
    PrintMuon(WPrimeUtil::Find(*wCand_.daughter(0), muons_));
//    PrintMuon    (*wCand_.muon(), PDGW);
  }
}

/////////////////Accessors///////////////////////

/////////////////Modifiers///////////////////////

/////////////////Cuts///////////////////////
bool
HadronicVWAnalyzer::PassCuts(const float& weight){
  if (debugme) cout<<"In Pass Cuts\n";
  
  for(int i=0; i<NCuts_; ++i){
    if(Cuts_[i] == "ValidV"){
      CalcVVariables();
    }else if(Cuts_[i] == "ValidW"){
      CalcWVariables();
      CalcEventVariables();
    }else if(Cuts_[i] == "ValidWVCand"){
      CalcWVVariables();
    }

    if(!(this->*CutFns_[i])()) return false;
    Tabulate_Me(i,weight); 
  }
  return true;
}

//Trigger requirements
//-----------------------------------------------------------
bool HadronicVWAnalyzer::PassTriggersCut() const{
  if(debugme) cout<<"Trigger requirements"<<endl;
  //Apply the trigger if running on data or MC 
  //If MC, apply if no V or if V exists, zCand == PDGMuon)
  if(wprimeUtil_->runningOnData() || !vCand_ || vCand_.flavor() == PDGMUON){
    return WPrimeUtil::PassTriggersCut(triggerEvent_,triggersToUse_);
  }else{
    return true;//Cory: This is not good, but will pass HLT in the meantime.
  }
  return false;
}//--- PassTriggersCut()


inline bool HadronicVWAnalyzer::PassValidWVCut() const{
  return wvCand_ && wvCand_.mass("minPz")>0.;
}

////////////////////////////////
/////////Check Electron Properties/////
////////////////////////////////
bool HadronicVWAnalyzer::PassTriggerEmulation(const heep::Ele& elec, const float minPt) const{
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

///////////////////////////////////

inline float HadronicVWAnalyzer::Calc_Q() const{
  return wvCand_.mass("minPz") - vCand_.mass() - WMASS;
}

inline int HadronicVWAnalyzer::Calc_EvtType() const{
  return (vCand_ && wCand_) ?  2 * (vCand_.flavor() != 11) + (wCand_.flavor() != 11) : -999;
}

inline bool HadronicVWAnalyzer::inEE(const TeVMuon& mu) const{
  return fabs(mu.eta()) >= 1.05;
}

///////////////Utilities//////////////////
//--------------------------------------------------------------

inline void
HadronicVWAnalyzer::ClearEvtVariables(){
  looseJets_.clear();
  electrons_.clear();
  looseElectrons_.clear();
  tightElectrons_.clear();
  muons_.clear();
  looseMuons_.clear();
  tightMuons_.clear();
  met_ = pat::MET();
  vCand_ = WCandidate();
  wCand_ = WCandidate();
  wvCand_ = WVCandidate();
  evtType_ = -999;
  WVMass_ = -999;
  Vpt_ = -999;
  Wpt_ = -999;
  Q_ = -999;
  weight_ = 0;
  rhoFastJet_ = 0;
}

float
HadronicVWAnalyzer::WLepPt() const{
  if(wCand_.flavor() == PDGELEC){
    return WPrimeUtil::Find(*wCand_.daughter(0), electrons_).patEle().pt();
  }else if(wCand_.flavor() == PDGMUON){
    return WPrimeUtil::Find(*wCand_.daughter(0), muons_).pt();
  }
  return -999.;
}

inline float
  HadronicVWAnalyzer::ElecPU(const heep::Ele & e) const{
  return rhoFastJet_*effectiveElecArea_[e.patEle().isEE()];
}

inline float
  HadronicVWAnalyzer::MuonPU(const TeVMuon & m) const{
  return rhoFastJet_*effectiveMuonArea_[inEE(m)];
}
