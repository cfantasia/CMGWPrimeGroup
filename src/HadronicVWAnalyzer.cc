#include "UserCode/CMGWPrimeGroup/interface/HadronicVWAnalyzer.h"

using namespace std;

HadronicVWAnalyzer::HadronicVWAnalyzer(){}
HadronicVWAnalyzer::HadronicVWAnalyzer(const edm::ParameterSet & cfg, int fileToRun) :
  AnalyzerBase(cfg, fileToRun){
  setupCutOrder();
  if(debugme) printf("Using %i cuts\n",NCuts_);

  effectiveElecArea_ = cfg.getParameter<vector<double> >("effectiveElecArea");
  effectiveMuonArea_ = cfg.getParameter<vector<double> >("effectiveMuonArea");
 
// +++++++++++++++++++General Cut values

// +++++++++++++++++++W Cuts

// +++++++++++++++++++V Cuts
}

HadronicVWAnalyzer::~HadronicVWAnalyzer(){
}

void HadronicVWAnalyzer::setupCutOrder(){
  mFnPtrs_["NoCuts"] = boost::bind(&HadronicVWAnalyzer::passNoCut, this);
  mFnPtrs_["HLT"] = boost::bind(&HadronicVWAnalyzer::passTriggersCut, this);
  mFnPtrs_["MinNLeptons"] = boost::bind(&HadronicVWAnalyzer::passMinNLeptonsCut, this);
  mFnPtrs_["MaxNLeptons"] = boost::bind(&HadronicVWAnalyzer::passMaxNLeptonsCut, this);
  mFnPtrs_["MinNJets"] = boost::bind(&HadronicVWAnalyzer::passMinNJetsCut, this);
  mFnPtrs_["ValidW"] = boost::bind(&HadronicVWAnalyzer::passValidWCut, this);
  mFnPtrs_["ValidV"] = boost::bind(&HadronicVWAnalyzer::passValidVCut, this);
  mFnPtrs_["ValidVWCand"] = boost::bind(&HadronicVWAnalyzer::passValidVWCut, this);
  mFnPtrs_["VMass"] = boost::bind(&HadronicVWAnalyzer::passVMassCut, this);
  mFnPtrs_["WTransMass"] = boost::bind(&HadronicVWAnalyzer::passWtransMassCut, this);
  mFnPtrs_["MET"] = boost::bind(&HadronicVWAnalyzer::passMinMETCut, this);
  mFnPtrs_["Vpt"] = boost::bind(&HadronicVWAnalyzer::passVptCut, this);
  mFnPtrs_["Wpt"] = boost::bind(&HadronicVWAnalyzer::passWptCut, this);
  mFnPtrs_["AllCuts"] = boost::bind(&HadronicVWAnalyzer::passNoCut, this);

  fillCuts();
}

//--------------------------------------------------------------
void HadronicVWAnalyzer::defineHistos(const TFileDirectory & dir){
  if(debugme) printf("Declare histos\n");
  AnalyzerBase::defineHistos(dir);

  defineHistoset("hVWMass", "Reconstructed VW Invariant Mass",
                  "M_{VW} (GeV)", 1200, 0, 1200, "GeV", hVWMass,dir);
  defineHistoset("hVW3e0muMass", "Reconstructed VW(3e0#mu) Invariant Mass",
                  "M_{VW}^{3e0#mu} (GeV)", 1200, 0, 1200, "GeV", hVW3e0muMass,dir);
  defineHistoset("hVW2e1muMass", "Reconstructed VW(2e1#mu) Invariant Mass",
                  "M_{VW}^{2e1#mu} (GeV)", 1200, 0, 1200, "GeV", hVW2e1muMass,dir);
  defineHistoset("hVW1e2muMass", "Reconstructed VW(1e2#mu) Invariant Mass",
                  "M_{VW}^{1e2#mu} (GeV)", 1200, 0, 1200, "GeV", hVW1e2muMass,dir);
  defineHistoset("hVW0e3muMass", "Reconstructed VW(0e3#mu) Invariant Mass",
                  "M_{VW}^{0e3#mu} (GeV)", 1200, 0, 1200, "GeV", hVW0e3muMass,dir);

//Q=M_{VW} - M_W - M_V
  defineHistoset("hQ", "Q=M_{VW} - M_{W} - M_{V}",
                  "Q (GeV)", 50, 0, 500, "GeV", hQ,dir);
  defineHistoset("hVWTransMass", "Reconstructed VW Transverse Mass",
                  "M_{VW}^{T} (GeV)", 100, 0, 1000, "GeV", hVWTransMass,dir);
//VWpt Histos
  defineHistoset("hVWpt", "Reconstructed VW Transverse Momentum",
                  "p_{VW}^{T} (GeV)", 50, 0, 500, "GeV", hVWpt,dir);

  defineHistoset("hEvtType", "Event Type",
                  "N_{#mu}", 4, 0, 4, "NONE", hEvtType,dir);
  defineHistoset("hEvtTypeP", "Event Type for Q=+1",
                  "N_{#mu},W^{+}", 4, 0, 4, "NONE", hEvtTypeP,dir);
  defineHistoset("hEvtTypeM", "Event Type for Q=-1",
                  "N_{#mu},W^{-}", 4, 0, 4, "NONE", hEvtTypeM,dir);

///////////////////////////
//V Mass Histos
  defineHistoset("hVMass" , "Reconstructed Mass of V",
                  "M_{V} (GeV)", 30, 60, 120, "GeV", hVMass,dir);
  defineHistoset("hVeeMass","Reconstructed Mass of Vee",
                  "M_{V}^{ee} (GeV)", 30, 60, 120, "GeV", hVeeMass,dir);
  defineHistoset("hVmmMass","Reconstructed Mass of V#mu#mu",
                  "M_{V}^{#mu#mu} (GeV)", 30, 60, 120, "GeV", hVmmMass,dir);
  defineHistoset("hV3e0muMass" , "Reconstructed Mass of V(3e0#mu)",
                  "M_{V}^{3e0#mu} (GeV)", 30, 60, 120, "GeV", hV3e0muMass,dir);
  defineHistoset("hV2e1muMass" , "Reconstructed Mass of V(2e1#mu)",
                  "M_{V}^{2e1#mu} (GeV)", 30, 60, 120, "GeV", hV2e1muMass,dir);
  defineHistoset("hV1e2muMass" , "Reconstructed Mass of V(1e2#mu)",
                  "M_{V}^{1e2#mu} (GeV)", 30, 60, 120, "GeV", hV1e2muMass,dir);
  defineHistoset("hV0e3muMass" , "Reconstructed Mass of V(0e3#mu)",
                  "M_{V}^{0e3#mu} (GeV)", 30, 60, 120, "GeV", hV0e3muMass,dir);

//Vpt Histos
  defineHistoset("hVpt", "p_{T}^{V}", 
                  "p_{T}^{V} (GeV)", 40, 0, 400, "GeV", hVpt,dir);
  defineHistoset("hVeept", "p_{T}^{V#rightarrowee}", 
                  "p_{T}^{V#rightarrowee} (GeV)", 40, 0, 400, "GeV", hVeept,dir);
  defineHistoset("hVmmpt", "p_{T}^{V#rightarrow#mu#mu}", 
                  "p_{T}^{V#rightarrow#mu#mu} (GeV)", 40, 0, 400, "GeV", hVmmpt,dir);
//MET Histos
  defineHistoset("hMET", "MET",
                  "#slash{E}_{T} (GeV)", 30, 0, 300, "GeV", hMET,dir);
  defineHistoset("hMETee", "MET",
                  "#slash{E}_{T}^{ee} (GeV)", 30, 0, 300, "GeV", hMETee,dir);
  defineHistoset("hMETmm", "MET",
                  "#slash{E}_{T}^{#mu#mu} (GeV)", 30, 0, 300, "GeV", hMETmm,dir);
  defineHistoset("hMET3e0mu", "MET",
                  "#slash{E}_{T}^{3e0#mu} (GeV)", 30, 0, 300, "GeV", hMET3e0mu,dir);
  defineHistoset("hMET2e1mu", "MET",
                  "#slash{E}_{T}^{2e1#mu} (GeV)", 30, 0, 300, "GeV", hMET2e1mu,dir);
  defineHistoset("hMET1e2mu", "MET",
                  "#slash{E}_{T}^{1e2#mu} (GeV)", 30, 0, 300, "GeV", hMET1e2mu,dir);
  defineHistoset("hMET0e3mu", "MET",
                  "#slash{E}_{T}^{0e3#mu} (GeV)", 30, 0, 300, "GeV", hMET0e3mu,dir);

//W Trans Mass Histos
  defineHistoset("hWTransMass", "Reconstructed Transverse Mass of W",
                  "M_{T} (GeV)", 20, 0, 100, "GeV", hWTransMass,dir);
  defineHistoset("hWenuTransMass", "Reconstructed Transverse Mass of We\\nu",
                  "M_{T}^{e#nu} (GeV)", 20, 0, 100, "GeV", hWenuTransMass,dir);
  defineHistoset("hWmnuTransMass", "Reconstructed TransverseMass of W#mu\\nu",
                  "M_{T}^{#mu#nu} (GeV)", 20, 0, 100, "GeV", hWmnuTransMass,dir);
  defineHistoset("hW3e0muTransMass", "Reconstructed Transverse Mass of W(3e0#mu)",
                  "M_{T}^{3e0#mu} (GeV)", 20, 0, 100, "GeV", hW3e0muTransMass,dir);
  defineHistoset("hW2e1muTransMass", "Reconstructed Transverse Mass of W(2e1#mu)",
                  "M_{T}^{2e1#mu} (GeV)", 20, 0, 100, "GeV", hW2e1muTransMass,dir);
  defineHistoset("hW1e2muTransMass", "Reconstructed Transverse Mass of W(1e2#mu)",
                  "M_{T}^{1e2#mu} (GeV)", 20, 0, 100, "GeV", hW1e2muTransMass,dir);
  defineHistoset("hW0e3muTransMass", "Reconstructed Transverse Mass of W(0e3#mu)",
                  "M_{T}^{0e3#mu} (GeV)", 20, 0, 100, "GeV", hW0e3muTransMass,dir);

//Wpt Histos
  defineHistoset("hWpt", "p_{T}^{W}", 
                  "p_{T}^{W} (GeV)", 40, 0, 400, "GeV", hWpt,dir);
  defineHistoset("hWptVee", "p_{T}^{W,V#rightarrowee}", 
                  "p_{T}^{W,V#rightarrowee} (GeV)", 40, 0, 400, "GeV", hWptVee,dir);
  defineHistoset("hWptVmm", "p_{T}^{W,V#rightarrow#mu#mu}", 
                  "p_{T}^{W,V#rightarrow#mu#mu} (GeV)", 40, 0, 400, "GeV", hWptVmm,dir);

//W Charge Histos
  defineHistoset("hWQ", "Reconstructed Charge of W",
                  "q_{W}", 3, -1, 1, "", hWQ,dir);
  defineHistoset("hWenuQ", "Reconstructed Charge of We\\nu",
                  "q_{W}^{e#nu}", 3, -1.5, 1.5, "", hWenuQ,dir);
  defineHistoset("hWmnuQ", "Reconstructed TransverseMass of W#mu\\nu",
                  "q_{W}^{#mu#nu}", 3, -1.5, 1.5, "", hWmnuQ,dir);
  defineHistoset("hW3e0muQ", "Reconstructed Charge of W(3e0#mu)",
                  "q_{W}^{3e0#mu}", 3, -1.5, 1.5, "", hW3e0muQ,dir);
  defineHistoset("hW2e1muQ", "Reconstructed Charge of W(2e1#mu)",
                  "q_{W}^{2e1#mu}", 3, -1.5, 1.5, "", hW2e1muQ,dir);
  defineHistoset("hW1e2muQ", "Reconstructed Charge of W(1e2#mu)",
                  "q_{W}^{1e2#mu}", 3, -1.5, 1.5, "", hW1e2muQ,dir);
  defineHistoset("hW0e3muQ", "Reconstructed Charge of W(0e3#mu)",
                  "q_{W}^{0e3#mu}", 3, -1.5, 1.5, "", hW0e3muQ,dir);

  defineHistoset("hNLElec", "Number of Loose Electrons in Event",
                  "N_{e}^{Loose}", 10, 0, 10, "NONE", hNLElec,dir);
  defineHistoset("hNLMuon", "Number of Loose Muons in Event",
                  "N_{#mu}^{Loose}", 10, 0, 10, "NONE", hNLMuon,dir);
  defineHistoset("hNLLeps", "Number of Loose Leptons in Event",
                  "N_{l}^{Loose}", 10, 0, 10, "NONE", hNLLeps,dir);
  defineHistoset("hNLLepsVee", "Number of Loose Leptons in Event, V#rightarrowee",
                  "N_{l}^{Loose,V#rightarrowee}", 10, 0, 10, "NONE", hNLLepsVee,dir);
  defineHistoset("hNLLepsVmm", "Number of Loose Leptons in Event",
                  "N_{l}^{Loose,V#rightarrow#mu#mu}", 10, 0, 10, "NONE", hNLLepsVmm,dir);

  defineHistoset("hNTElec", "Number of Tight Electrons in Event",
                  "N_{e}", 10, 0, 10, "NONE", hNTElec,dir);
  defineHistoset("hNTMuon", "Number of Tight Muons in Event",
                  "N_{#mu}", 10, 0, 10, "NONE", hNTMuon,dir);
  defineHistoset("hNTLeps", "Number of Tight Leptons in Event",
                  "N_{l}", 10, 0, 10, "NONE", hNTLeps,dir);

  defineHistoset("hNJets", "Number of Jets in Event",
                  "N_{Jets}", 10, 0, 10, "NONE", hNJets,dir);
  defineHistoset("hNJetsVee", "Number of Jets in Event, V#rightarrowee",
                  "N_{Jets}^{V#rightarrowee}", 10, 0, 10, "NONE", hNJetsVee,dir);
  defineHistoset("hNJetsVmm", "Number of Jets in Event, V#rightarrow#mu#mu",
                  "N_{Jets}^{V#rightarrow#mu#mu}", 10, 0, 10, "NONE", hNJetsVmm,dir);

  defineHistoset("hNVtxs", "Number of Vertexs in Event",
                  "N_{Vtx}", 50, 0, 50, "NONE", hNVtxs,dir);
  defineHistoset("hNVtxsVee", "Number of Vertexs in Event, V#rightarrowee",
                  "N_{Vtx}^{V#rightarrowee}", 50, 0, 50, "NONE", hNVtxsVee,dir);
  defineHistoset("hNVtxsVmm", "Number of Vertexs in Event, V#rightarrow#mu#mu",
                  "N_{Vtx}^{V#rightarrow#mu#mu}", 50, 0, 50, "NONE", hNVtxsVmm,dir);

  defineHistoset("hWenuCombRelIso", "Comb Rel Iso of W Electron",
                  "Electron Combined Relative Isolation", 20, 0, 0.2, "NONE", hWenuCombRelIso,dir);
  defineHistoset("hWmnuCombRelIso", "Comb Rel Iso of W Muon",
                  "Muon Combined Relative Isolation", 20, 0, 0.2, "NONE", hWmnuCombRelIso,dir);
  

  tVWCand = dir.make<TTree>("tVWCand", "Analysis Variables after VWCand");//Only 1 for now;
  tVWCand->Branch("VWMass", &VWMass_);
  tVWCand->Branch("EvtType", &evtType_);
  tVWCand->Branch("Vpt", &Vpt_);
  tVWCand->Branch("Wpt", &Wpt_);
  tVWCand->Branch("weight", &weight_);

}//defineHistos

//fill Histograms
void HadronicVWAnalyzer::fillHistos(const int& index, const float& weight){
  if(debugme) printf("filling Histos\n");
  if(wCand_ && vCand_){
    hVWMass[index]->Fill(vwCand_.mass("minPz"), weight);
    if     (evtType_ == 0) hVW3e0muMass[index]->Fill(vwCand_.mass("minPz"), weight);
    else if(evtType_ == 1) hVW2e1muMass[index]->Fill(vwCand_.mass("minPz"), weight);
    else if(evtType_ == 2) hVW1e2muMass[index]->Fill(vwCand_.mass("minPz"), weight);
    else if(evtType_ == 3) hVW0e3muMass[index]->Fill(vwCand_.mass("minPz"), weight);
    hQ[index]->Fill(Q_, weight); 
    hVWTransMass[index]->Fill(vwCand_.transMass(), weight);
    hVWpt[index]->Fill(vwCand_.pt(), weight);
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
    if(CutNames_[index] == "ValidVWCand"){
      tVWCand->Fill();
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
      hWenuCombRelIso[index]->Fill(calcCombRelIso(e.patEle(), ElecPU(e)), weight);
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

}//fillHistos

inline void
HadronicVWAnalyzer::calcVVariables(){
  if (debugme) cout<<"In calc V Variables\n";
  vCand_ = getWCand(looseJets_);
  Vpt_ = vCand_.pt();
  if(debugme) printf("    Contains: %i tight V candidate(s)\n", (bool)vCand_);
  if(debugme){
    printEventLeptons(); 
    printEventDetails();
  }
}

inline void
HadronicVWAnalyzer::calcWVariables(){
  if (debugme) cout<<"In calc W Variables\n";
  wCand_ = getWCand(tightElectrons_, tightMuons_, met_);
  Wpt_ = wCand_.pt();
  if(debugme) printf("    Contains: %i tight W candidate(s)\n", (bool)wCand_);
  if(debugme){
    printEventLeptons(); 
    printEventDetails();
  }
}

inline void
HadronicVWAnalyzer::calcWElecVariables(){
  if (debugme) cout<<"In calc W Elec Variables\n";
  wCand_ = getWCand(tightElectrons_, met_);
  if(debugme) printf("    Contains: %i tight W candidate(s)\n", (bool)wCand_);
}

inline void
HadronicVWAnalyzer::calcWMuonVariables(){
  if (debugme) cout<<"In calc W Muon Variables\n";
  wCand_ = getWCand(tightMuons_, met_);
  if(debugme) printf("    Contains: %i tight W candidate(s)\n", (bool)wCand_);
}

inline void
HadronicVWAnalyzer::calcVWVariables(){
  if (debugme) cout<<"In calc VW Variables\n";
  vwCand_ = (vCand_ && wCand_) ? XWLeptonic(vCand_, wCand_) : XWLeptonic();
  VWMass_ = vwCand_.mass("minPz");
  Q_ = (vCand_ && wCand_) ? calcQ() : -999.;
  if(debugme) printEventDetails();
}

void
HadronicVWAnalyzer::calcEventVariables(){
  if (debugme) cout<<"In calc Event Variables\n";
  evtType_ = (vCand_ && wCand_) ? calcEvtType() : -999;
  if(debugme) printf("evt Type: %i, V Flav: %i, W Flav: %i\n", evtType_, (int)vCand_.flavor(), (int)wCand_.flavor());
}

void 
HadronicVWAnalyzer::eventLoop(edm::EventBase const & event){
  clearEvtVariables();
  if(debugme) WPrimeUtil::printEvent(event);

  // Preselection - skip events that don't look promising
  if (doPreselect_){
    if(debugme) cout<<"Testing Preselection...\n";
  }

  //get Jets
  event.getByLabel(jetsLabel_,patJetsH_);
  if(debugme) printf("    Contains: %i pat jets(s)\n",
                     (int)patJetsH_->size());
  for (size_t i = 0; i < patJetsH_->size(); i++) {
    const pat::Jet & jet = (*patJetsH_.product())[i];
    if (looseJet_(jet, jetLooseResult_))
      looseJets_.push_back(jet);
  }
  if(looseJets_.size() == 0) return;


  // get leptons
  ////////Cory:Not using electrons yet
  //event.getByLabel(electronsLabel_,patElectronsH_);
  event.getByLabel(muonsLabel_,patMuonsH_);
  //Cory: Didn't save pfMET
  //event.getByLabel(metLabel_, metH_);
  if(useAdjustedMET_) event.getByLabel(pfCandsLabel_, pfCandidatesH_);
  WPrimeUtil::getMuonsMET(patMuonsH_, muReconstructor_, allMuons_,
                          metH_, useAdjustedMET_, met_,
                          pfCandidatesH_);
  
/*
  WPrimeUtil::getLeptonsMET(patElectronsH_, allElectrons_,
                            patMuonsH_, muReconstructor_, allMuons_,
                            metH_, useAdjustedMET_, met_,
                            pfCandidatesH_);
  if(debugme) printf("    Contains: %i electron(s), %i muon(s)\n",
                          (int)allElectrons_.size(), (int)allMuons_.size());

  rhoFastJet_ = getProduct<double>(event,"kt6PFJets:rho");


  // Make vectors of leptons passing various criteria
  for (size_t i = 0; i < allElectrons_.size(); i++) {
    if(Overlap(allElectrons_[i].patEle(), allMuons_)) continue;
    const float pu = ElecPU(allElectrons_[i]);
    if (looseElectron_(allElectrons_[i].patEle(), electronLooseResult_, pu))
      looseAllElectrons_.push_back(allElectrons_[i]);

    if (tightElectron_(allElectrons_[i].patEle(), electronTightResult_, pu))
      tightElectrons_.push_back(allElectrons_[i]);
  }
*/
  for (size_t i = 0; i < allMuons_.size(); i++) {
    const float pu = MuonPU(allMuons_[i]);
    if (looseMuon_(allMuons_[i], muonLooseResult_,pu))
      looseMuons_.push_back(allMuons_[i]);

    if (tightMuon_(allMuons_[i], muonTightResult_,pu))
      tightMuons_.push_back(allMuons_[i]);
  }
  if(looseElectrons_.size() + looseMuons_.size() == 0) return;

  if(debugme){
    printLeptons();
    printf("    Contains: %i loose electron(s), %i loose muon(s), %i loose jet(s)\n",
           (int)looseElectrons_.size(), (int)looseMuons_.size(), (int)looseJets_.size());
    printf("    Contains: %i tight electron(s), %i tightmuon(s)\n",
           (int)tightElectrons_.size(), (int)tightMuons_.size());
  }

  ///////////////////
  //Cory: Remove overlap btw jets and leptons! 
  //(maybe I don't have to for pfJets)
  ///////////////////


  //get Trigger 
  triggerEvent_ = getProduct<pat::TriggerEvent>(event,hltEventLabel_); 

  //get Vertex
  vertices_ = getProduct<vector<reco::Vertex> >(event,vertexLabel_);

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
  }//MC Only If

  if(wprimeUtil_->DebugEvent(event)){
    cout<<"This is a debug event\n";
    printPassingEvent(event);
    printDebugEvent();
  }

  weight_ = wprimeUtil_->getWeight();
  if(!passCuts(weight_)) return;
  if(wprimeUtil_->runningOnData()){
    cout<<" The following data events passed All Cuts!!!\n";
    printPassingEvent(event);
    if(debugme) printEventLeptons();
    cout<<" ------------------\n";
  }
  if(debugme) printEventLeptons();
}

void HadronicVWAnalyzer::printDebugEvent() const{
  WPrimeUtil::printPassingTriggers(triggerEvent_,triggersToUse_);
  printEventDetails();
  printEventLeptons();
  printLeptons();
}

void HadronicVWAnalyzer::printEventDetails() const{
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
  if(vCand_ && wCand_ && vwCand_.mass("minPz")>0.){
    cout<<" VW Mass: "<<vwCand_.mass("minPz")
        <<" Neu Pz: "<<vwCand_.neutrinoPz("minPz")
        <<" Vpt: "<<vCand_.pt()
        <<" Wpt: "<<wCand_.pt()
        <<endl;
  }
  return;
}

void
HadronicVWAnalyzer::printEventLeptons() const{
  if     (vCand_.flavor() == PDGELEC){
//    printElectron(*vCand_.elec1(), PDGZ);
//    printElectron(*vCand_.elec2(), PDGZ);
    printElectron(WPrimeUtil::Find(*vCand_.daughter(0), allElectrons_));
    printElectron(WPrimeUtil::Find(*vCand_.daughter(1), allElectrons_));
  }else if(vCand_.flavor() == PDGMUON){
    printMuon(WPrimeUtil::Find(*vCand_.daughter(0), allMuons_));
    printMuon(WPrimeUtil::Find(*vCand_.daughter(1), allMuons_));
//    printMuon(*vCand_.muon1(), PDGZ);
//    printMuon(*vCand_.muon2(), PDGZ);
  }

  if     (wCand_.flavor() == PDGELEC){   
    printElectron(WPrimeUtil::Find(*wCand_.daughter(0), allElectrons_));
//    printElectron(*wCand_.elec(), PDGW);
  }else if(wCand_.flavor() == PDGMUON){
    printMuon(WPrimeUtil::Find(*wCand_.daughter(0), allMuons_));
//    printMuon    (*wCand_.muon(), PDGW);
  }
}

/////////////////Accessors///////////////////////

/////////////////Modifiers///////////////////////

/////////////////Cuts///////////////////////
bool
HadronicVWAnalyzer::passCuts(const float& weight){
  if (debugme) cout<<"In pass Cuts\n";
  
  for(int i=0; i<NCuts_; ++i){
    if(CutNames_[i] == "ValidV"){
      calcVVariables();
    }else if(CutNames_[i] == "ValidW"){
      calcWVariables();
      calcEventVariables();
    }else if(CutNames_[i] == "ValidVWCand"){
      calcVWVariables();
    }

    if(!CutFns_[i]()) return false;
    tabulateEvent(i,weight); 
  }
  return true;
}

//Trigger requirements
//-----------------------------------------------------------
bool HadronicVWAnalyzer::passTriggersCut() const{
  if(debugme) cout<<"Trigger requirements"<<endl;
  //Apply the trigger if running on data or MC 
  //If MC, apply if no V or if V exists, zCand == PDGMuon)
  if(wprimeUtil_->runningOnData() || !vCand_ || vCand_.flavor() == PDGMUON){
    return WPrimeUtil::passTriggersCut(triggerEvent_,triggersToUse_);
  }else{
    return true;//Cory: This is not good, but will pass HLT in the meantime.
  }
  return false;
}//--- passTriggersCut()


inline bool HadronicVWAnalyzer::passValidVWCut() const{
  return vwCand_ && vwCand_.mass("minPz")>0.;
}

////////////////////////////////
/////////Check Electron Properties/////
////////////////////////////////
bool HadronicVWAnalyzer::passTriggerEmulation(const heep::Ele& elec, const float minPt) const{
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

inline float HadronicVWAnalyzer::calcQ() const{
  return vwCand_.mass("minPz") - vCand_.mass() - WMASS;
}

inline int HadronicVWAnalyzer::calcEvtType() const{
  return (vCand_ && wCand_) ?  2 * (vCand_.flavor() != 11) + (wCand_.flavor() != 11) : -999;
}

inline bool HadronicVWAnalyzer::inEE(const TeVMuon& mu) const{
  return fabs(mu.eta()) >= 1.05;
}

///////////////Utilities//////////////////
//--------------------------------------------------------------

inline void
HadronicVWAnalyzer::clearEvtVariables(){
  AnalyzerBase::clearEvtVariables();
  met_ = pat::MET();
  vwCand_ = XWLeptonic();
  evtType_ = -999;
  VWMass_ = -999;
  Vpt_ = -999;
  Wpt_ = -999;
  Q_ = -999;
  weight_ = 0;
  rhoFastJet_ = 0;
}

float
HadronicVWAnalyzer::WLepPt() const{
  if(wCand_.flavor() == PDGELEC){
    return WPrimeUtil::Find(*wCand_.daughter(0), allElectrons_).patEle().pt();
  }else if(wCand_.flavor() == PDGMUON){
    return WPrimeUtil::Find(*wCand_.daughter(0), allMuons_).pt();
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
