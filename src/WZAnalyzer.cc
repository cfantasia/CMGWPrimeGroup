#include "UserCode/CMGWPrimeGroup/interface/WZAnalyzer.h"

using namespace std;

WZAnalyzer::WZAnalyzer(){}
WZAnalyzer::WZAnalyzer(const edm::ParameterSet & cfg, WPrimeUtil * wprimeUtil) :
  AnalyzerBase(cfg, wprimeUtil){
  setupCutOrder();
  if(debugme) printf("Using %i cuts\n",NCuts_);

  maxZMassDiff_ = cfg.getParameter<double>("maxZMassDiff");
  minDeltaR_ = cfg.getParameter<double>("minDeltaR");

  effectiveElecArea_ = cfg.getParameter<vector<double> >("effectiveElecArea");
  effectiveMuonArea_ = cfg.getParameter<vector<double> >("effectiveMuonArea");
  
// +++++++++++++++++++General Cut values
  maxNumZs_ = cfg.getParameter<uint>("maxNumZs");
  minLeadPt_ = cfg.getParameter<double>("minLeadPt");

// +++++++++++++++++++Ht Cuts
  minHt_ = cfg.getParameter<double>("minHt");

// +++++++++++++++++++W Cuts
  minWlepPt_ = cfg.getParameter<double>("minWlepPt");

// +++++++++++++++++++Z Cuts
  minZeePt1_ = cfg.getParameter<double>("minZeePt1");
  minZeePt2_ = cfg.getParameter<double>("minZeePt2");
  minZmmPt1_ = cfg.getParameter<double>("minZmmPt1");
  minZmmPt2_ = cfg.getParameter<double>("minZmmPt2");
}

WZAnalyzer::~WZAnalyzer(){
}

void WZAnalyzer::setupCutOrder(){
  mFnPtrs_["NoCuts"] = &WZAnalyzer::passNoCut;
  mFnPtrs_["HLT"] = &WZAnalyzer::passTriggersCut;
  mFnPtrs_["MinNLeptons"] = &WZAnalyzer::passMinNLeptonsCut;
  mFnPtrs_["MaxNLeptons"] = &WZAnalyzer::passMaxNLeptonsCut;
  mFnPtrs_["ValidW"] = &WZAnalyzer::passValidWCut;
  mFnPtrs_["ValidWElec"] = &WZAnalyzer::passValidWElecCut;
  mFnPtrs_["ValidWMuon"] = &WZAnalyzer::passValidWMuonCut;
  mFnPtrs_["ValidZ"] = &WZAnalyzer::passValidZCut;
  mFnPtrs_["ValidWZCand"] = &WZAnalyzer::passValidWZCut;
  mFnPtrs_["LeadLepPt"] = &WZAnalyzer::passLeadingLeptonPtCut;
  mFnPtrs_["NumZs"] = &WZAnalyzer::passNumberOfZsCut;
  mFnPtrs_["ZMass"] = &WZAnalyzer::passZMassCut;
  mFnPtrs_["WTransMass"] = &WZAnalyzer::passWtransMassCut;
  mFnPtrs_["MET"] = &WZAnalyzer::passMinMETCut;
  mFnPtrs_["Ht"] = &WZAnalyzer::passHtCut;
  mFnPtrs_["Zpt"] = &WZAnalyzer::passZptCut;
  mFnPtrs_["Wpt"] = &WZAnalyzer::passWptCut;
  mFnPtrs_["ZLepPt"] = &WZAnalyzer::passZLepPtCut;
  mFnPtrs_["WLepPt"] = &WZAnalyzer::passWLepPtCut;
  mFnPtrs_["ZLepTrigMatch"] = &WZAnalyzer::passZLepTriggerMatchCut;
  mFnPtrs_["AllCuts"] = &WZAnalyzer::passNoCut;

  mFnPtrs_["WLepTight"] = &WZAnalyzer::passWLepTightCut;
  mFnPtrs_["WFlavorElec"] = &WZAnalyzer::passWFlavorElecCut;
  mFnPtrs_["WFlavorMuon"] = &WZAnalyzer::passWFlavorMuonCut;
  mFnPtrs_["FakeEvt"] = &WZAnalyzer::passFakeEvtCut;
  mFnPtrs_["FakeLepTag"] = &WZAnalyzer::passFakeLeptonTagCut;
  mFnPtrs_["FakeLepProbeTight"] = &WZAnalyzer::passFakeLeptonProbeTightCut;

  CutFns_.resize(NCuts_);
  for(int i=0; i<NCuts_; ++i){
    if(mFnPtrs_.find(Cuts_[i]) == mFnPtrs_.end()){
      cout<<"Didn't find Cut named "<<Cuts_[i]<<endl;
      abort();
    }
    CutFns_[i] = mFnPtrs_.find(Cuts_[i])->second;
  } 

}

void WZAnalyzer::defineHistos(const TFileDirectory & dir){
  if(debugme) printf("Declare histos\n");
  AnalyzerBase::defineHistos(dir);

  defineHistoset("hWZMass", "Reconstructed WZ Invariant Mass",
                  "M_{WZ} (GeV)", 250, 0, 2500, "GeV", hWZMass,dir);
  defineHistoset("hWZ3e0muMass", "Reconstructed WZ(3e0#mu) Invariant Mass",
                  "M_{WZ}^{3e0#mu} (GeV)", 250, 0, 2500, "GeV", hWZ3e0muMass,dir);
  defineHistoset("hWZ2e1muMass", "Reconstructed WZ(2e1#mu) Invariant Mass",
                  "M_{WZ}^{2e1#mu} (GeV)", 250, 0, 2500, "GeV", hWZ2e1muMass,dir);
  defineHistoset("hWZ1e2muMass", "Reconstructed WZ(1e2#mu) Invariant Mass",
                  "M_{WZ}^{1e2#mu} (GeV)", 250, 0, 2500, "GeV", hWZ1e2muMass,dir);
  defineHistoset("hWZ0e3muMass", "Reconstructed WZ(0e3#mu) Invariant Mass",
                  "M_{WZ}^{0e3#mu} (GeV)", 250, 0, 2500, "GeV", hWZ0e3muMass,dir);

//Q=M_{WZ} - M_W - M_Z
  defineHistoset("hQ", "Q=M_{WZ} - M_{W} - M_{Z}",
                  "Q (GeV)", 50, 0, 500, "GeV", hQ,dir);
  defineHistoset("hWZTransMass", "Reconstructed WZ Transverse Mass",
                  "M_{WZ}^{T} (GeV)", 250, 0, 2500, "GeV", hWZTransMass,dir);
//WZpt Histos
  defineHistoset("hWZpt", "Reconstructed WZ Transverse Momentum",
                  "p_{WZ}^{T} (GeV)", 50, 0, 500, "GeV", hWZpt,dir);
  defineHistoset("hWZTheta", "Theta of WZ",
                  "#theta_{WZ}", 50, 0, PI, "NONE", hWZTheta,dir);

//Ht Histos
  defineHistoset("hHt", "H_{T}", 
                  "Lepton Pt Sum: H_{T} (GeV)", 80, 0, 800, "GeV", hHt,dir);
  defineHistoset("hTriLepMass", "hTriLepMass",
                  "Trilepton Invariant Mass", 100, 0., 1000., "GeV", hTriLepMass, dir);
  defineHistoset("hEvtType", "Event Type",
                  "N_{#mu}", 4, 0, 4, "NONE", hEvtType,dir);
  defineHistoset("hEvtTypeP", "Event Type for Q=+1",
                  "N_{#mu},W^{+}", 4, 0, 4, "NONE", hEvtTypeP,dir);
  defineHistoset("hEvtTypeM", "Event Type for Q=-1",
                  "N_{#mu},W^{-}", 4, 0, 4, "NONE", hEvtTypeM,dir);
//Lead Lepton Pt
  defineHistoset("hLeadPt", "Leading Lepton Pt",
                  "p_{T}^{Max}", 50, 0, 1000., "GeV", hLeadPt,dir);
  defineHistoset("hLeadPtZee", "Leading Lepton Pt Zee",
                  "p_{T}^{Max, ee}", 50, 0, 1000., "GeV", hLeadPtZee,dir);
  defineHistoset("hLeadPtZmm", "Leading Lepton Pt Zmm",
                  "p_{T}^{Max #mu#mu}", 50, 0, 1000., "GeV", hLeadPtZmm,dir);
  defineHistoset("hLeadElecPt", "Leading Electron Pt",
                  "p_{T}^{Max e}", 50, 0, 1000., "GeV", hLeadElecPt,dir);
  defineHistoset("hLeadMuonPt", "Leading Muon Pt",
                  "p_{T}^{Max #mu}", 50, 0, 1000., "GeV", hLeadMuonPt,dir);

///////////////////////////
//Z Mass Histos
  defineHistoset("hZMass" , "Reconstructed Mass of Z",
                  "M_{Z} (GeV)", 30, 60, 120, "GeV", hZMass,dir);
  defineHistoset("hZeeMass","Reconstructed Mass of Zee",
                  "M_{Z}^{ee} (GeV)", 30, 60, 120, "GeV", hZeeMass,dir);
  defineHistoset("hZmmMass","Reconstructed Mass of Z#mu#mu",
                  "M_{Z}^{#mu#mu} (GeV)", 30, 60, 120, "GeV", hZmmMass,dir);
  defineHistoset("hZ3e0muMass" , "Reconstructed Mass of Z(3e0#mu)",
                  "M_{Z}^{3e0#mu} (GeV)", 30, 60, 120, "GeV", hZ3e0muMass,dir);
  defineHistoset("hZ2e1muMass" , "Reconstructed Mass of Z(2e1#mu)",
                  "M_{Z}^{2e1#mu} (GeV)", 30, 60, 120, "GeV", hZ2e1muMass,dir);
  defineHistoset("hZ1e2muMass" , "Reconstructed Mass of Z(1e2#mu)",
                  "M_{Z}^{1e2#mu} (GeV)", 30, 60, 120, "GeV", hZ1e2muMass,dir);
  defineHistoset("hZ0e3muMass" , "Reconstructed Mass of Z(0e3#mu)",
                  "M_{Z}^{0e3#mu} (GeV)", 30, 60, 120, "GeV", hZ0e3muMass,dir);
  defineHistoset("hZeeMassTT","Reconstructed MassTT of ZeeTT",
                  "M_{Z}^{ee,TT} (GeV)", 30, 60, 120, "GeV", hZeeMassTT,dir);
  defineHistoset("hZeeMassTF","Reconstructed Mass of ZeeTF",
                  "M_{Z}^{ee,TF} (GeV)", 30, 60, 120, "GeV", hZeeMassTF,dir);
  defineHistoset("hZmmMassTT","Reconstructed Mass of Z#mu#muTT",
                  "M_{Z}^{#mu#mu,TT} (GeV)", 30, 60, 120, "GeV", hZmmMassTT,dir);
  defineHistoset("hZmmMassTF","Reconstructed Mass of Z#mu#muTF",
                  "M_{Z}^{#mu#mu,TF} (GeV)", 30, 60, 120, "GeV", hZmmMassTF,dir);

//Zpt Histos
  defineHistoset("hZpt", "p_{T}^{Z}", 
                  "p_{T}^{Z} (GeV)", 50, 0, 1000, "GeV", hZpt,dir);
  defineHistoset("hZeept", "p_{T}^{Z#rightarrowee}", 
                  "p_{T}^{Z#rightarrowee} (GeV)", 50, 0, 1000, "GeV", hZeept,dir);
  defineHistoset("hZmmpt", "p_{T}^{Z#rightarrow#mu#mu}", 
                  "p_{T}^{Z#rightarrow#mu#mu} (GeV)", 50, 0, 1000, "GeV", hZmmpt,dir);
//MET Histos
  defineHistoset("hMET", "MET",
                  "#slash{E}_{T} (GeV)", 50, 0, 500, "GeV", hMET,dir);
  defineHistoset("hMETee", "MET",
                  "#slash{E}_{T}^{ee} (GeV)", 50, 0, 500, "GeV", hMETee,dir);
  defineHistoset("hMETmm", "MET",
                  "#slash{E}_{T}^{#mu#mu} (GeV)", 50, 0, 500, "GeV", hMETmm,dir);
  defineHistoset("hMET3e0mu", "MET",
                  "#slash{E}_{T}^{3e0#mu} (GeV)", 50, 0, 500, "GeV", hMET3e0mu,dir);
  defineHistoset("hMET2e1mu", "MET",
                  "#slash{E}_{T}^{2e1#mu} (GeV)", 50, 0, 500, "GeV", hMET2e1mu,dir);
  defineHistoset("hMET1e2mu", "MET",
                  "#slash{E}_{T}^{1e2#mu} (GeV)", 50, 0, 500, "GeV", hMET1e2mu,dir);
  defineHistoset("hMET0e3mu", "MET",
                  "#slash{E}_{T}^{0e3#mu} (GeV)", 50, 0, 500, "GeV", hMET0e3mu,dir);
  defineHistoset("hMETSig", "MET Significance",
                  "#slash{E}_{T}^{Signif}", 50, 0, 500, "NONE", hMETSig,dir);

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
                  "p_{T}^{W} (GeV)", 50, 0, 1000, "GeV", hWpt,dir);
  defineHistoset("hWptZee", "p_{T}^{W,Z#rightarrowee}", 
                  "p_{T}^{W,Z#rightarrowee} (GeV)", 50, 0, 1000, "GeV", hWptZee,dir);
  defineHistoset("hWptZmm", "p_{T}^{W,Z#rightarrow#mu#mu}", 
                  "p_{T}^{W,Z#rightarrow#mu#mu} (GeV)", 50, 0, 1000, "GeV", hWptZmm,dir);

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
  defineHistoset("hWTheta", "Theta of W",
                  "#theta_{W}", 50, 0, PI, "NONE", hWTheta,dir);

  defineHistoset("hNLElec", "Number of Loose Electrons in Event",
                  "N_{e}^{Loose}", 10, 0, 10, "NONE", hNLElec,dir);
  defineHistoset("hNLMuon", "Number of Loose Muons in Event",
                  "N_{#mu}^{Loose}", 10, 0, 10, "NONE", hNLMuon,dir);
  defineHistoset("hNLLeps", "Number of Loose Leptons in Event",
                  "N_{Lep}^{Loose}", 10, 0, 10, "NONE", hNLLeps,dir);
  defineHistoset("hNLLepsZee", "Number of Loose Leptons in Event, Z#rightarrowee",
                  "N_{Lep}^{Loose,Z#rightarrowee}", 10, 0, 10, "NONE", hNLLepsZee,dir);
  defineHistoset("hNLLepsZmm", "Number of Loose Leptons in Event",
                  "N_{Lep}^{Loose,Z#rightarrow#mu#mu}", 10, 0, 10, "NONE", hNLLepsZmm,dir);

  defineHistoset("hNTElec", "Number of Tight Electrons in Event",
                  "N_{e}^{Tight}", 10, 0, 10, "NONE", hNTElec,dir);
  defineHistoset("hNTMuon", "Number of Tight Muons in Event",
                  "N_{#mu}^{Tight}", 10, 0, 10, "NONE", hNTMuon,dir);
  defineHistoset("hNTLeps", "Number of Tight Leptons in Event",
                  "N_{Lep}^{Tight}", 10, 0, 10, "NONE", hNTLeps,dir);

  defineHistoset("hNJets", "Number of Jets in Event",
                  "N_{Jets}", 15, 0, 15, "NONE", hNJets,dir);
  defineHistoset("hNJetsZee", "Number of Jets in Event, Z#rightarrowee",
                  "N_{Jets}^{Z#rightarrowee}", 15, 0, 15, "NONE", hNJetsZee,dir);
  defineHistoset("hNJetsZmm", "Number of Jets in Event, Z#rightarrow#mu#mu",
                  "N_{Jets}^{Z#rightarrow#mu#mu}", 15, 0, 15, "NONE", hNJetsZmm,dir);

  defineHistoset("hNVtxs", "Number of Vertexs in Event",
                  "N_{Vtx}", 50, 0, 50, "NONE", hNVtxs,dir);
  defineHistoset("hNVtxsZee", "Number of Vertexs in Event, Z#rightarrowee",
                  "N_{Vtx}^{Z#rightarrowee}", 50, 0, 50, "NONE", hNVtxsZee,dir);
  defineHistoset("hNVtxsZmm", "Number of Vertexs in Event, Z#rightarrow#mu#mu",
                  "N_{Vtx}^{Z#rightarrow#mu#mu}", 50, 0, 50, "NONE", hNVtxsZmm,dir);

  defineHistoset("hWenuCombRelIso", "Comb Rel Iso of W Electron",
                  "Electron Combined Relative Isolation", 20, 0, 0.2, "NONE", hWenuCombRelIso,dir);
  defineHistoset("hWmnuCombRelIso", "Comb Rel Iso of W Muon",
                  "Muon Combined Relative Isolation", 20, 0, 0.2, "NONE", hWmnuCombRelIso,dir);
  

  tWZCand = dir.make<TTree>("tWZCand", "Analysis Variables after WZCand");//Only 1 for now;
  tWZCand->Branch("Run", &runNumber_);
  tWZCand->Branch("Lumi", &lumiNumber_);
  tWZCand->Branch("Event", &evtNumber_);
  tWZCand->Branch("WZMass", &WZMass_);
  tWZCand->Branch("EvtType", &evtType_);
  tWZCand->Branch("Ht", &Ht_);
  tWZCand->Branch("Zpt", &Zpt_);
  tWZCand->Branch("Wpt", &Wpt_);
  tWZCand->Branch("MET", &MET_);
  tWZCand->Branch("Q", &Q_);
  tWZCand->Branch("TriLepMass", &TriLepMass_);
  tWZCand->Branch("weight", &weight_);

}//defineHistos

//fill Histograms
void WZAnalyzer::fillHistos(const int& index, const float& weight){
  if(debugme) printf("filling Histos\n");
  if(wCand_ && zCand_){
    hWZMass[index]->Fill(wzCand_.mass("minPz"), weight);
    if     (evtType_ == 0) hWZ3e0muMass[index]->Fill(wzCand_.mass("minPz"), weight);
    else if(evtType_ == 1) hWZ2e1muMass[index]->Fill(wzCand_.mass("minPz"), weight);
    else if(evtType_ == 2) hWZ1e2muMass[index]->Fill(wzCand_.mass("minPz"), weight);
    else if(evtType_ == 3) hWZ0e3muMass[index]->Fill(wzCand_.mass("minPz"), weight);
    hQ[index]->Fill(Q_, weight); 
    hWZTransMass[index]->Fill(wzCand_.transMass(), weight);
    hWZpt[index]->Fill(wzCand_.pt(), weight);
    /////////////hWZTheta[index]->Fill(wzCand_.theta(), weight);
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
      hNJetsZee[index]->Fill(allJets_.size(), weight);
      hNVtxsZee[index]->Fill(vertices_.size(), weight);
    }else if(zCand_.flavor() == PDGMUON){ 
      hLeadPtZmm[index]->Fill(LeadPt_, weight);
      hWptZmm[index]->Fill(wCand_.pt(), weight);
      hNLLepsZmm[index]->Fill(looseElectrons_.size()+looseMuons_.size(), weight);
      hNJetsZmm[index]->Fill(allJets_.size(), weight);
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
    ///////////hWTheta[index]->Fill(wCand_.theta(), weight);
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
  hMETSig[index]->Fill(met_.significance(), weight);

  hNLElec[index]->Fill(looseElectrons_.size(), weight);
  hNLMuon[index]->Fill(looseMuons_    .size(), weight);
  hNLLeps[index]->Fill(looseElectrons_.size()+looseMuons_.size(), weight);

  hNTElec[index]->Fill(tightElectrons_.size(), weight);
  hNTMuon[index]->Fill(tightMuons_    .size(), weight);
  hNTLeps[index]->Fill(tightElectrons_.size()+tightMuons_.size(), weight);

  hNJets[index]->Fill(allJets_.size(), weight);
  hNVtxs[index]->Fill(vertices_.size(), weight);

}//fillHistos

int
WZAnalyzer::countZCands(ZCandV & Zs) const{
  int count =0;
  for(uint i=0; i<Zs.size(); ++i) 
    if( (Zs[i].mass() > minZmass_) && (Zs[i].mass() < maxZmass_)) 
      count++;
  return count;
}

void
WZAnalyzer::calcZVariables(){
  if (debugme) cout<<"In calc Z Variables\n";
  // Reconstruct the Z
  float matchptcut = 0.;
  bool minHighPt = false;
  
  matchptcut = 10.;
  ElectronV zElectrons;
  for (size_t i=0; i < looseElectrons_.size(); i++)
    if(passTriggerEmulation(looseElectrons_[i], matchptcut))
      zElectrons.push_back(looseElectrons_[i]);
  matchptcut = 20.;
  minHighPt = false;
  for (size_t i=0; i < zElectrons.size(); i++){
    if(passTriggerEmulation(looseElectrons_[i], matchptcut)){
      minHighPt= true; 
      break;  
    }
  }
  if(!minHighPt) zElectrons.clear();
  
  matchptcut = 8.;
  MuonV zMuons;
  for (size_t i=0; i < looseMuons_.size(); i++)
    if(passTriggerMatch(looseMuons_[i], matchptcut, triggersToUse_))
      zMuons.push_back(looseMuons_[i]);
  
  matchptcut = 13.;
  minHighPt = false;
  for (size_t i=0; i < zMuons.size(); i++){
    if(passTriggerMatch(zMuons[i], matchptcut, triggersToUse_)){
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
        if(!tight1 && WPrimeUtil::Match(tightElectrons_[i], *zCand_.daughter(0))) tight1 = true;
        if(!tight2 && WPrimeUtil::Match(tightElectrons_[i], *zCand_.daughter(1))) tight2 = true;
      }
    }else if(zCand_.flavor() == PDGMUON){
      for(uint i=0; i<tightMuons_.size(); ++i){
        if(!tight1 && WPrimeUtil::Match(tightMuons_[i], *zCand_.daughter(0))) tight1 = true;
        if(!tight2 && WPrimeUtil::Match(tightMuons_[i], *zCand_.daughter(1))) tight2 = true;
      } 
    }
    TT = tight1 && tight2;
    TF = (tight1 && !tight2) || (!tight1 && tight2);
  }
  zCandsAll.insert(zCandsAll.begin(), zCand_);
  removeOverlapping(zCandsAll);
  numZs_ = countZCands(zCandsAll); 
  Zpt_ = zCand_.pt();

  if(debugme){
    printf("    Contains: %i Z candidate(s)\n", (int)zCandsAll.size());
    printEventLeptons();
    printEventDetails();
  }
}

inline void
WZAnalyzer::calcWVariables(){
  if (debugme) cout<<"In calc W Variables\n";
  wCand_ = getWCand(tightElectrons_, tightMuons_, met_, zCand_, minDeltaR_);
  Wpt_ = wCand_.pt();
  if(debugme) printf("    Contains: %i tight W candidate(s)\n", (bool)wCand_);
  if(debugme){
    printEventLeptons(); 
    printEventDetails();
  }
}

inline void
WZAnalyzer::calcWElecVariables(){
  if (debugme) cout<<"In calc W Elec Variables\n";
  wCand_ = getWCand(tightElectrons_, met_);
  if(debugme) printf("    Contains: %i tight W candidate(s)\n", (bool)wCand_);
}

inline void
WZAnalyzer::calcWMuonVariables(){
  if (debugme) cout<<"In calc W Muon Variables\n";
  wCand_ = getWCand(tightMuons_, met_);
  if(debugme) printf("    Contains: %i tight W candidate(s)\n", (bool)wCand_);
}

inline void
WZAnalyzer::calcWZVariables(){
  if (debugme) cout<<"In calc WZ Variables\n";
  wzCand_ = (zCand_ && wCand_) ? XWLeptonic(zCand_, wCand_) : XWLeptonic();
  WZMass_ = wzCand_.mass("minPz");
  Q_ = (zCand_ && wCand_) ? calcQ() : -999.;
  if(debugme) printEventDetails();
}

void
WZAnalyzer::calcEventVariables(){
  if (debugme) cout<<"In calc Event Variables\n";
  evtType_ = (zCand_ && wCand_) ? calcEvtType() : -999;
  if(debugme) printf("evt Type: %i, Z Flav: %i, W Flav: %i\n", evtType_, (int)zCand_.flavor(), (int)wCand_.flavor());
  LeadPt_ = calcLeadPt(); 
  LeadElecPt_ = calcLeadPt(PDGELEC);
  LeadMuonPt_ = calcLeadPt(PDGMUON);
  Ht_ = (zCand_ && wCand_) ? calcHt() : -999.;
  TriLepMass_ = (zCand_ && wCand_) ? calcTriLepMass() : -999.;
  MET_ = met_.et();
}

void 
WZAnalyzer::eventLoop(edm::EventBase const & event){
  clearEvtVariables();
  runNumber_ = event.id().run();
  lumiNumber_ = event.id().luminosityBlock();
  evtNumber_ = event.id().event();
/*
  if(!wprimeUtil_->runningOnData()){//Don't do this for data
    GenParticleV genParticles = getProduct<GenParticleV>(event, "genParticles");
    for (size_t i = 0; i < genParticles.size(); i++){
      if (abs(genParticles[i].pdgId()) == PDGTAU){
        return;
      }
    }
  }
*/  
  if(debugme) WPrimeUtil::printEvent(event);

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
//  // get information about flavorHistory
//  const uint flavorType = getProduct<uint>(event, "flavorHistoryFilter", 0);
//  if(debugme) printf("    FlavorType: %i", flavorType);

  // get leptons
  event.getByLabel(electronsLabel_,patElectronsH_);
  //const vector<pat::Electron> patElectrons = getProduct<vector<pat::Electron> >(event, electronsLabel_);
  event.getByLabel(muonsLabel_,patMuonsH_);
  //const vector<pat::Muon    > patMuons     = getProduct<vector<pat::Muon    > >(event, muonsLabel_);
  event.getByLabel(metLabel_, metH_);

  if(useAdjustedMET_) event.getByLabel(pfCandsLabel_, pfCandidatesH_);
  WPrimeUtil::getLeptonsMET(patElectronsH_, allElectrons_,
                            patMuonsH_, muReconstructor_, allMuons_,
                            metH_, useAdjustedMET_, met_,
                            pfCandidatesH_);
  if(debugme) printf("    Contains: %i electron(s), %i muon(s)\n",
                          (int)allElectrons_.size(), (int)allMuons_.size());

  rhoFastJet_ = getProduct<double>(event,"kt6PFJets:rho");

  // Make vectors of leptons passing various criteria
  for (size_t i = 0; i < allElectrons_.size(); i++) {
    if(Overlap(allElectrons_[i].patEle(), *patMuonsH_.product(), 0.01)) continue;
    const float pu = ElecPU(allElectrons_[i]);
    if (looseElectron_(allElectrons_[i].patEle(), electronResult_, pu))
      looseElectrons_.push_back(allElectrons_[i]);

    if (tightElectron_(allElectrons_[i].patEle(), electronResult_, pu))
      tightElectrons_.push_back(allElectrons_[i]);
  }

  for (size_t i = 0; i < allMuons_.size(); i++) {
    const float pu = MuonPU(allMuons_[i]);
    if (looseMuon_(allMuons_[i], muonResult_,pu))
      looseMuons_.push_back(allMuons_[i]);

    if (tightMuon_(allMuons_[i], muonResult_,pu))
      tightMuons_.push_back(allMuons_[i]);
  }
  if(looseElectrons_.size() + looseMuons_.size() == 0) return;

  if(debugme){
    printLeptons();
    printf("    Contains: %i loose electron(s), %i loose muon(s)\n",
           (int)looseElectrons_.size(), (int)looseMuons_.size());
    printf("    Contains: %i tight electron(s), %i tightmuon(s)\n",
           (int)tightElectrons_.size(), (int)tightMuons_.size());
  }

  //get Jets
  allJets_  = getProduct< JetV>(event,jetsLabel_);

  //get Trigger 
  triggerEvent_ = getProduct<pat::TriggerEvent>(event,hltEventLabel_); 

  //get Vertex
  vertices_ = getProduct<vector<reco::Vertex> >(event,vertexLabel_);

  if(!wprimeUtil_->runningOnData()){//Don't do this for data
    if(debugme){
      GenParticleV genParticles = getProduct<GenParticleV>(event, "genParticles");
      const reco::Candidate * genZ = 0;
      const reco::Candidate * genW = 0;
      for (size_t i = 0; i < genParticles.size(); i++){
        if (abs(genParticles[i].pdgId()) == PDGZ){
          genZ = & genParticles[i];
          cout<<"Mass of gen Z is "<<genZ->mass()<<endl;
        }else if (abs(genParticles[i].pdgId()) == PDGW){ 
          genW = & genParticles[i];
        cout<<"Mass of gen W is "<<genW->mass()<<endl;
        }
      }
    }
  }//MC Only If

  if(debugme && wprimeUtil_->DebugEvent(event)){
    cout<<"This is a debug event\n";
    printPassingEvent(event);
    printDebugEvent();
  }

  weight_ = wprimeUtil_->getWeight();
  if(!passCuts(weight_)) return;
  if(wprimeUtil_->runningOnData()){
    cout<<" The following data event passed All Cuts!!!\n";
    printPassingEvent(event);
    if(1 || debugme) printEventLeptons();
    cout<<" ------------------\n";
  }
  if(debugme) printEventLeptons();

  if(0 && !wprimeUtil_->runningOnData()){//Don't do this for data
    GenParticleV genParticles = getProduct<GenParticleV>(event, "genParticles");
    const reco::Candidate * genE = 0;
    for (size_t i = 0; i < genParticles.size(); i++){
      if (abs(genParticles[i].pdgId()) == PDGELEC){
        genE = & genParticles[i];
        cout<<"Gen e : "
            <<" pt: "<<genE->pt()
            <<" eta: "<<genE->eta()
            <<" phi: "<<genE->phi()
            <<endl;
      }
    }
  }
/*
  const reco::Candidate * genE1 = zCand_.genLepton1();
  const reco::Candidate * genE2 = zCand_.genLepton2();
  cout<<"Gen e 2: "
      <<" pt: "<<genE2->pt()
      <<" eta: "<<genE2->eta()
      <<" phi: "<<genE2->phi()
      <<endl;
*/
}

void WZAnalyzer::printDebugEvent() const{
  WPrimeUtil::printPassingTriggers(triggerEvent_,triggersToUse_);
  printEventDetails();
  printEventLeptons();
  printLeptons();
}

void WZAnalyzer::printEventDetails() const{
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
WZAnalyzer::printEventLeptons() const{
  if     (zCand_.flavor() == PDGELEC){
    cout<<"------- Electron 1 from Z -------\n";
//    printElectron(*zCand_.elec1(), PDGZ);
    printElectron(WPrimeUtil::Find(*zCand_.daughter(0), allElectrons_));
    cout<<"------- Electron 2 from Z -------\n";
//    printElectron(*zCand_.elec2(), PDGZ);
    printElectron(WPrimeUtil::Find(*zCand_.daughter(1), allElectrons_));
  }else if(zCand_.flavor() == PDGMUON){
    cout<<"------- Muon 1 from Z -------\n";
//    printMuon(*zCand_.muon1(), PDGZ);
    printMuon(WPrimeUtil::Find(*zCand_.daughter(0), allMuons_));
    cout<<"------- Muon 2 from Z -------\n";
//    printMuon(*zCand_.muon2(), PDGZ);
    printMuon(WPrimeUtil::Find(*zCand_.daughter(1), allMuons_));
  }

  if     (wCand_.flavor() == PDGELEC){   
    cout<<"------- Electron from W -------\n";
    printElectron(WPrimeUtil::Find(*wCand_.daughter(0), allElectrons_));
//    printElectron(*wCand_.elec(), PDGW);
  }else if(wCand_.flavor() == PDGMUON){
    cout<<"------- Muon from W -------\n";
    printMuon(WPrimeUtil::Find(*wCand_.daughter(0), allMuons_));
//    printMuon    (*wCand_.muon(), PDGW);
  }
}

float
WZAnalyzer::calcLeadPt(int type) const{
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
  return TMath::Max(calcLeadPt(PDGELEC), calcLeadPt(PDGMUON));
}

/////////////////Accessors///////////////////////

/////////////////Modifiers///////////////////////

/////////////////Cuts///////////////////////
bool
WZAnalyzer::passCuts(const float& weight){
  if (debugme) cout<<"In pass Cuts\n";

  
  for(int i=0; i<NCuts_; ++i){
    if(Cuts_[i] == "ValidZ") calcZVariables();  
    else if(Cuts_[i] == "ValidW"){  
      calcWVariables();  
      calcEventVariables();  
    }else if(Cuts_[i] == "ValidWZCand") calcWZVariables();  
    else if(Cuts_[i] == "ValidWElec") calcWElecVariables();  
    else if(Cuts_[i] == "ValidWMuon") calcWMuonVariables();

    if(!(this->*CutFns_[i])()) return false;
    tabulateEvent(i,weight); 
  }
  return true;
}

inline bool WZAnalyzer::passWLepTightCut() const{
  if(wCand_.flavor() == PDGELEC){
    const heep::Ele & e = *wCand_.elec();
    return WPrimeUtil::Contains(e.patEle(), tightElectrons_);
  }else if(wCand_.flavor() == PDGMUON){
    const TeVMuon & m = *wCand_.muon();
    return WPrimeUtil::Contains(m, tightMuons_);
}
  return false;
}

//Trigger requirements
//-----------------------------------------------------------
bool WZAnalyzer::passTriggersCut() const{
  if(debugme) cout<<"Trigger requirements"<<endl;
  //Apply the trigger if running on data or MC 
  //If MC, apply if no Z or if Z exists, zCand == PDGMuon)
  if(wprimeUtil_->runningOnData() || !zCand_ || zCand_.flavor() == PDGMUON){
    return WPrimeUtil::passTriggersCut(triggerEvent_,triggersToUse_);
  }else{
    return true;//Cory: This is not good, but will pass HLT in the meantime.
  }
  return false;
}//--- passTriggersCut()

inline bool
WZAnalyzer::passValidWElecCut() const{
  return wCand_ && wCand_.mt()>0.;
}

inline bool
WZAnalyzer::passValidWMuonCut() const{
  return wCand_ && wCand_.mt()>0.;
}

inline bool
WZAnalyzer::passNumberOfZsCut() const{
  return numZs_ <= maxNumZs_;
}

inline bool
WZAnalyzer::passLeadingLeptonPtCut() const{
  return LeadPt_ > minLeadPt_;
}

////////////////////////////////
/////////Check Z Properties/////
////////////////////////////////
bool
WZAnalyzer::passZLepPtCut() const{
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
WZAnalyzer::passZLepTriggerMatchCut() const{
  if     (zCand_.flavor() == PDGELEC){ 
    const heep::Ele& e1 = *zCand_.elec1();
    const heep::Ele& e2 = *zCand_.elec2();
    return (passTriggerEmulation(e1) && passTriggerEmulation(e2));
  }else if(zCand_.flavor() == PDGMUON){
    return true;//Trigger matching now before making zs
    //TeVMuon& m1 = WPrimeUtil::Find(*zCand_.daughter(0));
    //TeVMuon& m2 = WPrimeUtil::Find(*zCand_.daughter(1));
    //return passTriggerMatch(m1, m2); 
  }
  return false;
}

bool 
WZAnalyzer::passTriggerMatch(const pat::Electron & p, const float cut, const vstring& triggers) const{
  for (size_t i=0; i < triggers.size(); ++i){
    if (p.triggerObjectMatchesByPath(triggers[i], true, false).size() > 0){
      const pat::TriggerObjectStandAlone * trigRef = p.triggerObjectMatchByPath(triggers[i], true, false);
      if(trigRef->et() > cut) return true;
    }
  }
  return false;
}

bool 
WZAnalyzer::passTriggerMatch(const TeVMuon & p, const float cut, const vstring& triggers) const{
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
WZAnalyzer::passTriggerMatch(const heep::Ele& e1, const heep::Ele& e2) const{
  return (e1.patEle().triggerObjectMatches().size() > 0 &&
          e2.patEle().triggerObjectMatches().size() > 0 &&
          max(e1.patEle().triggerObjectMatches()[0].et(), e2.patEle().triggerObjectMatches()[0].et()) > 17. &&
          min(e1.patEle().triggerObjectMatches()[0].et(), e2.patEle().triggerObjectMatches()[0].et()) > 8.);
}

inline bool
WZAnalyzer::passTriggerMatch(const TeVMuon& m1, const TeVMuon& m2) const{
  return (m1.triggerObjectMatches().size() > 0 &&
          m2.triggerObjectMatches().size() > 0 &&
          max(m1.triggerObjectMatches()[0].pt(), m2.triggerObjectMatches()[0].pt()) > 13. &&
          min(m1.triggerObjectMatches()[0].pt(), m2.triggerObjectMatches()[0].pt()) > 8.);
}


////////////////////////////////
/////////Check W Properties/////
////////////////////////////////

//Check W Transverse Mass
inline bool
WZAnalyzer::passWLepPtCut() const{
  return WLepPt() > minWlepPt_;
}

/////////Check WZ Properties/////
inline bool WZAnalyzer::passValidWZCut() const{
  return wzCand_ && wzCand_.mass("minPz")>0.;
}

////////////////////////////////
/////////Check Electron Properties/////
////////////////////////////////
bool WZAnalyzer::passTriggerEmulation(const heep::Ele& elec, const float minPt) const{
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
inline bool WZAnalyzer::passHtCut() const{
//-----------------------------------------------------------
  if(debugme) cout<<"Check Ht Cuts"<<endl;
  return Ht_ > minHt_;   
}//--- passHtCut

///////////////////////////////////
//Fake Rate Cuts
inline bool WZAnalyzer::passWFlavorElecCut() const{
  return wCand_ && wCand_.flavor() == PDGELEC;
}

inline bool WZAnalyzer::passWFlavorMuonCut() const{
  return wCand_ && wCand_.flavor() == PDGMUON;
}

bool WZAnalyzer::passFakeEvtCut() const{
  if(looseElectrons_.size() != 1) return false;
  if(looseMuons_    .size() != 1) return false;
  if(looseMuons_[0].charge() != looseElectrons_[0].charge()) return false;
  return true;
}

//pass Fake Lepton Tag Cut
//-----------------------------------------------------------
bool WZAnalyzer::passFakeLeptonTagCut() const{
  if(wCand_.flavor() == PDGELEC){
    return WPrimeUtil::Contains(looseElectrons_[0].patEle(), tightElectrons_);
  }else if(wCand_.flavor() == PDGMUON){
    return WPrimeUtil::Contains(looseMuons_[0],tightMuons_);
  }
  return true;
}//--- Tag Cut

bool WZAnalyzer::passFakeLeptonProbeTightCut() const{
  if(wCand_.flavor() == PDGELEC){
    return WPrimeUtil::Contains(looseMuons_[0], tightMuons_); //Check the other lepton
  }else if(wCand_.flavor() == PDGMUON){
    return WPrimeUtil::Contains(looseElectrons_[0].patEle(), tightElectrons_);
  }
  return true;
}//--- Probe Cut

///////////////////////////////////

//calc Ht
//-----------------------------------------------------------
inline float WZAnalyzer::calcHt() const{
  return WLepPt() + ZLepPt(0) + ZLepPt(1);
}//--- calcHt

inline float WZAnalyzer::calcTriLepMass() const{
  return (zCand_.daughter(0)->p4() +
          zCand_.daughter(1)->p4() + 
          wCand_.daughter(0)->p4()).M();
}//--- calcTriLepMass

inline float WZAnalyzer::calcQ() const{
  return wzCand_.mass("minPz") - zCand_.mass() - WMASS;
}

inline int WZAnalyzer::calcEvtType() const{
  return (zCand_ && wCand_) ?  2 * (zCand_.flavor() != 11) + (wCand_.flavor() != 11) : -999;
}

inline bool WZAnalyzer::inEE(const TeVMuon& mu) const{
  return fabs(mu.eta()) >= 1.05;
}

///////////////Utilities//////////////////
//--------------------------------------------------------------

inline void
WZAnalyzer::clearEvtVariables(){
  AnalyzerBase::clearEvtVariables();
  met_ = pat::MET();
  wzCand_ = XWLeptonic();
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
  runNumber_ = 0;
  lumiNumber_ = 0;
  evtNumber_ = 0;
  MET_ = 0;
  weight_ = 0;
}

float
WZAnalyzer::WLepPt() const{
  if(wCand_.flavor() == PDGELEC){
    return WPrimeUtil::Find(*wCand_.daughter(0), allElectrons_).patEle().pt();
  }else if(wCand_.flavor() == PDGMUON){
    return WPrimeUtil::Find(*wCand_.daughter(0), allMuons_).pt();
  }
  return -999.;
}

inline float
WZAnalyzer::ZLepPt(int idx) const{
  if(zCand_.flavor() == PDGELEC)
    return WPrimeUtil::Find(*zCand_.daughter(idx), allElectrons_).patEle().pt();
  else if(zCand_.flavor() == PDGMUON)
    return WPrimeUtil::Find(*zCand_.daughter(idx), allMuons_).pt();
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
