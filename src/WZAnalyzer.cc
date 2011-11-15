#include "UserCode/CMGWPrimeGroup/interface/WZAnalyzer.h"

using namespace std;

WZAnalyzer::WZAnalyzer(){}
WZAnalyzer::WZAnalyzer(const edm::ParameterSet & cfg, int fileToRun) :
  AnalyzerBase(cfg, fileToRun){
  setupCutOrder();
  if(debugme) printf("Using %i cuts\n",NCuts_);

  wzAlgo_ = (NuAlgos)cfg.getUntrackedParameter<int>("NuAlgo", kMinPz);
  doSystematics_ = cfg.getUntrackedParameter<bool>("doSystematics", false);

  effectiveElecArea_ = cfg.getParameter<vector<double> >("effectiveElecArea");
  effectiveMuonArea_ = cfg.getParameter<vector<double> >("effectiveMuonArea");
  
// +++++++++++++++++++General Cut values
  maxNumZs_ = cfg.getParameter<uint>("maxNumZs");
  minLeadPt_ = cfg.getParameter<double>("minLeadPt");
  minNLeptons_ = cfg.getUntrackedParameter<uint>("minNLeptons", 0);
 
  minMET_ = cfg.getUntrackedParameter<double>("minMET", 0.);

// +++++++++++++++++++Ht Cuts
  minHt_ = cfg.getUntrackedParameter<double>("minHt", -1);

// +++++++++++++++++++W Cuts
  minWlepPt_ = cfg.getParameter<double>("minWlepPt");
  minWtransMass_ = cfg.getUntrackedParameter<double>("minWtransMass", 0.);
  minWpt_ = cfg.getUntrackedParameter<double>("minWpt", 0.);

// +++++++++++++++++++Z Cuts
  minZeePt1_ = cfg.getParameter<double>("minZeePt1");
  minZeePt2_ = cfg.getParameter<double>("minZeePt2");
  minZmmPt1_ = cfg.getParameter<double>("minZmmPt1");
  minZmmPt2_ = cfg.getParameter<double>("minZmmPt2");
  minZmass_ = cfg.getUntrackedParameter<double>("minZmass", 0.);
  maxZmass_ = cfg.getUntrackedParameter<double>("maxZmass", 9e9);
  minZpt_ = cfg.getUntrackedParameter<double>("minZpt", 0.);

  maxZMassDiff_ = cfg.getParameter<double>("maxZMassDiff");
  minDeltaR_ = cfg.getParameter<double>("minDeltaR");


  //Selectors
  Pset eSelectorPset = cfg.getParameter<Pset>("electronSelectors");
  string looseElectronType = cfg.getUntrackedParameter<string>("LooseElectronType", "wp95");
  string tightElectronType = cfg.getUntrackedParameter<string>("TightElectronType", "wp95");
  looseElectron_ = ElectronSelector(eSelectorPset, looseElectronType);
  tightElectron_ = ElectronSelector(eSelectorPset, tightElectronType);
  if(debugme) cout<<"Using "<<looseElectronType<<" for loose electrons and "
                  <<tightElectronType<<" for tight electrons\n";

  Pset mSelectorPset = cfg.getParameter<Pset>("muonSelectors");
  string looseMuonType = cfg.getUntrackedParameter<string>("LooseMuonType", "exotica");
  string tightMuonType = cfg.getUntrackedParameter<string>("TightMuonType", "exotica");
  looseMuon_ = MuonSelector(mSelectorPset, looseMuonType);
  tightMuon_ = MuonSelector(mSelectorPset, tightMuonType);
  if(debugme) cout<<"Using "<<looseMuonType<<" for loose muons and "
                  <<tightMuonType<<" for tight muons\n";

}

WZAnalyzer::~WZAnalyzer(){
}

void WZAnalyzer::setupCutOrder(){
  map<string,fnCut > mFnPtrs;
  mFnPtrs["NoCuts"] = boost::bind(&WZAnalyzer::passNoCut, this);
  mFnPtrs["HLT"] = boost::bind(&WZAnalyzer::passTriggersCut, this);
  mFnPtrs["MinNLeptons"] = boost::bind(&WZAnalyzer::passMinNLeptonsCut, this, boost::cref(looseElectrons_), boost::cref(looseMuons_), boost::cref(minNLeptons_));
  mFnPtrs["ValidW"] = boost::bind(&WZAnalyzer::passValidWCut, this, boost::ref(wCand_));
  mFnPtrs["ValidZ"] = boost::bind(&WZAnalyzer::passValidZCut, this, boost::ref(zCand_));
  mFnPtrs["ValidWZCand"] = boost::bind(&WZAnalyzer::passValidWZCut, this, boost::ref(wzCand_));
  mFnPtrs["LeadLepPt"] = boost::bind(&WZAnalyzer::passLeadingLeptonPtCut, this);
  mFnPtrs["NumZs"] = boost::bind(&WZAnalyzer::passNumberOfZsCut, this);
  mFnPtrs["ZMass"] = boost::bind(&WZAnalyzer::passZMassCut, this, boost::cref(zCand_), boost::cref(minZmass_), boost::cref(maxZmass_));
  mFnPtrs["WTransMass"] = boost::bind(&WZAnalyzer::passWtransMassCut, this, boost::cref(wCand_), boost::cref(minWtransMass_));
  mFnPtrs["MET"] = boost::bind(&WZAnalyzer::passMinMETCut, this, boost::cref(met_), boost::cref(minMET_));
  mFnPtrs["Ht"] = boost::bind(&WZAnalyzer::passHtCut, this);
  mFnPtrs["Zpt"] = boost::bind(&WZAnalyzer::passZptCut, this, boost::cref(zCand_), boost::cref(minZpt_));
  mFnPtrs["Wpt"] = boost::bind(&WZAnalyzer::passWptCut, this, boost::cref(wCand_), boost::cref(minWpt_));
  mFnPtrs["AllCuts"] = boost::bind(&WZAnalyzer::passNoCut, this);

  mFnPtrs["WLepTight"] = boost::bind(&WZAnalyzer::passWLepTightCut, this);
  mFnPtrs["WFlavorElec"] = boost::bind(&WZAnalyzer::passWFlavorElecCut, this);
  mFnPtrs["WFlavorMuon"] = boost::bind(&WZAnalyzer::passWFlavorMuonCut, this);
  mFnPtrs["FakeEvt"] = boost::bind(&WZAnalyzer::passFakeEvtCut, this);
  mFnPtrs["FakeLepProbe"] = boost::bind(&WZAnalyzer::passFakeLeptonProbeCut, this);

  fillCuts(mFnPtrs);

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

  if(doSystematics_){//These histograms are only for systematic studies (should put in a flag)
    defineHistoset("hZeeMassTT","Reconstructed MassTT of ZeeTT",
                   "M_{Z}^{ee,TT} (GeV)", 30, 60, 120, "GeV", hZeeMassTT,dir);
    defineHistoset("hZeeMassTF","Reconstructed Mass of ZeeTF",
                   "M_{Z}^{ee,TF} (GeV)", 30, 60, 120, "GeV", hZeeMassTF,dir);
    defineHistoset("hZmmMassTT","Reconstructed Mass of Z#mu#muTT",
                   "M_{Z}^{#mu#mu,TT} (GeV)", 30, 60, 120, "GeV", hZmmMassTT,dir);
    defineHistoset("hZmmMassTF","Reconstructed Mass of Z#mu#muTF",
                   "M_{Z}^{#mu#mu,TF} (GeV)", 30, 60, 120, "GeV", hZmmMassTF,dir);
  }

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
  defineHistoset("hMETSig", "MET Significance",
                  "#slash{E}_{T}^{Signif}", 50, 0, 500, "NONE", hMETSig,dir);

//W Trans Mass Histos
  defineHistoset("hWTransMass", "Reconstructed Transverse Mass of W",
                  "M_{T} (GeV)", 20, 0, 100, "GeV", hWTransMass,dir);
  defineHistoset("hWenuTransMass", "Reconstructed Transverse Mass of We\\nu",
                  "M_{T}^{e#nu} (GeV)", 20, 0, 100, "GeV", hWenuTransMass,dir);
  defineHistoset("hWmnuTransMass", "Reconstructed TransverseMass of W#mu\\nu",
                  "M_{T}^{#mu#nu} (GeV)", 20, 0, 100, "GeV", hWmnuTransMass,dir);

//Wpt Histos
  defineHistoset("hWpt", "p_{T}^{W}", 
                  "p_{T}^{W} (GeV)", 50, 0, 1000, "GeV", hWpt,dir);

//W Charge Histos
  defineHistoset("hWQ", "Reconstructed Charge of W",
                  "q_{W}", 3, -1, 1, "", hWQ,dir);

  defineHistoset("hNLElec", "Number of Loose Electrons in Event",
                  "N_{e}^{Loose}", 10, 0, 10, "NONE", hNLElec,dir);
  defineHistoset("hNLMuon", "Number of Loose Muons in Event",
                  "N_{#mu}^{Loose}", 10, 0, 10, "NONE", hNLMuon,dir);
  defineHistoset("hNLLeps", "Number of Loose Leptons in Event",
                  "N_{Lep}^{Loose}", 10, 0, 10, "NONE", hNLLeps,dir);

  defineHistoset("hNTElec", "Number of Tight Electrons in Event",
                  "N_{e}^{Tight}", 10, 0, 10, "NONE", hNTElec,dir);
  defineHistoset("hNTMuon", "Number of Tight Muons in Event",
                  "N_{#mu}^{Tight}", 10, 0, 10, "NONE", hNTMuon,dir);
  defineHistoset("hNTLeps", "Number of Tight Leptons in Event",
                  "N_{Lep}^{Tight}", 10, 0, 10, "NONE", hNTLeps,dir);

  defineHistoset("hNJets", "Number of Jets in Event",
                  "N_{Jets}", 15, 0, 15, "NONE", hNJets,dir);

  defineHistoset("hNVtxs", "Number of Vertexs in Event",
                  "N_{Vtx}", 50, 0, 50, "NONE", hNVtxs,dir);

  defineHistoset("hWenuCombRelIso", "Comb Rel Iso of W Electron",
                  "Electron Combined Relative Isolation", 20, 0, 0.2, "NONE", hWenuCombRelIso,dir);
  defineHistoset("hWmnuCombRelIso", "Comb Rel Iso of W Muon",
                  "Muon Combined Relative Isolation", 20, 0, 0.2, "NONE", hWmnuCombRelIso,dir);

  hVtxMatch = dir.make<TH1F>("hVtxMatch","Mask of leptons in PV", 10, 0, 10);

  tWZCand = dir.make<TTree>("tWZCand", "Analysis Variables after WZCand");//Only 1 for now;
  tWZCand->Branch("Run", &runNumber_);
  tWZCand->Branch("Lumi", &lumiNumber_);
  tWZCand->Branch("Event", &evtNumber_);
  tWZCand->Branch("WZMass", &WZMass_);
  tWZCand->Branch("EvtType", &evtType_);
  tWZCand->Branch("Ht", &Ht_);
  tWZCand->Branch("Zpt", &Zpt_);
  tWZCand->Branch("ZMass", &ZMass_);
  tWZCand->Branch("Wpt", &Wpt_);
  tWZCand->Branch("WTransMass", &WTransMass_);
  tWZCand->Branch("MET", &MET_);
  tWZCand->Branch("METSig", &METSig_);
  tWZCand->Branch("Q", &Q_);
  tWZCand->Branch("TriLepMass", &TriLepMass_);
  tWZCand->Branch("NVtxs", &NVtxs_);
  tWZCand->Branch("weight", &weight_);

}//defineHistos

//fill Histograms
void WZAnalyzer::fillHistos(const int& index, const float& weight){
  if(debugme) printf("filling Histos\n");
  if(wCand_ && zCand_){
    hWZMass[index]->Fill(WZMass_, weight);
    if     (evtType_ == 0) hWZ3e0muMass[index]->Fill(WZMass_, weight);
    else if(evtType_ == 1) hWZ2e1muMass[index]->Fill(WZMass_, weight);
    else if(evtType_ == 2) hWZ1e2muMass[index]->Fill(WZMass_, weight);
    else if(evtType_ == 3) hWZ0e3muMass[index]->Fill(WZMass_, weight);
    hQ[index]->Fill(Q_, weight); 
    hWZTransMass[index]->Fill(wzCand_().mt(), weight);
    hWZpt[index]->Fill(wzCand_().pt(), weight);
    hHt[index]->Fill(Ht_, weight);
    hTriLepMass[index]->Fill(TriLepMass_, weight);
    hEvtType[index]->Fill(evtType_, weight);
    if     (wCand_.charge() > 0) hEvtTypeP[index]->Fill(evtType_, weight);
    else if(wCand_.charge() < 0) hEvtTypeM[index]->Fill(evtType_, weight);
    hLeadPt[index]->Fill(LeadPt_, weight);
    hLeadElecPt[index]->Fill(LeadElecPt_, weight);
    hLeadMuonPt[index]->Fill(LeadMuonPt_, weight);
    if     (zCand_.flavor() == PDG_ID_ELEC){
      hLeadPtZee[index]->Fill(LeadPt_, weight);
    }else if(zCand_.flavor() == PDG_ID_MUON){ 
      hLeadPtZmm[index]->Fill(LeadPt_, weight);
    }
    if(CutNames_[index] == "ValidWZCand"){//All, Wpt, Zpt, Ht + 1 for starting @ 0
      tWZCand->Fill();
    }
  }
  if(zCand_){
    hZMass[index]->Fill(zCand_.mass(), weight);
    hZpt[index]->Fill(zCand_.pt(), weight);
    if      (zCand_.flavor() == PDG_ID_ELEC){
      hZeeMass[index]->Fill(zCand_.mass(), weight);
      if(doSystematics_){
        if(TT) hZeeMassTT[index]->Fill(zCand_.mass(), weight);
        if(TF) hZeeMassTF[index]->Fill(zCand_.mass(), weight);
      }
      hZeept[index]->Fill(zCand_.pt(), weight);
      hMETee[index]->Fill(met_.et(), weight);
    }else if (zCand_.flavor() == PDG_ID_MUON){
      hZmmMass[index]->Fill(zCand_.mass(), weight);
      if(doSystematics_){
        if(TT) hZmmMassTT[index]->Fill(zCand_.mass(), weight);
        if(TF) hZmmMassTF[index]->Fill(zCand_.mass(), weight);
      }
      hMETmm[index]->Fill(met_.et(), weight);
      hZmmpt[index]->Fill(zCand_.pt(), weight);
    }
  }
  if(wCand_){
    hWTransMass[index]->Fill(wCand_.mt(), weight);
    hWpt[index]->Fill(wCand_.pt(), weight);
    hWQ[index]->Fill(wCand_.charge(), weight);
    ///////////hWTheta[index]->Fill(wCand_.theta(), weight);
    if      (wCand_.flavor() == PDG_ID_ELEC){
      hWenuTransMass[index]->Fill(wCand_.mt(), weight);
      const heep::Ele& e = *wCand_.elec();
      hWenuCombRelIso[index]->Fill(calcCombRelIso(e.patEle(), ElecPU(e)), weight);
    }else if (wCand_.flavor() == PDG_ID_MUON){
      hWmnuTransMass[index]->Fill(wCand_.mt(), weight);
      const TeVMuon& m = *wCand_.muon();
      hWmnuCombRelIso[index]->Fill(m.combRelIsolation03(MuonPU(m)), weight);
    }
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
    if(passTriggerMatch(looseElectrons_[i], matchptcut, triggersToUse_))
      zElectrons.push_back(looseElectrons_[i]);
  matchptcut = 20.;
  minHighPt = false;
  for (size_t i=0; i < zElectrons.size(); i++){
    if(passTriggerMatch(looseElectrons_[i], matchptcut, triggersToUse_)){
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

    if(doSystematics_){
      bool tight1=false, tight2=false;
      if(zCand_.flavor() == PDG_ID_ELEC){
        for(uint i=0; i<tightElectrons_.size(); ++i){
          if(!tight1 && WPrimeUtil::Match(tightElectrons_[i], *zCand_.daughter(0))) tight1 = true;
          if(!tight2 && WPrimeUtil::Match(tightElectrons_[i], *zCand_.daughter(1))) tight2 = true;
        }
      }else if(zCand_.flavor() == PDG_ID_MUON){
        for(uint i=0; i<tightMuons_.size(); ++i){
          if(!tight1 && WPrimeUtil::Match(tightMuons_[i], *zCand_.daughter(0))) tight1 = true;
          if(!tight2 && WPrimeUtil::Match(tightMuons_[i], *zCand_.daughter(1))) tight2 = true;
        } 
      }
      TT = tight1 && tight2;
      TF = (tight1 && !tight2) || (!tight1 && tight2);
    }
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
WZAnalyzer::calcWZVariables(){
  if (debugme) cout<<"In calc WZ Variables\n";
  wzCand_ = (zCand_ && wCand_) ? XWLeptonic(zCand_, wCand_) : XWLeptonic();
  WZMass_ = wzCand_(wzAlgo_).mass();
  Q_ = (zCand_ && wCand_) ? calcQ() : -999.;
  if(debugme) printEventDetails();
}

void
WZAnalyzer::calcEventVariables(){
  if (debugme) cout<<"In calc Event Variables\n";
  evtType_ = (zCand_ && wCand_) ? calcEvtType() : -999;
  if(debugme) printf("evt Type: %i, Z Flav: %i, W Flav: %i\n", evtType_, (int)zCand_.flavor(), (int)wCand_.flavor());
  LeadPt_ = calcLeadPt(); 
  LeadElecPt_ = calcLeadPt(PDG_ID_ELEC);
  LeadMuonPt_ = calcLeadPt(PDG_ID_MUON);
  Ht_ = (zCand_ && wCand_) ? calcHt() : -999.;
  TriLepMass_ = (zCand_ && wCand_) ? calcTriLepMass() : -999.;
  ZMass_ = zCand_.mass();
  WTransMass_ = wCand_.mt();
  MET_ = met_.et();
  METSig_ = met_.significance();
  NVtxs_ = vertices_.size();
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
    if (looseElectron_(allElectrons_[i].patEle(), pu))
      looseElectrons_.push_back(allElectrons_[i]);

    if (tightElectron_(allElectrons_[i].patEle(), pu))
      tightElectrons_.push_back(allElectrons_[i]);
  }

  for (size_t i = 0; i < allMuons_.size(); i++) {
    const float pu = MuonPU(allMuons_[i]);
    if (looseMuon_(allMuons_[i], pu))
      looseMuons_.push_back(allMuons_[i]);

    if (tightMuon_(allMuons_[i], pu))
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
        if (abs(genParticles[i].pdgId()) == PDG_ID_Z){
          genZ = & genParticles[i];
          cout<<"Mass of gen Z is "<<genZ->mass()<<endl;
        }else if (abs(genParticles[i].pdgId()) == PDG_ID_W){ 
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

  /*
  const reco::Vertex& pv = findPV(vertices_, *zCand_.daughter(0));

  int vtxMask = 0;
  if(!sameVertex(pv, *zCand_.daughter(0))) vtxMask += 1;
  if(!sameVertex(pv, *zCand_.daughter(1))) vtxMask += 2;
  if(!sameVertex(pv, *wCand_.daughter(0))) vtxMask += 4;
  */
  /*
  float cutoff = 0.03;
  int vtxMask = 0;
  if(fabs(zCand_.daughter(0)->vz() - zCand_.daughter(0)->vz()) > cutoff) vtxMask += 1;
  if(fabs(zCand_.daughter(1)->vz() - zCand_.daughter(0)->vz()) > cutoff) vtxMask += 2;
  if(fabs(wCand_.daughter(0)->vz() - zCand_.daughter(0)->vz()) > cutoff) vtxMask += 4;
  
  if(vtxMask) cout<<"Someone is not a match to the pv! = "<<vtxMask<<" and z mass "<<zCand_.mass()<<endl;
  hVtxMatch->Fill(vtxMask, weight_);
  */
  if(0 && !wprimeUtil_->runningOnData()){//Don't do this for data
    GenParticleV genParticles = getProduct<GenParticleV>(event, "genParticles");
    const reco::Candidate * genE = 0;
    for (size_t i = 0; i < genParticles.size(); i++){
      if (abs(genParticles[i].pdgId()) == PDG_ID_ELEC){
        genE = & genParticles[i];
        cout<<"Gen e : "
            <<" pt: "<<genE->pt()
            <<" eta: "<<genE->eta()
            <<" phi: "<<genE->phi()
            <<endl;
      }
    }
  }
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
  if(zCand_ && wCand_ && wzCand_(wzAlgo_).mass()>0.){
    cout<<" WZ Mass: "<<wzCand_(wzAlgo_).mass()
        <<" Neu Pz: "<<wzCand_.neutrinoPz(wzAlgo_)
        <<" Ht: "<<Ht_
        <<" Zpt: "<<zCand_.pt()
        <<" Wpt: "<<wCand_.pt()
        <<endl;
  }
  return;
}

void
WZAnalyzer::printEventLeptons() const{
  if     (zCand_.flavor() == PDG_ID_ELEC){
    cout<<"------- Electron 1 from Z -------\n";
//    printElectron(*zCand_.elec1(), PDG_ID_Z);
    printElectron(WPrimeUtil::Find(*zCand_.daughter(0), allElectrons_));
    cout<<"------- Electron 2 from Z -------\n";
//    printElectron(*zCand_.elec2(), PDG_ID_Z);
    printElectron(WPrimeUtil::Find(*zCand_.daughter(1), allElectrons_));
  }else if(zCand_.flavor() == PDG_ID_MUON){
    cout<<"------- Muon 1 from Z -------\n";
//    printMuon(*zCand_.muon1(), PDG_ID_Z);
    printMuon(WPrimeUtil::Find(*zCand_.daughter(0), allMuons_));
    cout<<"------- Muon 2 from Z -------\n";
//    printMuon(*zCand_.muon2(), PDG_ID_Z);
    printMuon(WPrimeUtil::Find(*zCand_.daughter(1), allMuons_));
  }

  if     (wCand_.flavor() == PDG_ID_ELEC){   
    cout<<"------- Electron from W -------\n";
    printElectron(WPrimeUtil::Find(*wCand_.daughter(0), allElectrons_));
//    printElectron(*wCand_.elec(), PDG_ID_W);
  }else if(wCand_.flavor() == PDG_ID_MUON){
    cout<<"------- Muon from W -------\n";
    printMuon(WPrimeUtil::Find(*wCand_.daughter(0), allMuons_));
//    printMuon    (*wCand_.muon(), PDG_ID_W);
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
  return TMath::Max(calcLeadPt(PDG_ID_ELEC), calcLeadPt(PDG_ID_MUON));
}

/////////////////Accessors///////////////////////

/////////////////Modifiers///////////////////////

/////////////////Cuts///////////////////////
inline bool WZAnalyzer::passWLepTightCut() const{
  if(wCand_.flavor() == PDG_ID_ELEC){
    const heep::Ele & e = *wCand_.elec();
    return WPrimeUtil::Contains(e.patEle(), tightElectrons_);
  }else if(wCand_.flavor() == PDG_ID_MUON){
    const TeVMuon & m = *wCand_.muon();
    return WPrimeUtil::Contains(m, tightMuons_);
}
  return false;
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
  if     (zCand_.flavor() == PDG_ID_ELEC){
    return ( max(ZLepPt(0),ZLepPt(1)) > minZeePt1_ && 
             min(ZLepPt(0),ZLepPt(1)) > minZeePt2_ );
  }else if(zCand_.flavor() == PDG_ID_MUON){

    return ( max(ZLepPt(0),ZLepPt(1)) > minZmmPt1_ && 
             min(ZLepPt(0),ZLepPt(1)) > minZmmPt2_ );
  }
  return false;
}

bool 
WZAnalyzer::passTriggerMatch(const heep::Ele& e, const float cut, const vstring& triggers) const{
  const pat::Electron& p = e.patEle();
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


////////////////////////////////
/////////Check W Properties/////
////////////////////////////////
inline bool
WZAnalyzer::passValidWCut(WCandidate& w){
  calcWVariables();//Cory: These should take args
  calcEventVariables();  
  return AnalyzerBase::passValidWCut(w);
}

//Check W Transverse Mass
inline bool
WZAnalyzer::passWLepPtCut() const{
  return WLepPt() > minWlepPt_;
}

////////////////////////////////
/////////Check Z Properties/////
////////////////////////////////

inline bool
WZAnalyzer::passValidZCut(ZCandidate& z){
  calcZVariables();//Cory: These should take args
  return AnalyzerBase::passValidZCut(z);
}

/////////Check WZ Properties/////
inline bool 
WZAnalyzer::passValidWZCut(XWLeptonic & xw){
  calcWZVariables();
  return AnalyzerBase::passValidXWCut(xw);
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
  return wCand_ && wCand_.flavor() == PDG_ID_ELEC;
}

inline bool WZAnalyzer::passWFlavorMuonCut() const{
  return wCand_ && wCand_.flavor() == PDG_ID_MUON;
}

bool WZAnalyzer::passFakeEvtCut() const{
  if(looseElectrons_.size() != 1) return false;
  if(looseMuons_    .size() != 1) return false;
  if(looseMuons_[0].charge() != looseElectrons_[0].charge()) return false;
  return true;
}

bool WZAnalyzer::passFakeLeptonProbeCut() const{
  if(wCand_.flavor() == PDG_ID_ELEC){
    return WPrimeUtil::Contains(looseMuons_[0], tightMuons_); //Check the other lepton
  }else if(wCand_.flavor() == PDG_ID_MUON){
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
  return wzCand_(wzAlgo_).mass() - zCand_.mass() - WMASS;
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
  zCand_ = ZCandidate();
  wCand_ = WCandidate();
  wzCand_ = XWLeptonic();
  evtType_ = -999;
  LeadPt_ = -999;
  LeadElecPt_ = -999;
  LeadMuonPt_ = -999;
  WZMass_ = -999;
  Ht_= -999;
  TriLepMass_ = -999;
  Zpt_ = -999;
  ZMass_ = -999;
  Wpt_ = -999;
  WTransMass_ = -999;
  Q_ = -999;
  TT = TF = false;
  runNumber_ = 0;
  lumiNumber_ = 0;
  evtNumber_ = 0;
  MET_ = 0;
  METSig_ = 0;
  NVtxs_ = 0;
  weight_ = 0;
}

float
WZAnalyzer::WLepPt() const{
  if(wCand_.flavor() == PDG_ID_ELEC){
    return WPrimeUtil::Find(*wCand_.daughter(0), allElectrons_).patEle().pt();
  }else if(wCand_.flavor() == PDG_ID_MUON){
    return WPrimeUtil::Find(*wCand_.daughter(0), allMuons_).pt();
  }
  return -999.;
}

inline float
WZAnalyzer::ZLepPt(int idx) const{
  if(zCand_.flavor() == PDG_ID_ELEC)
    return WPrimeUtil::Find(*zCand_.daughter(idx), allElectrons_).patEle().pt();
  else if(zCand_.flavor() == PDG_ID_MUON)
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

