#include "UserCode/CMGWPrimeGroup/interface/WZAnalyzer.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
using namespace std;
WZAnalyzer::WZAnalyzer(){}
WZAnalyzer::WZAnalyzer(const edm::ParameterSet & cfg, WPrimeUtil * wprimeUtil){

  wprimeUtil_ = wprimeUtil;
  assert(wprimeUtil_);

   Cuts_          = cfg.getParameter<vstring>("Cuts");
  NCuts_          = Cuts_.size();
   LooseElecCuts_ = cfg.getParameter<vstring>("LooseElecCuts");
  NLooseElecCuts_ = LooseElecCuts_.size();
   LooseMuonCuts_ = cfg.getParameter<vstring>("LooseMuonCuts");
  NLooseMuonCuts_ = LooseMuonCuts_.size();
   TightElecCuts_ = cfg.getParameter<vstring>("TightElecCuts");
  NTightElecCuts_ = TightElecCuts_.size();
   TightMuonCuts_ = cfg.getParameter<vstring>("TightMuonCuts");
  NTightMuonCuts_ = TightMuonCuts_.size();

  FillCutFns();
  SetLogFile(cfg.getParameter<string>("LogFile"));
  SetCandEvtFile(cfg.getParameter<string>("CandEvtFile"));

  intOptions_["report"] = cfg.getParameter<uint>("reportAfter");
  intOptions_["verbose"] = cfg.getParameter<bool>("debugme");
  doPreselect_ = cfg.getParameter<bool>("preselect");
  intOptions_["events"] = 0;

  debugme = cfg.getParameter<bool>("debugme");

  electronsLabel_ = cfg.getParameter<string>("electrons");
  muonsLabel_ = cfg.getParameter<string>("muons");
  metLabel_ = cfg.getParameter<string>("met");

  muonAlgo_ = cfg.getParameter<uint>("muonAlgo");
  
  hltEventLabel_ = cfg.getParameter<string>("hltEventTag");
  pileupLabel_ = cfg.getParameter<string>("pileupTag");

  triggersToUse_          = cfg.getParameter<vstring>("triggersToUse");

  minDeltaR_ = cfg.getParameter<double>("minDeltaR");

  ResetCounters();

  verbose("Using %i cuts\n",NCuts_);
  verbose("Using %i loose elec cuts\n",NLooseElecCuts_);
  verbose("Using %i loose muon cuts\n",NLooseMuonCuts_);
  verbose("Using %i tight elec cuts\n",NTightElecCuts_);
  verbose("Using %i tight muon cuts\n",NTightMuonCuts_);

  PDGMUON = 13;
  PDGELEC = 11;
  PDGW = 24;
  PDGZ = 23;
  PDGWPRIME = 34;

  PI    = 2.0 * TMath::ACos(0.);
  TWOPI = 2.0 * PI;
  NOCUT = 9e9;
 
// +++++++++++++++++++General Cut values
  maxNumZs_ = cfg.getParameter<uint>("maxNumZs");
  minNLeptons_ = cfg.getParameter<uint>("minNLeptons");
  minLeadPt_ = cfg.getParameter<double>("minLeadPt");
  minMET_ = cfg.getParameter<double>("minMET");

// +++++++++++++++++++Ht Cuts
  minHt_ = cfg.getParameter<double>("minHt");

// +++++++++++++++++++W Cuts
  minWtransMass_ = cfg.getParameter<double>("minWtransMass");
  minWpt_ = cfg.getParameter<double>("minWpt");

  maxWmunuCombRelIso_ = cfg.getParameter<double>("maxWmunuCombRelIso");
  cutWenuWPRelIsoMask_ = cfg.getParameter<int>("cutWenuWPRelIsoMask");
  cutElecWPTightType_ = cfg.getParameter<string>("cutElecWPTightType");

// +++++++++++++++++++Z Cuts
  minZpt_ = cfg.getParameter<double>("minZpt");
  minZmass_ = cfg.getParameter<double>("minZmass");
  maxZmass_ = cfg.getParameter<double>("maxZmass");

// +++++++++++++++++++Electron General Cuts
//VBTF Recommended Cuts
  minElecLooseEt_ = cfg.getParameter<double>("minElecLooseEt");
  minElecTightEt_ = cfg.getParameter<double>("minElecTightEt");
  cutElecWPLooseMask_ = cfg.getParameter<int>("cutElecWPLooseMask");
  cutElecWPLooseType_ = cfg.getParameter<string>("cutElecWPLooseType");

  maxElecSigmaiEtaiEta_ = cfg.getParameter<vector<double> >("maxElecSigmaiEtaiEta");
  maxElecDeltaPhiIn_ = cfg.getParameter<vector<double> >("maxElecDeltaPhiIn");
  maxElecDeltaEtaIn_ = cfg.getParameter<vector<double> >("maxElecDeltaEtaIn");
  maxElecHOverE_     = cfg.getParameter<vector<double> >("maxElecHOverE");

// +++++++++++++++++++Muon General Cuts
  maxMuonEta_ = cfg.getParameter<double>("maxMuonEta");
  minMuonLoosePt_ = cfg.getParameter<double>("minMuonLoosePt");
  minMuonTightPt_ = cfg.getParameter<double>("minMuonTightPt");
//VBTF Recommended Cuts
  maxMuonDxy_ = cfg.getParameter<double>("maxMuonDxy");
  maxMuonNormChi2_ = cfg.getParameter<double>("maxMuonNormChi2");
  minMuonNPixHit_ = cfg.getParameter<int>("minMuonNPixHit");
  minMuonNTrkHit_ = cfg.getParameter<int>("minMuonNTrkHit");
  minMuonStations_ = cfg.getParameter<int>("minMuonStations");
  minMuonHitsUsed_ = cfg.getParameter<int>("minMuonHitsUsed");
}
WZAnalyzer::~WZAnalyzer(){
  outCandEvt.close(); 
  outLogFile.close(); 
}

void WZAnalyzer::FillCutFns(){
  mFnPtrs_["NoCuts"] = &WZAnalyzer::PassNoCut;
  mFnPtrs_["HLT"] = &WZAnalyzer::PassTriggersCut;
  mFnPtrs_["NLeptons"] = &WZAnalyzer::PassNLeptonsCut;
  mFnPtrs_["ValidW"] = &WZAnalyzer::PassValidWCut;
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
  mFnPtrs_["AllCuts"] = &WZAnalyzer::PassNoCut;

  //Electron cuts
  mElecFnPtrs_["ElecLoose"] = &WZAnalyzer::PassElecLooseCut;
  mElecFnPtrs_["ElecTight"] = &WZAnalyzer::PassElecTightCut;
  mElecFnPtrs_["ElecLooseEt"] = &WZAnalyzer::PassElecLooseEtCut;
  mElecFnPtrs_["ElecTightEt"] = &WZAnalyzer::PassElecTightEtCut;
  mElecFnPtrs_["ElecLooseID"] = &WZAnalyzer::PassElecLooseWPCut;
  mElecFnPtrs_["ElecIso"] = &WZAnalyzer::PassElecWPRelIsoCut;

  //Muon cuts
  mMuonFnPtrs_["MuonLoose"] = &WZAnalyzer::PassMuonLooseCut;
  mMuonFnPtrs_["MuonTight"] = &WZAnalyzer::PassMuonTightCut;
  mMuonFnPtrs_["MuonLoosePt"] = &WZAnalyzer::PassMuonLoosePtCut;
  mMuonFnPtrs_["MuonTightPt"] = &WZAnalyzer::PassMuonTightPtCut;
  mMuonFnPtrs_["MuonIso"] = &WZAnalyzer::PassMuonCombRelIsoCut;
  mMuonFnPtrs_["MuonEta"] = &WZAnalyzer::PassMuonEtaCut;
  mMuonFnPtrs_["MuonGlobal"] = &WZAnalyzer::PassMuonGlobalCut;
  mMuonFnPtrs_["MuonDxy"] = &WZAnalyzer::PassMuonDxyCut;
  mMuonFnPtrs_["MuonNpxl"] = &WZAnalyzer::PassMuonNpixhitCut;
  mMuonFnPtrs_["MuonNtrk"] = &WZAnalyzer::PassMuonNtrkhitCut;
  mMuonFnPtrs_["MuonNormChi2"] = &WZAnalyzer::PassMuonNormChi2Cut;
  mMuonFnPtrs_["MuonHitsUsed"] = &WZAnalyzer::PassMuonHitsUsedCut;
  mMuonFnPtrs_["MuonStations"] = &WZAnalyzer::PassMuonStationsCut;

  CutFns_.resize(NCuts_);
  for(int i=0; i<NCuts_; ++i){
    CutFns_[i] = mFnPtrs_.find(Cuts_[i])->second;
  }

  LooseElecCutFns_.resize(NLooseElecCuts_);
  for(int i=0; i<NLooseElecCuts_; ++i){
    LooseElecCutFns_[i] = mElecFnPtrs_.find(LooseElecCuts_[i])->second;
  }

  LooseMuonCutFns_.resize(NLooseMuonCuts_);
  for(int i=0; i<NLooseMuonCuts_; ++i){
    LooseMuonCutFns_[i] = mMuonFnPtrs_.find(LooseMuonCuts_[i])->second;
  }

  TightElecCutFns_.resize(NTightElecCuts_);
  for(int i=0; i<NTightElecCuts_; ++i){
    TightElecCutFns_[i] = mElecFnPtrs_.find(TightElecCuts_[i])->second;
  }

  TightMuonCutFns_.resize(NTightMuonCuts_);
  for(int i=0; i<NTightMuonCuts_; ++i){
    TightMuonCutFns_[i] = mMuonFnPtrs_.find(TightMuonCuts_[i])->second;
  }
  
}

void
WZAnalyzer::ResetCounters(){
  verbose("Reset Counters\n");
  Num_surv_cut_.clear();
  Num_surv_cut_.resize(NCuts_,0.);
  eventNum = 0;
  runNumber = -1;
  lumiID = -1;
}

//--------------------------------------------------------------
void WZAnalyzer::Declare_Histos(TFileDirectory & dir)
{
  verbose("Declare histos\n");

  DeclareHistoSet("hEvtType", "Event Type",
                  "N_{#mu}", 4, 0, 4, hEvtType,dir);

  DeclareHistoSet("hWZInvMass", "Reconstructed WZ Invariant Mass",
                  "m_{WZ} (GeV)", 100, 0, 1000, hWZInvMass,dir);
  DeclareHistoSet("hWZ3e0muInvMass", "Reconstructed WZ(3e0\\mu) Invariant Mass",
                  "m_{WZ} (GeV)", 100, 0, 1000, hWZ3e0muInvMass,dir);
  DeclareHistoSet("hWZ2e1muInvMass", "Reconstructed WZ(2e1\\mu) Invariant Mass",
                  "m_{WZ} (GeV)", 100, 0, 1000, hWZ2e1muInvMass,dir);
  DeclareHistoSet("hWZ1e2muInvMass", "Reconstructed WZ(1e2\\mu) Invariant Mass",
                  "m_{WZ} (GeV)", 100, 0, 1000, hWZ1e2muInvMass,dir);
  DeclareHistoSet("hWZ0e3muInvMass", "Reconstructed WZ(0e3\\mu) Invariant Mass",
                  "m_{WZ} (GeV)", 100, 0, 1000, hWZ0e3muInvMass,dir);

  DeclareHistoSet("hWZTransMass", "Reconstructed WZ Transverse Mass",
                  "m_{WZ} (GeV)", 100, 0, 1000, hWZTransMass,dir);
//Ht Histos
  DeclareHistoSet("hHt", "H_{T}", 
                  "Lepton Pt Sum: H_{T} (GeV)", 80, 0, 800, hHt,dir);
//Wpt Histos
  DeclareHistoSet("hWpt", "p_{T}^{W}", 
                  "p_{T}^{W} (GeV)", 40, 0, 400, hWpt,dir);
//Zpt Histos
  DeclareHistoSet("hZpt", "p_{T}^{Z}", 
                  "p_{T}^{Z} (GeV)", 40, 0, 400, hZpt,dir);
//MET Histos
  DeclareHistoSet("hMET", "MET",
                  "MET (GeV)", 30, 0, 300, hMET,dir);

//Z Mass Histos
  DeclareHistoSet("hZMass" , "Reconstructed Mass of Z",
                  "m_{Z}^{Reco} (GeV)", 30, 60, 120, hZMass,dir);
  DeclareHistoSet("hZeeMass","Reconstructed Mass of Zee",
                  "m_{Z}^{Reco} (GeV)", 30, 60, 120, hZeeMass,dir);
  DeclareHistoSet("hZmumuMass","Reconstructed Mass of Zmumu",
                  "m_{Z}^{Reco} (GeV)", 30, 60, 120, hZmumuMass,dir);
  DeclareHistoSet("hZ3e0muMass" , "Reconstructed Mass of Z(3e0\\mu)",
                  "m_{Z}^{Reco} (GeV)", 30, 60, 120, hZ3e0muMass,dir);
  DeclareHistoSet("hZ2e1muMass" , "Reconstructed Mass of Z(2e1\\mu)",
                  "m_{Z}^{Reco} (GeV)", 30, 60, 120, hZ2e1muMass,dir);
  DeclareHistoSet("hZ1e2muMass" , "Reconstructed Mass of Z(1e2\\mu)",
                  "m_{Z}^{Reco} (GeV)", 30, 60, 120, hZ1e2muMass,dir);
  DeclareHistoSet("hZ0e3muMass" , "Reconstructed Mass of Z(0e3\\mu)",
                  "m_{Z}^{Reco} (GeV)", 30, 60, 120, hZ0e3muMass,dir);
  DeclareHistoSet("hZeeMassTT","Reconstructed MassTT of ZeeTT",
                  "m_{Z}^{Reco} (GeV)", 30, 60, 120, hZeeMassTT,dir);
  DeclareHistoSet("hZeeMassTF","Reconstructed Mass of ZeeTF",
                  "m_{Z}^{Reco} (GeV)", 30, 60, 120, hZeeMassTF,dir);
  DeclareHistoSet("hZmumuMassTT","Reconstructed Mass of ZmumuTT",
                  "m_{Z}^{Reco} (GeV)", 30, 60, 120, hZmumuMassTT,dir);
  DeclareHistoSet("hZmumuMassTF","Reconstructed Mass of ZmumuTF",
                  "m_{Z}^{Reco} (GeV)", 30, 60, 120, hZmumuMassTF,dir);

//W Trans Mass Histos
  DeclareHistoSet("hWTransMass", "Reconstructed TransMass of W",
                  "m_{T}^{Reco} (GeV)", 20, 0, 100, hWTransMass,dir);
  DeclareHistoSet("hWenuTransMass", "Reconstructed TransMass of Wenu",
                  "m_{T}^{Reco} (GeV)", 20, 0, 100, hWenuTransMass,dir);
  DeclareHistoSet("hWmunuTransMass", "Reconstructed TransMass of Wmunu",
                  "m_{T}^{Reco} (GeV)", 20, 0, 100, hWmunuTransMass,dir);
  DeclareHistoSet("hW3e0muTransMass", "Reconstructed TransMass of W(3e0\\mu)",
                  "m_{T}^{Reco} (GeV)", 20, 0, 100, hW3e0muTransMass,dir);
  DeclareHistoSet("hW2e1muTransMass", "Reconstructed TransMass of W(2e1\\mu)",
                  "m_{T}^{Reco} (GeV)", 20, 0, 100, hW2e1muTransMass,dir);
  DeclareHistoSet("hW1e2muTransMass", "Reconstructed TransMass of W(1e2\\mu)",
                  "m_{T}^{Reco} (GeV)", 20, 0, 100, hW1e2muTransMass,dir);
  DeclareHistoSet("hW0e3muTransMass", "Reconstructed TransMass of W(0e3\\mu)",
                  "m_{T}^{Reco} (GeV)", 20, 0, 100, hW0e3muTransMass,dir);

//Q=M_{WZ} - M_W - M_Z
  DeclareHistoSet("hQ", "Q=M_{WZ} - M_{W} - M_{Z}",
                  "Q (GeV)", 50, 0, 500, hQ,dir);

///Eff Plots///////
  string title = Form("Expected # of Events / %.0f pb^{-1}",  wprimeUtil_->getLumi_ipb());
  hNumEvts = NULL; hNumEvts = dir.make<TH1F>("hNumEvts",title.c_str(),NCuts_,0,NCuts_);
  hEffRel  = NULL; hEffRel  = dir.make<TH1F>("hEffRel","Relative Efficiency",NCuts_,0,NCuts_);
  hEffAbs  = NULL; hEffAbs  = dir.make<TH1F>("hEffAbs","Absolute Efficiency",NCuts_,0,NCuts_);

//  DeclareHisto("hNumEvts",title.c_str(),       NCuts_,0,NCuts_,hNumEvts,dir);
//  DeclareHisto("hEffRel","Relative Efficiency",NCuts_,0,NCuts_,hEffRel,dir);
//  DeclareHisto("hEffAbs","Absolute Efficiency",NCuts_,0,NCuts_,hEffAbs,dir);

  for(int i=0; i<NCuts_; ++i) hNumEvts->GetXaxis()->SetBinLabel(i+1,Cuts_[i].c_str());
  for(int i=0; i<NCuts_; ++i) hEffRel ->GetXaxis()->SetBinLabel(i+1,Cuts_[i].c_str());
  for(int i=0; i<NCuts_; ++i) hEffAbs ->GetXaxis()->SetBinLabel(i+1,Cuts_[i].c_str());

  listOfHists.push_back(hNumEvts);
  listOfHists.push_back(hEffRel);
  listOfHists.push_back(hEffAbs);

}//Declare_Histos

//Scale Histograms
void
WZAnalyzer::ScaleHistos(){
  verbose("Scaling Histos\n");
  for(uint i=0; i != listOfHists.size(); ++i){
    listOfHists[i]->Scale(wprimeUtil_->getWeight());
  }  
  for(int i=0; i<NCuts_; ++i) Num_surv_cut_[i] *= wprimeUtil_->getWeight();
}

//Fill Histograms
//-----------------------------------------------------------
void WZAnalyzer::Fill_Histos(int index, float weight)
{
//-----------------------------------------------------------
  verbose("Filling Histos\n");
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
  }
  if(zCand_){
    hZpt[index]->Fill(zCand_.pt(), weight);
    hZMass[index]->Fill(zCand_.mass(), weight);
    if      (zCand_.flavor() == PDGELEC){
      hZeeMass[index]->Fill(zCand_.mass(), weight);
    }else if (zCand_.flavor() == PDGMUON){
      hZmumuMass[index]->Fill(zCand_.mass(), weight);
    }
    if     (evtType_ == 0) hZ3e0muMass[index]->Fill(zCand_.mass(), weight);
    else if(evtType_ == 1) hZ2e1muMass[index]->Fill(zCand_.mass(), weight);
    else if(evtType_ == 2) hZ1e2muMass[index]->Fill(zCand_.mass(), weight);
    else if(evtType_ == 3) hZ0e3muMass[index]->Fill(zCand_.mass(), weight);
  }
  if(wCand_){
    hWpt[index]->Fill(wCand_.pt(), weight);
    hWTransMass[index]->Fill(wCand_.mt(), weight);
    if      (wCand_.flavor() == PDGELEC) hWenuTransMass[index]->Fill(wCand_.mt(), weight);
    else if (wCand_.flavor() == PDGMUON) hWmunuTransMass[index]->Fill(wCand_.mt(), weight);
    if(evtType_ == 0) hW3e0muTransMass[index]->Fill(wCand_.mt(), weight);
    if(evtType_ == 1) hW2e1muTransMass[index]->Fill(wCand_.mt(), weight);
    if(evtType_ == 2) hW1e2muTransMass[index]->Fill(wCand_.mt(), weight);
    if(evtType_ == 3) hW0e3muTransMass[index]->Fill(wCand_.mt(), weight);
  }  
  hMET[index]->Fill(met_.et(), weight);
}//Fill_Histos

void
WZAnalyzer::CalcZVariables(){
  if (debugme) cout<<"In Calc Z Variables\n";
  // Reconstruct the Z
  ZCandV looseZCands = getZCands(looseElectrons_, looseMuons_, 12.5);
  zCand_ = looseZCands.size() ? looseZCands[0] : ZCandidate();
  verbose("    Contains: %i loose Z candidate(s)", looseZCands.size());
  numZs_ = looseZCands.size(); 
}

void
WZAnalyzer::CalcWVariables(){
  if (debugme) cout<<"In Calc W Variables\n";
  // Reconstruct the W
  wCand_ = getWCand(tightElectrons_, tightMuons_, met_, zCand_, minDeltaR_);
  verbose("    Contains: %i tight W candidate(s)", (bool)wCand_);
}

void
WZAnalyzer::CalcWZVariables(){
  if (debugme) cout<<"In Calc WZ Variables\n";
  //Calculate Important Quantities for each event
  wzCand_ = (zCand_ && wCand_) ? WZCandidate(zCand_, wCand_) : WZCandidate();
  
  Ht_ = (zCand_ && wCand_) ? Calc_Ht() : -999.;
  Q_ = (zCand_ && wCand_) ? Calc_Q() : -999.;
  verbose("evt Type: %i, Z Flav: %i, W Flav: %i\n", evtType_, zCand_.flavor(), wCand_.flavor());

}

void
WZAnalyzer::CalcEventVariables(){
  if (debugme) cout<<"In Calc Event Variables\n";
  evtType_ = (zCand_ && wCand_) ? Calc_EvtType() : -999;
  LeadPt_ = CalcLeadPt();
  LeadElecPt_ = CalcLeadPt(PDGELEC);
  LeadMuonPt_ = CalcLeadPt(PDGMUON);
//  TT = TF = false;
//  if(zCand_.flavor()){
//    bool tight1 = PassTightCut(zCand_.daughter(0), zCand_.flavor());
//    bool tight2 = PassTightCut(zCand_.daughter(1), zCand_.flavor());
//    //cout<<"tight1: "<<tight1<<" tight2: "<<tight2<<endl;
//    TT = tight1 && tight2;
//    TF = (tight1 && !tight2) || (!tight1 && tight2);
//  }

}

void
WZAnalyzer::DeclareHistoSet(string n, string t, string xtitle,
                            int nbins, float min, float max,
                            vector<TH1F*>& h, TFileDirectory& d){
  ClearAndResize(h,NCuts_,NULL);

  float binWidth = (max-min)/nbins;
  for(int i=0; i<NCuts_; ++i){
    
    string name = n + "_" + Cuts_[i];
    string title = t + " (After " + Cuts_[i] + " Cut);" 
      + xtitle + Form(";Evts/%.0f GeV/%.0f pb^{-1}",binWidth, wprimeUtil_->getLumi_ipb());
//    h[i] = new TH1F(name.c_str(),title.c_str(),nbins,min,max);
    h[i] = d.make<TH1F>(name.c_str(),title.c_str(),nbins,min,max);

//    cout<<"Address in declarehistoset:"<<h[i]<<endl;
//    cout<<"Name in declarehistoset:"<<h[i]->GetName()<<endl;
    listOfHists.push_back(h[i]);
//    cout<<"Name in declarehistoset:"<<listOfHists[listOfHists.size()-1]->GetName()<<endl;
  }
}

//------------------------------------------------------------------------
void WZAnalyzer::deleteHistos()
{
//------------------------------------------------------------------------
  verbose("Delete Histos.....\n");

  for(uint i = 0; i != listOfHists.size(); ++i){
    delete listOfHists[i];
  }  
  listOfHists.clear();
  return;

}//deleteHistos

//Tabulate results after the cut has been passed
//-----------------------------------------------------------
void WZAnalyzer::Tabulate_Me(int& cut_index, const float& weight)
{
//-----------------------------------------------------------
  if(debugme) cout<<"Tabulating results for cut_index = "
                  <<cut_index<<endl;

//increase the number of events passing the cuts
  Num_surv_cut_[cut_index] += weight;
//fill the histograms
  Fill_Histos(cut_index,weight);
    
}//Tabulate_Me

//Writing results to a txt file
//--------------------------------------------------------------------------
void WZAnalyzer::printSummary(const string& dir) 
{ 
//------------------------------------------------------------------------
  if(debugme) cout<<"Writing results to a txt file"<<endl;

  outLogFile << setiosflags(ios::fixed) << setprecision(2);
  outLogFile<<"$$$$$$$$$$$$$$$$$$$$$$$ Type of sample: "<<dir<<endl;
//  outLogFile << " xsec*lumi expected events = " << Nthe_evt << endl;
//  outLogFile << " # of evt passing preselection = " << Nexp_evt << " per "<<Form("%.0f", wprimeUtil_->getLumi_ipb())<<" inv pb"<<endl;
        
  for(int i = 0; i < NCuts_; ++i){
    outLogFile<<right<<"Cut " << setw(2) << i << "("
              <<left<< setw(15) << Cuts_[i]
              <<right << "): " <<"expected evts = " << setw(10) << Num_surv_cut_[i];
    hNumEvts->SetBinContent(i+1,Num_surv_cut_[i]);

//calculate efficiencies
    float eff, deff;
    if(i == 0){
      getEff(eff, deff, Num_surv_cut_[i], Num_surv_cut_[0]);
    }else{
      getEff(eff, deff, Num_surv_cut_[i], Num_surv_cut_[i-1]);
    }
    hEffRel->SetBinContent(i+1,eff*100);
    outLogFile << setw(15) <<"\tRelative eff = "<<setw(6)<<eff*100 << " +/- " << setw(6)<<deff*100 << "%";
    getEff(eff, deff, Num_surv_cut_[i], Num_surv_cut_[0]);
    hEffAbs->SetBinContent(i+1,eff*100);
    outLogFile << setw(15) <<"\tAbsolute eff = "<<setw(6)<<eff*100 << " +/- " << setw(6)<<deff*100 << "%"
               << endl;
        
  } // loop over different cuts
}//printSummary


void 
WZAnalyzer::eventLoop(edm::EventBase const & event){
  //if (intOptions_["events"] && eventNum > intOptions_["events"]) return;
  //reportProgress(eventNum++);
  //updateEventCounts(event, nEvents, runNumber,  wprimeUtil_->getLumi_ipb());
  //datasetName = getDatasetName(event, datasetName);

  ClearEvtVariables();

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
//  verbose("    FlavorType: %i", flavorType);

  // Get leptons
  const vector<pat::Electron> patElectrons = getProduct<vector<pat::Electron> >(event, electronsLabel_);
  const vector<pat::Muon    > patMuons     = getProduct<vector<pat::Muon    > >(event, muonsLabel_);
  verbose("    Contains: %i electron(s), %i muon(s)",
          patElectrons.size(), patMuons.size());

  // Make vectors of leptons passing various criteria
  for (size_t i = 0; i < patElectrons.size(); i++) {
    electrons_.push_back(heep::Ele(patElectrons[i]));   
    if (PassElecLooseCut(&electrons_[i])){
      looseElectrons_.push_back(electrons_[i]);
      if (PassElecTightCut(&electrons_[i]))
        tightElectrons_.push_back(electrons_[i]);
    }
  }
  for (size_t i = 0; i < patMuons.size(); i++) {
    muons_.push_back(TeVMuon(patMuons[i],muonAlgo_));   
    if (PassMuonLooseCut(&muons_[i])){
      looseMuons_.push_back(muons_[i]);
      if (PassMuonTightCut(&muons_[i]))
        tightMuons_.push_back(muons_[i]);
    }
  }

  verbose("    Contains: %i loose electron(s), %i loose muon(s)",
          looseElectrons_.size(), looseMuons_.size());
  verbose("    Contains: %i tight electron(s), %i tightmuon(s)",
          tightElectrons_.size(), tightMuons_.size());

  // Get MET
  met_ = getProduct<METV>(event, metLabel_)[0];
  if(debugme) cout<<"Before met et: "<<met_.et()<<" met phi: "<<met_.phi()<<endl;
  met_ = AdjustedMET(looseElectrons_,looseMuons_, met_);
  if(debugme) cout<<"After  met et: "<<met_.et()<<" met phi: "<<met_.phi()<<endl;

  triggerEvent_ = getProduct<pat::TriggerEvent>(event,hltEventLabel_); 

/*
  GenParticleV genParticles = getUntrackedProduct<GenParticleV>(event, "genParticles");
  const Candidate * genZ = 0;
  const Candidate * genW = 0;
  for (size_t i = 0; i < genParticles.size(); i++)
    if (abs(genParticles[i].pdgId()) == 23) genZ = & genParticles[i];
    else if (abs(genParticles[i].pdgId()) == 24) genW = & genParticles[i];
*/

/*
  if(wprimeUtil_->getSampleName().find("data") == string:npos){
  
    PupInfo_ = getProduct<std::vector< PileupSummaryInfo > >(event, pileupLabel_);   
    
    std::vector<PileupSummaryInfo>::const_iterator PVI;                       
    for(PVI = PupInfo_.begin(); PVI != PupInfo_.end(); ++PVI) {               
      //Cory: What am I looping over? Primary Vtxs
      int PU_BunchCrossing = PVI->getBunchCrossing();                  
      int PU_NumInteractions = PVI->getPU_NumInteractions();           
      cout<<"N_BX: "<<PU_BunchCrossing
          <<" PU_NumInteractions: "<<PVI->getPU_NumInteractions()
          <<endl;
    }
  }
*/
/////////////////////
  if(!PassCuts(wprimeUtil_->getWeight())) return;
  if(!wprimeUtil_->getSampleName().find("data")){
    cout<<" The following data events passed All Cuts!!!\n\n";
    PrintEventFull(event);
    cout<<" ------------------\n";
  }
}

void WZAnalyzer::PrintEventToFile(edm::EventBase const & event){
  outCandEvt<<event.id().run()<<":"
            <<event.id().luminosityBlock()<<":"
            <<event.id().event()<<endl;
}

void WZAnalyzer::PrintEvent(edm::EventBase const & event){
  PrintEventToFile(event);
  cout<<"run #: "<<event.id().run()
      <<" lumi: "<<event.id().luminosityBlock()
      <<" eventID: "<<event.id().event()<<endl

      <<" Z Flavor: "<<zCand_.flavor()
      <<" Z Mass: "<<zCand_.mass()
      <<" W Flavor: "<<wCand_.flavor()
      <<" W MT: "<<wCand_.mt()
      <<endl;

  cout<<" Z lep1 pt "<<zCand_.daughter(0)->pt()
      <<" Z lep2 pt "<<zCand_.daughter(1)->pt()
      <<endl;

  
  cout<<" W lep pt "<<wCand_.daughter(0)->pt()
      <<" pfMet et: "<<met_.et()
      <<endl;

  cout<<" Ht: "<<Ht_
      <<" Zpt: "<<zCand_.pt()
      <<" Wpt: "<<wCand_.pt()<<endl
      <<" WZ Mass: "<<wzCand_.mass("minPz")<<endl
//      <<" # Elec: "<<elec.pt->size()
//      <<" # Muon: "<<mu.pt->size()
      <<endl<<endl;
  return;
}

void
WZAnalyzer::PrintEventFull(edm::EventBase const & event){
  PrintEvent(event);
  PrintTrigger();
  if     (zCand_.flavor() == PDGELEC){
    heep::Ele ze1 = heep::Ele(*(pat::Electron*)zCand_.daughter(0));
    heep::Ele ze2 = heep::Ele(*(pat::Electron*)zCand_.daughter(0));
    PrintElectron(&ze1, PDGZ);
    PrintElectron(&ze2, PDGZ);
//    PrintElectron((const heep::Ele*)zCand_.daughter(0), PDGZ);
//    PrintElectron((const heep::Ele*)zCand_.daughter(1), PDGZ);
  }else if(zCand_.flavor() == PDGMUON){
    PrintMuon((const TeVMuon*)zCand_.daughter(0), PDGZ);
    PrintMuon((const TeVMuon*)zCand_.daughter(1), PDGZ);
  }

  if     (wCand_.flavor() == PDGELEC){
    heep::Ele we = heep::Ele(*(pat::Electron*)wCand_.daughter(0));
    PrintElectron(&we, PDGW);
  }else if(wCand_.flavor() == PDGMUON) PrintMuon    ((const TeVMuon*    )wCand_.daughter(0), PDGW);
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
WZAnalyzer::PrintElectron(const heep::Ele* elec, int parent){
  if     (parent == PDGZ) cout<<"-----Electron from Z-------------------------"<<endl;
  else if(parent == PDGW) cout<<"-----Electron from W-------------------------"<<endl;
  else                    cout<<"-----Electron from ?-------------------------"<<endl;
  cout<<" Elec ScEt: "<<elec->et()<<endl; //ScEt
  if(!elec->isPatEle()) return;
  cout<<" Elec Pt: "<<elec->patEle().pt()<<endl
      <<" Elec Eta: "<<elec->patEle().eta()<<endl //Eta
      <<" Elec SigmaNN: "<<elec->patEle().sigmaIetaIeta()<<endl //sigmaNN
      <<" Elec dPhi: "<<elec->patEle().deltaPhiSuperClusterTrackAtVtx()<<endl //DeltaPhi
      <<" Elec dEta: "<<elec->patEle().deltaEtaSuperClusterTrackAtVtx()<<endl //DeltaEta
      <<" Elec HoverE: "<<elec->patEle().hadronicOverEm()<<endl// H/E
      <<" Elec EoverP: "<<elec->patEle().eSuperClusterOverP()<<endl;// E/P
  cout<<" Elec WP95: "<<elec->patEle().electronID("simpleEleId95relIso")<<endl
      <<" Elec WP90: "<<elec->patEle().electronID("simpleEleId90relIso")<<endl
      <<" Elec WP85: "<<elec->patEle().electronID("simpleEleId85relIso")<<endl
      <<" Elec WP80: "<<elec->patEle().electronID("simpleEleId80relIso")<<endl;
}

void
WZAnalyzer::PrintMuon(const TeVMuon* mu, int parent){
  if     (parent == PDGZ) cout<<"-----Muon from Z-------------------------"<<endl;
  else if(parent == PDGW) cout<<"-----Muon from W-------------------------"<<endl;
  else                    cout<<"-----Muon from ?-------------------------"<<endl;
  reco::TrackRef gm = mu->globalTrack();
  cout<<" Muon Pt: "  <<mu->pt()<<endl
      <<" Muon Eta: " <<mu->eta()<<endl
      <<" Muon Dxy: " <<mu->userFloat("d0")<<endl //Dxy
      <<" Muon NormX2: "<<gm->normalizedChi2()<<endl //NormX2
      <<" Muon NPix: "  <<gm->hitPattern().numberOfValidPixelHits()<<endl //Npixhit
      <<" Muon NTrk: "  <<gm->hitPattern().numberOfValidTrackerHits()<<endl //Ntrk hit
      <<" Muon NMatches: "<<mu->numberOfMatches()<<endl //MuonStations
      <<" Muon Hits: "  <<gm->hitPattern().numberOfValidMuonHits()<<endl; //Muon Hits

  if(parent == PDGW){
    cout<<" Muon RelIso: "<<Calc_MuonRelIso(mu)<<endl;// CombRelIso
  }
}

double
WZAnalyzer::CalcLeadPt(int type){
  if(type){
    double leadpt = -999.;
    if(wCand_ && wCand_.flavor() == type) 
      leadpt = TMath::Max(leadpt, wCand_.daughter(0)->pt());
    if(zCand_ && zCand_.flavor() == type){
      leadpt = TMath::Max(leadpt, zCand_.daughter(0)->pt());
      leadpt = TMath::Max(leadpt, zCand_.daughter(1)->pt());
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

/////////////////Modifies///////////////////////
void WZAnalyzer::CheckStream(ofstream& stream, string s){
  if(!stream) { 
    cout << "Cannot open file " << s << endl; 
    abort();
  } 
}

void WZAnalyzer::SetCandEvtFile(string s){
  outCandEvt.open(s.c_str());
  CheckStream(outCandEvt, s);
}

void WZAnalyzer::SetLogFile(string s){
  outLogFile.open(s.c_str());      
  CheckStream(outLogFile, s); 
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

bool WZAnalyzer::PassNoCut(){ 
  return true;
}

//Trigger requirements
//-----------------------------------------------------------
bool WZAnalyzer::PassTriggersCut()
{
//-----------------------------------------------------------
  if(debugme) cout<<"Trigger requirements"<<endl;
  
  const pat::TriggerPathRefVector acceptedPaths = triggerEvent_.acceptedPaths();
  if(debugme) cout<<"Using "<<acceptedPaths.size()<<" accepted paths from HLT"<<endl;
  for (size_t i = 0; i < acceptedPaths.size(); i++){
    string A = acceptedPaths[i]->name();
    for (size_t j = 0; j < triggersToUse_.size(); j++){
      //Need to acct for version numbers
//      cout<<"A: "<<acceptedPaths[i]->name()<< " B: "<<triggersToUse_[j]<<endl;
//      if(SameTrigger(A,B)){
      if(SameTrigger(A, triggersToUse_[j])){
//      if(SameTrigger(acceptedPaths[i]->name(), triggersToUse_[j])){
        if(debugme) cout<<"Match A: "<<acceptedPaths[i]->name()<<" B: "<<triggersToUse_[j]<<endl;
        if(acceptedPaths[i]->prescale() == 1  && acceptedPaths[i]->wasAccept()) return true;
        break;
      }
    }
  }
  return false;
}//--- PassTriggersCut()

bool
WZAnalyzer::PassNLeptonsCut(){
  return (looseElectrons_.size() + looseMuons_.size()) > minNLeptons_;
}

bool
WZAnalyzer::PassValidWandZCut(){
  if(!zCand_) CalcZVariables();
  if(!wCand_) CalcWVariables();
  CalcEventVariables();
  
  return PassValidZCut() && PassValidWCut();
}

bool
WZAnalyzer::PassValidWCut(){
  if(!wCand_) CalcWVariables();
  return wCand_ && wCand_.mt()>0;
}

bool
WZAnalyzer::PassValidZCut(){
  if(!zCand_) CalcZVariables();
  return zCand_ && zCand_.mass()>0.;
}

bool
WZAnalyzer::PassValidWZCandCut(){
  CalcWZVariables();
  return wzCand_.mass("minPz")>0.;
}

bool
WZAnalyzer::PassNumberOfZsCut(){
  return numZs_ < maxNumZs_;
}

bool
WZAnalyzer::PassLeadingLeptonPtCut(){
  return LeadPt_ > minLeadPt_;
}

bool
WZAnalyzer::PassMETCut(){
  return met_.et() > minMET_;
}

////////////////////////////////
/////////Check Z Properties/////
////////////////////////////////
bool
WZAnalyzer::PassZMassCut(){
  return (zCand_.mass() > minZmass_) && (zCand_.mass() < maxZmass_);  
}

bool
WZAnalyzer::PassZptCut(){
  return zCand_.pt() > minZpt_;
}
////////////////////////////////
/////////Check W Properties/////
////////////////////////////////

//Check W Transverse Mass
//-------------------s----------------------------------------
bool WZAnalyzer::PassWtransMassCut(){
  return wCand_.mt() > minWtransMass_;
}//--- PassWtransMassCut

bool
WZAnalyzer::PassWptCut(){
  return wCand_.pt() > minWpt_;
}

////////////////////////////////
/////////Check Electron Properties/////
////////////////////////////////
bool WZAnalyzer::PassElecLooseCut(const heep::Ele* elec){
  for(uint i=0; i<LooseElecCutFns_.size(); ++i){
    if(!(this->*LooseElecCutFns_[i])(elec)) return false;
  }
  return true;
}

bool WZAnalyzer::PassElecTightCut(const heep::Ele* elec){
  for(uint i=0; i<TightElecCutFns_.size(); ++i){
    if(!(this->*TightElecCutFns_[i])(elec)) return false;
  }
  return true;
}

bool WZAnalyzer::PassElecLooseEtCut(const heep::Ele* elec){
  if(debugme) cout<<"Check Electron Loose Et Cut"<<endl;
  return (elec->et() > minElecLooseEt_);
}//--- PassElecLooseEtCut

bool WZAnalyzer::PassElecTightEtCut(const heep::Ele* elec){
  if(debugme) cout<<"Check Electron Tight Et Cut"<<endl;
  return (elec->et() > minElecTightEt_);
}//--- PassElecTightEtCut

bool WZAnalyzer::PassElecLooseWPCut(const heep::Ele* elec){
//-----------------------------------------------------------
  if(debugme) cout<<"Check Electron WP Loose Cut"<<endl;
  return ((int)elec->patEle().electronID(cutElecWPLooseType_) & cutElecWPLooseMask_) 
    == cutElecWPLooseMask_;
}//--- PassElecLooseWPCut

bool WZAnalyzer::PassElecWPRelIsoCut(const heep::Ele* elec){
//-----------------------------------------------------------
  if(debugme) cout<<"Check Electron WP RelIso Cut"<<endl;
  return ((int)elec->patEle().electronID(cutElecWPTightType_) & cutWenuWPRelIsoMask_) 
    == cutWenuWPRelIsoMask_;
}//--- PassElecWPRelIsoElecCut

////////////////////////////////
/////////Check Muon Properties/////
////////////////////////////////
bool WZAnalyzer::PassMuonLooseCut(const TeVMuon* mu){
  for(uint i=0; i<LooseMuonCutFns_.size(); ++i){
    if(!(this->*LooseMuonCutFns_[i])(mu)) return false;
  }
  return true;
}

bool WZAnalyzer::PassMuonTightCut(const TeVMuon* mu){
  for(uint i=0; i<TightMuonCutFns_.size(); ++i){
    if(!(this->*TightMuonCutFns_[i])(mu)) return false;
  }
  return true;
}

bool WZAnalyzer::PassMuonLoosePtCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Loose Pt Cut"<<endl;
  return (mu->pt() > minMuonLoosePt_);
}

bool WZAnalyzer::PassMuonTightPtCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Tight Pt Cut"<<endl;
  return (mu->pt() > minMuonTightPt_);
}//--- PassMuonPtCut

bool WZAnalyzer::PassMuonGlobalCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Global Cut"<<endl;
  return (mu->isGlobalMuon()); 
}//--- PassMuonGlobalCut

bool WZAnalyzer::PassMuonNpixhitCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon NpixhitCut"<<endl;
  return (mu->globalTrack()->hitPattern().numberOfValidPixelHits() > minMuonNPixHit_);
}//--- PassMuonNpixhitCut

bool WZAnalyzer::PassMuonNtrkhitCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon NtrkhitCut"<<endl;
  return (mu->globalTrack()->hitPattern().numberOfValidTrackerHits() > minMuonNTrkHit_);
}//--- PassMuonNtrkhitCut

bool WZAnalyzer::PassMuonNormChi2Cut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Chi2 Cut"<<endl;
  return (mu->globalTrack()->normalizedChi2() < maxMuonNormChi2_);
}//--- PassMuonChi2Cut

bool WZAnalyzer::PassMuonHitsUsedCut(const TeVMuon* mu){
  //Num Valid Muon Hits
  if(debugme) cout<<"Check Muon Hits Used Cut"<<endl;
  return (mu->globalTrack()->hitPattern().numberOfValidMuonHits() > minMuonHitsUsed_);
}//--- PassMuonHits Used Cut

bool WZAnalyzer::PassMuonStationsCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Stations Cut"<<endl;
  return (mu->numberOfMatches() > minMuonStations_);
}//--- PassMuonStationsCut

bool WZAnalyzer::PassMuonEtaCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Eta Cut"<<endl;
  return (fabs(mu->eta()) < maxMuonEta_);
}//--- PassMuonEta Cut

bool WZAnalyzer::PassMuonCombRelIsoCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon CombRelIso Cut"<<endl;
  return (Calc_MuonRelIso(mu) < maxWmunuCombRelIso_);
}//--- PassMuonCombRelIsoCut

bool WZAnalyzer::PassMuonDxyCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Dxy Cut"<<endl;
  return (fabs(mu->userFloat("d0")) < maxMuonDxy_);
}//--- PassMuonDxyCut

////////////////////////////////
/////Check TeV Properties/////
////////////////////////////////

//Check Ht Properties
//-----------------------------------------------------------
bool WZAnalyzer::PassHtCut(){
//-----------------------------------------------------------
  if(debugme) cout<<"Check Ht Cuts"<<endl;
  return Ht_ > minHt_;   
}//--- PassHtCut

//Calc Ht (Cory: Does this take TeV values?)
//-----------------------------------------------------------
float WZAnalyzer::Calc_Ht(){
  return wCand_.daughter(0)->pt() +
    zCand_.daughter(0)->pt() +
    zCand_.daughter(1)->pt();
}//--- CalcHt

float WZAnalyzer::Calc_Q(){
  return wzCand_.mass("minPz") - zCand_.mass() - WMASS;
}

int WZAnalyzer::Calc_EvtType(){
  return (zCand_ && wCand_) ?  2 * (zCand_.flavor() != 11) + (wCand_.flavor() != 11) : -999;
}

float
WZAnalyzer::CalcElecSc(const heep::Ele* elec){
  float scEt=-999.;
  if (elec->patEle().core().isNonnull()) {
    const reco::SuperClusterRef sc = elec->patEle().superCluster();
    scEt = sc.get()->energy() / cosh(sc.get()->eta());
  }
  return scEt;
}

float
WZAnalyzer::Calc_MuonRelIso(const TeVMuon* mu){
  return (mu->isolationR03().emEt + mu->isolationR03().hadEt + mu->isolationR03().sumPt)
    / mu->pt();
}

///////////////Utilities//////////////////
//--------------------------------------------------------------
void WZAnalyzer::getEff(float & eff, float & deff, float Num, float Denom)
{
//--------------------------------------------------------------
  if(Denom){
    eff = Num/Denom;
    deff = TMath::Sqrt(eff * (1-eff)/Denom);
  }else{
    eff = deff = 0;
  }

}//getEff

void
WZAnalyzer::ClearEvtVariables(){
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
}

void WZAnalyzer::ClearAndResize(vector<TH1F*>& h, int& size, TH1F* ptr){
  h.clear();
  h.resize(size, ptr);
}

void WZAnalyzer::reportProgress(int eventNum) {
  if (eventNum % intOptions_["report"] == 0) {
    printf("\rWe've processed %i events so far...", eventNum);
    cout.flush();
    verbose("\n");
  }
  verbose("Event number: %i", ++eventNum);
}

/// Print to screen (like printf), but only if --verbose option is on
void WZAnalyzer::verbose(const char *string, ...)
{
  if (intOptions_["verbose"]) {
    va_list ap;
    va_start(ap, string);
    vprintf(string, ap);
    va_end(ap);
    cout << endl;
  }
}

//Cory: This should really be beginSample!!!
void WZAnalyzer::beginFile(std::vector<wprime::InputFile>::const_iterator fi){
  TFileDirectory dir = wprimeUtil_->getFileService()->mkdir(fi->samplename); 
  Declare_Histos(dir);
  ResetCounters();
}

// operations to be done when closing input file 
// (e.g. print summary)
void WZAnalyzer::endFile(std::vector<wprime::InputFile>::const_iterator fi,
                         ofstream & out){
  //ScaleHistos();//Already scaled
  printSummary(fi->samplename);  
  //deleteHistos();
  listOfHists.clear();
}

void WZAnalyzer::endAnalysis(ofstream & out){
}
