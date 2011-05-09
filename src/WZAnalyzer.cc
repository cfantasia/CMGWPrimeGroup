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

  intOptions_["report"] = cfg.getParameter<uint>("reportAfter");
  intOptions_["verbose"] = cfg.getParameter<bool>("debugme");
  doPreselect_ = cfg.getParameter<bool>("preselect");
  intOptions_["events"] = 0;

  debugme = cfg.getParameter<bool>("debugme");

  electronsLabel_ = cfg.getParameter<string>("electrons");
  muonsLabel_ = cfg.getParameter<string>("muons");
  metLabel_ = cfg.getParameter<string>("met");
  
  hltEventLabel_ = cfg.getParameter<string>("hltEventTag");
  pileupLabel_ = cfg.getParameter<string>("pileupTag");

  triggersToUse_          = cfg.getParameter<vstring>("triggersToUse");

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

  PDGZMASS = 91.1876; //GeV
  W_mass = 80.398;

  PI    = 2.0 * TMath::ACos(0.);
  TWOPI = 2.0 * PI;
  NOCUT = 9e9;

  hEvtType.resize(NCuts_,NULL);

  hWZInvMass     .resize(NCuts_,NULL);
  hWZ3e0muInvMass.resize(NCuts_,NULL);
  hWZ2e1muInvMass.resize(NCuts_,NULL);
  hWZ1e2muInvMass.resize(NCuts_,NULL);
  hWZ0e3muInvMass.resize(NCuts_,NULL);

  hWZTransMass.resize(NCuts_,NULL);
  hHt.resize(NCuts_,NULL);
  hWpt.resize(NCuts_,NULL);
  hZpt.resize(NCuts_,NULL);
  hMET.resize(NCuts_,NULL);

  hZMass     .resize(NCuts_,NULL);
  hZeeMass   .resize(NCuts_,NULL);
  hZmumuMass .resize(NCuts_,NULL);
  hZ3e0muMass.resize(NCuts_,NULL);
  hZ2e1muMass.resize(NCuts_,NULL);
  hZ1e2muMass.resize(NCuts_,NULL);
  hZ0e3muMass.resize(NCuts_,NULL);
  hZeeMassTT.resize(NCuts_,NULL);
  hZeeMassTF.resize(NCuts_,NULL);
  hZmumuMassTT.resize(NCuts_,NULL);
  hZmumuMassTF.resize(NCuts_,NULL);

  hWTransMass     .resize(NCuts_,NULL);
  hWenuTransMass  .resize(NCuts_,NULL);
  hWmunuTransMass .resize(NCuts_,NULL);
  hW3e0muTransMass.resize(NCuts_,NULL);
  hW2e1muTransMass.resize(NCuts_,NULL);
  hW1e2muTransMass.resize(NCuts_,NULL);
  hW0e3muTransMass.resize(NCuts_,NULL);

  hQ.resize(NCuts_,NULL);

  hLeadPt.resize(NCuts_,NULL);
  hLeadElecPt.resize(NCuts_,NULL);
  hLeadMuonPt.resize(NCuts_,NULL);

  hElecPt.resize(NCuts_,NULL);
  hElecEt.resize(NCuts_,NULL);
  hElecdEta.resize(NCuts_,NULL);
  hElecdPhi.resize(NCuts_,NULL);
  hElecSigmann.resize(NCuts_,NULL);
  hElecEP.resize(NCuts_,NULL);
  hElecHE.resize(NCuts_,NULL);
  hElecTrkRelIso.resize(NCuts_,NULL);
  hElecECalRelIso.resize(NCuts_,NULL);
  hElecHCalRelIso.resize(NCuts_,NULL);

  hMuonPt.resize(NCuts_,NULL);
  hMuonDxy.resize(NCuts_,NULL);
  hMuonNormChi2.resize(NCuts_,NULL);
  hMuonNPix.resize(NCuts_,NULL);
  hMuonNTrk.resize(NCuts_,NULL);
  hMuonRelIso.resize(NCuts_,NULL);
  hMuonStation.resize(NCuts_,NULL);
  hMuonSip.resize(NCuts_,NULL);

// +++++++++++++++++++General Cut values
  maxNumZs = cfg.getParameter<int>("maxNumZs");
  minNumLeptons = cfg.getParameter<int>("minNumLeptons");
  minMET = cfg.getParameter<double>("minMET");

// +++++++++++++++++++Ht Cuts
  minHt = cfg.getParameter<double>("minHt");

// +++++++++++++++++++W Cuts
  minWtransMass = cfg.getParameter<double>("minWtransMass");
  minWpt = cfg.getParameter<double>("minWpt");

  maxWmunuCombRelIso = cfg.getParameter<double>("maxWmunuCombRelIso");
  cutWenuWPRelIsoMask = cfg.getParameter<int>("cutWenuWPRelIsoMask");

// +++++++++++++++++++Z Cuts
  minZpt = cfg.getParameter<double>("minZpt");
  minZmass = cfg.getParameter<double>("minZmass");
  maxZmass = cfg.getParameter<double>("maxZmass");

// +++++++++++++++++++Electron General Cuts
//VBTF Recommended Cuts
  minElecLooseEt = cfg.getParameter<double>("minElecLooseEt");
  minElecTightEt = cfg.getParameter<double>("minElecTightEt");
  cutElecWPLooseMask = cfg.getParameter<int>("cutElecWPLooseMask");

  maxElecSigmaiEtaiEta = cfg.getParameter<vector<double> >("maxElecSigmaiEtaiEta");
  maxElecDeltaPhiIn  = cfg.getParameter<vector<double> >("maxElecDeltaPhiIn");
  maxElecDeltaEtaIn  = cfg.getParameter<vector<double> >("maxElecDeltaEtaIn");
  maxElecHOverE      = cfg.getParameter<vector<double> >("maxElecHOverE");

// +++++++++++++++++++Muon General Cuts
  maxMuonEta = cfg.getParameter<double>("maxMuonEta");
  minMuonLoosePt = cfg.getParameter<double>("minMuonLoosePt");
  minMuonTightPt = cfg.getParameter<double>("minMuonTightPt");
//VBTF Recommended Cuts
  maxMuonDxy = cfg.getParameter<double>("maxMuonDxy");
  maxMuonNormChi2 = cfg.getParameter<double>("maxMuonNormChi2");
  minMuonNPixHit = cfg.getParameter<int>("minMuonNPixHit");
  minMuonNTrkHit = cfg.getParameter<int>("minMuonNTrkHit");
  minMuonStations = cfg.getParameter<int>("minMuonStations");
  minMuonHitsUsed = cfg.getParameter<int>("minMuonHitsUsed");
}
WZAnalyzer::~WZAnalyzer(){
  outCandEvt.close(); 
  outLogFile.close(); 
}

void WZAnalyzer::FillCutFns(){
  mFnPtrs_["NoCuts"] = &WZAnalyzer::PassNoCut;
  mFnPtrs_["HLT"] = &WZAnalyzer::PassTriggersCut;
  mFnPtrs_["ValidW"] = &WZAnalyzer::PassValidWCut;
  mFnPtrs_["ValidZ"] = &WZAnalyzer::PassValidZCut;
  mFnPtrs_["ValidWandZ"] = &WZAnalyzer::PassValidWandZCut;
  mFnPtrs_["ValidWZCand"] = &WZAnalyzer::PassValidWZCandCut;
  mFnPtrs_["NumZs"] = &WZAnalyzer::PassNumberOfZsCut;
//  mFnPtrs_["LooseElec"] = &WZAnalyzer::PassWZElecLooseCut;
//  mFnPtrs_["LooseMuon"] = &WZAnalyzer::PassWZMuonLooseCut;
//  mFnPtrs_["LooseZElec"] = &WZAnalyzer::PassZElecLooseCut;
//  mFnPtrs_["LooseZMuon"] = &WZAnalyzer::PassZMuonLooseCut;
//  mFnPtrs_["WLepIso"] = &WZAnalyzer::PassWLepIsoCut;
  mFnPtrs_["ZMass"] = &WZAnalyzer::PassZMassCut;
  mFnPtrs_["WTransMass"] = &WZAnalyzer::PassWtransMassCut;
  mFnPtrs_["MET"] = &WZAnalyzer::PassMETCut;
//  mFnPtrs_["ZLepPt"] = &WZAnalyzer::PassZLepPtCut;
//  mFnPtrs_["WLepPt"] = &WZAnalyzer::PassWLepPtCut;
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
  Num_surv_cut_.resize(NCuts_,0.);
  eventNum = 0;
  runNumber = -1;
  lumiID = -1;
}

//--------------------------------------------------------------
void WZAnalyzer::Declare_Histos(TFileDirectory & dir)
{
  verbose("Declare histos\n");

  DeclareHistoSet("hEvtType", "Event Type (Number of Electrons)",
                  "N_{Elec}", 4, 0, 4, hEvtType,dir);

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
  hNumEvts = dir.make<TH1F>("hNumEvts",title.c_str(),NCuts_,0,NCuts_);
  hEffRel = dir.make<TH1F>("hEffRel","Relative Efficiency",NCuts_,0,NCuts_);
  hEffAbs = dir.make<TH1F>("hEffAbs","Absolute Efficiency",NCuts_,0,NCuts_);

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
  if(wCand && zCand){
    hWZInvMass[index]->Fill(wzCand.mass("minPz"));
    if     (evtType == 0) hWZ3e0muInvMass[index]->Fill(wzCand.mass("minPz"));
    else if(evtType == 1) hWZ2e1muInvMass[index]->Fill(wzCand.mass("minPz"));
    else if(evtType == 2) hWZ1e2muInvMass[index]->Fill(wzCand.mass("minPz"));
    else if(evtType == 3) hWZ0e3muInvMass[index]->Fill(wzCand.mass("minPz"));
    hEvtType[index]->Fill(evtType);
    hQ[index]->Fill(Q); 
    hWZTransMass[index]->Fill(wzCand.transMass());
    hHt[index]->Fill(Ht);
  }
  if(zCand){
    hZpt[index]->Fill(zCand.pt());
    hZMass[index]->Fill(zCand.mass());
    if      (zCand.flavor() == PDGELEC){
      hZeeMass[index]->Fill(zCand.mass());
    }else if (zCand.flavor() == PDGMUON){
      hZmumuMass[index]->Fill(zCand.mass());
    }
    if     (evtType == 0) hZ3e0muMass[index]->Fill(zCand.mass());
    else if(evtType == 1) hZ2e1muMass[index]->Fill(zCand.mass());
    else if(evtType == 2) hZ1e2muMass[index]->Fill(zCand.mass());
    else if(evtType == 3) hZ0e3muMass[index]->Fill(zCand.mass());
  }
  if(wCand){
    hWpt[index]->Fill(wCand.pt());
    hWTransMass[index]->Fill(wCand.mt());
    if      (wCand.flavor() == PDGELEC) hWenuTransMass[index]->Fill(wCand.mt());
    else if (wCand.flavor() == PDGMUON) hWmunuTransMass[index]->Fill(wCand.mt());
    if(evtType == 0) hW3e0muTransMass[index]->Fill(wCand.mt());
    if(evtType == 1) hW2e1muTransMass[index]->Fill(wCand.mt());
    if(evtType == 2) hW1e2muTransMass[index]->Fill(wCand.mt());
    if(evtType == 3) hW0e3muTransMass[index]->Fill(wCand.mt());
  }  
  hMET[index]->Fill(met.et());
}//Fill_Histos

void
WZAnalyzer::CalcEventVariables(){
  if (debugme) cout<<"In CalcEventVariables\n";
////Calculate Important Quantities for each event
  Ht = (zCand && wCand) ? Calc_Ht() : -999.;
  Q = (zCand && wCand) ? Calc_Q() : -999.;
  evtType = (zCand && wCand) ? Calc_EvtType() : -999;
  LeadPt = CalcLeadPt();
  LeadElecPt = CalcLeadPt(PDGELEC);
  LeadMuonPt = CalcLeadPt(PDGMUON);
//  TT = TF = false;
//  if(zCand.flavor()){
//    bool tight1 = PassTightCut(zCand.daughter(0), zCand.flavor());
//    bool tight2 = PassTightCut(zCand.daughter(1), zCand.flavor());
//    //cout<<"tight1: "<<tight1<<" tight2: "<<tight2<<endl;
//    TT = tight1 && tight2;
//    TF = (tight1 && !tight2) || (!tight1 && tight2);
//  }
}

void
WZAnalyzer::DeclareHistoSet(string n, string t, string xtitle,
                            int nbins, float min, float max,
                            vector<TH1F*>& h, TFileDirectory& d){
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
  ++Num_surv_cut_[cut_index];
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
    hNumEvts->Fill(i,Num_surv_cut_[i]);

//calculate efficiencies
    float eff, deff;
    if(i == 0){
      getEff(eff, deff, Num_surv_cut_[i], Num_surv_cut_[0]);
    }else{
      getEff(eff, deff, Num_surv_cut_[i], Num_surv_cut_[i-1]);
    }
    hEffRel->Fill(i,eff*100);
    outLogFile << setw(15) <<"\tRelative eff = "<<setw(6)<<eff*100 << " +/- " << setw(6)<<deff*100 << "%";
    if(Num_surv_cut_[i-1] && Num_surv_cut_[i-1] < 1.) printf("eff:%.2f deff:%.2f num:%.2f den:%.2f\n",eff,deff,Num_surv_cut_[i],Num_surv_cut_[i-1]);    
    getEff(eff, deff, Num_surv_cut_[i], Num_surv_cut_[0]);
    hEffAbs->Fill(i,eff*100);
    outLogFile << setw(15) <<"\tAbsolute eff = "<<setw(6)<<eff*100 << " +/- " << setw(6)<<deff*100 << "%"
               << endl;
        
  } // loop over different cuts
}//printSummary


void 
WZAnalyzer::eventLoop(edm::EventBase const & event){
  if (intOptions_["events"] && eventNum > intOptions_["events"]) return;
  reportProgress(eventNum++);
  //updateEventCounts(event, nEvents, runNumber,  wprimeUtil_->getLumi_ipb());
  datasetName = getDatasetName(event, datasetName);

  // Preselection - skip events that don't look promising
  if (doPreselect_)
    if (getProduct<double>(event, 
                           "wzPreselectionProducer:ZMassDiff") > 12.5 ||
        getProduct<double>(event, 
                           "wzPreselectionProducer:highestLeptonPt") < 10 ||
        getProduct<vector<uint> >(event, 
                                  "wzPreselectionProducer:nLeptonsEid")[5] < 3)
      return;
//  // Get information about flavorHistory
//  const uint flavorType = getProduct<uint>(event, "flavorHistoryFilter", 0);
//  verbose("    FlavorType: %i", flavorType);

  // Get leptons
  const ElectronV electrons = getProduct<ElectronV>(event,electronsLabel_);
  const MuonV muons = getProduct<MuonV>(event, muonsLabel_);
  verbose("    Contains: %i electron(s), %i muon(s)",
          electrons.size(), muons.size());

  // Get various types of MET
//  METV mets;
//  VInputTag metTags = pset_.getParameter<VInputTag>("mets");
//  for (size_t i = 0; i < metTags.size(); i++)
//    mets.push_back(getProduct<METV>(event, string(metTags[i].label()))[0]);
//  met = mets[2];
  met = getProduct<METV>(event, metLabel_)[0];

  // Make vectors of leptons passing various criteria
  ElectronV looseElectrons, tightElectrons;
  for (size_t i = 0; i < electrons.size(); i++) {
    if (PassElecLooseCut(&electrons[i])){
      looseElectrons.push_back(electrons[i]);
      if (PassElecTightCut(&electrons[i]))
        tightElectrons.push_back(electrons[i]);
    }
  }

  MuonV looseMuons, tightMuons;
  for (size_t i = 0; i < muons.size(); i++) {
    if (PassMuonLooseCut(&muons[i])){
      looseMuons.push_back(muons[i]);
      if (PassMuonTightCut(&muons[i]))
        tightMuons.push_back(muons[i]);
    }
  }
/*
  GenParticleV genParticles = getUntrackedProduct<GenParticleV>(event, "genParticles");
  const Candidate * genZ = 0;
  const Candidate * genW = 0;
  for (size_t i = 0; i < genParticles.size(); i++)
    if (abs(genParticles[i].pdgId()) == 23) genZ = & genParticles[i];
    else if (abs(genParticles[i].pdgId()) == 24) genW = & genParticles[i];
*/
  // Reconstruct the Z
  ZCandV looseZCands = getZCands(looseElectrons, looseMuons, 12.5);
  zCand = looseZCands.size() ? looseZCands[0] : ZCandidate();
  numZs = looseZCands.size();
  verbose("    Contains: %i loose Z candidate(s)", looseZCands.size());

  // Reconstruct the W
  wCand = getWCand(tightElectrons, tightMuons, met, zCand);

  wzCand = (zCand && wCand) ? WZCandidate(zCand, wCand) : WZCandidate();

  triggerEvent_ = getProduct<pat::TriggerEvent>(event,hltEventLabel_); 

///////////////////////
/*
  PupInfo_ = getProduct<std::vector< PileupSummaryInfo > >(event, pileupLabel_);   
  
  std::vector<PileupSummaryInfo>::const_iterator PVI;                       
  for(PVI = PupInfo_.begin(); PVI != PupInfo_.end(); ++PVI) {               
    //Cory: What am I looping over?
    int PU_BunchCrossing = PVI->getBunchCrossing();                  
    int PU_NumInteractions = PVI->getPU_NumInteractions();           
  }
*/
/////////////////////
  CalcEventVariables();

  if(!PassCuts()) return;
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
      <<" zCand.flavor(): "<<zCand.flavor()
      <<" zCand.mass(): "<<zCand.mass()
      <<" wCand.flavor(): "<<wCand.flavor()
      <<" wCand.mt(): "<<wCand.mt()
      <<endl;

  cout<<" Z lep1 pt "<<zCand.daughter(0)->pt()
      <<" Z lep2.pt "<<zCand.daughter(1)->pt()
      <<endl;

  cout<<" W lep pt "<<wCand.daughter(0)->pt()
      <<" pfMet_et: "<<met.et()
      <<endl;

  cout<<" Ht: "<<Ht
      <<" Zpt: "<<zCand.pt()
      <<" Wpt: "<<wCand.pt()<<endl
      <<" WZ Mass: "<<wzCand.mass("minPz")<<endl
//      <<" # Elec: "<<elec.pt->size()
//      <<" # Muon: "<<mu.pt->size()
      <<endl<<endl;
  return;
}

void
WZAnalyzer::PrintEventFull(edm::EventBase const & event){
  PrintEvent(event);
  if     (zCand.flavor() == PDGELEC){
    PrintElectron((const pat::Electron*)zCand.daughter(0), PDGZ);
    PrintElectron((const pat::Electron*)zCand.daughter(1), PDGZ);
  }else if(zCand.flavor() == PDGMUON){
    PrintMuon((const pat::Muon*)zCand.daughter(0), PDGZ);
    PrintMuon((const pat::Muon*)zCand.daughter(1), PDGZ);
  }

  if     (wCand.flavor() == PDGELEC) PrintElectron((const pat::Electron*)wCand.daughter(0), PDGW);
  else if(wCand.flavor() == PDGMUON) PrintMuon    ((const pat::Muon*    )wCand.daughter(0), PDGW);
}

void
WZAnalyzer::PrintElectron(const pat::Electron* elec, int parent){
  if     (parent == PDGZ) cout<<"-----Electron from Z-------------------------"<<endl;
  else if(parent == PDGW) cout<<"-----Electron from W-------------------------"<<endl;
  else                    cout<<"-----Electron from ?-------------------------"<<endl;
  cout<<" Elec Pt: "<<elec->pt()<<endl
      <<" Elec ScEt: "<<CalcElecSc(elec)<<endl //ScEt
      <<" Elec Eta: "<<elec->eta()<<endl //Eta
      <<" Elec SigmaNN: "<<elec->sigmaIetaIeta()<<endl //sigmaNN
      <<" Elec dPhi: "<<elec->deltaPhiSuperClusterTrackAtVtx()<<endl //DeltaPhi
      <<" Elec dEta: "<<elec->deltaEtaSuperClusterTrackAtVtx()<<endl //DeltaEta
      <<" Elec HoverE: "<<elec->hadronicOverEm()<<endl// H/E
      <<" Elec EoverP: "<<elec->eSuperClusterOverP()<<endl// E/P
      <<" Elec WP95: "<<elec->electronID("simpleEleId95relIso")<<endl
      <<" Elec WP90: "<<elec->electronID("simpleEleId90relIso")<<endl
      <<" Elec WP85: "<<elec->electronID("simpleEleId85relIso")<<endl
      <<" Elec WP80: "<<elec->electronID("simpleEleId80relIso")<<endl;
}

void
WZAnalyzer::PrintMuon(const pat::Muon* mu, int parent){
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
    if(wCand && wCand.flavor() == type) 
      leadpt = TMath::Max(leadpt, wCand.daughter(0)->pt());
    if(zCand && zCand.flavor() == type){
      leadpt = TMath::Max(leadpt, zCand.daughter(0)->pt());
      leadpt = TMath::Max(leadpt, zCand.daughter(1)->pt());
    }
    return leadpt;
  }
  return TMath::Max(CalcLeadPt(PDGELEC), CalcLeadPt(PDGMUON));
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
  for (size_t i = 0; i < acceptedPaths.size(); i++)
    for (size_t j = 0; j < triggersToUse_.size(); j++)
      if (!acceptedPaths[i]->name().compare(triggersToUse_[j])){
//        string pat::TriggerPath::name()
//        unsigned pat::TriggerPath::prescale()
//        bool pat::TriggerPath::wasAccept()
        if(acceptedPaths[i]->prescale() == 1  && acceptedPaths[i]->wasAccept()) return true;
        break;
      }
  return false;
}//--- PassTriggersCut()

bool
WZAnalyzer::PassValidWandZCut(){
  return PassValidZCut() && PassValidWCut();
}

bool
WZAnalyzer::PassValidWCut(){
  if(wCand) return true;
  else      return false;
}

bool
WZAnalyzer::PassValidZCut(){
  return zCand.mass()>0.;
}

bool
WZAnalyzer::PassValidWZCandCut(){
  return wzCand.mass("minPz")>0.;
}

bool
WZAnalyzer::PassNumberOfZsCut(){
  return numZs < maxNumZs;
}

bool
WZAnalyzer::PassMETCut(){
  return met.et() > minMET;
}

////////////////////////////////
/////////Check Z Properties/////
////////////////////////////////
bool
WZAnalyzer::PassZMassCut(){
  return (zCand.mass() > minZmass) && (zCand.mass() < maxZmass);  
}

bool
WZAnalyzer::PassZptCut(){
  return zCand.pt() > minZpt;
}
////////////////////////////////
/////////Check W Properties/////
////////////////////////////////

//Check W Transverse Mass
//-------------------s----------------------------------------
bool WZAnalyzer::PassWtransMassCut(){
  return wCand.mt() > minWtransMass;
}//--- PassWtransMassCut

bool
WZAnalyzer::PassWptCut(){
  return wCand.pt() > minWpt;
}

////////////////////////////////
/////////Check Electron Properties/////
////////////////////////////////
bool WZAnalyzer::PassElecLooseCut(const pat::Electron* elec){
  for(uint i=0; i<LooseElecCutFns_.size(); ++i){
    if(!(this->*LooseElecCutFns_[i])(elec)) return false;
  }
  return true;
}

bool WZAnalyzer::PassElecTightCut(const pat::Electron* elec){
  for(uint i=0; i<LooseElecCutFns_.size(); ++i){
    if(!(this->*LooseElecCutFns_[i])(elec)) return false;
  }
  return true;
}

bool WZAnalyzer::PassElecLooseEtCut(const pat::Electron* elec){
  if(debugme) cout<<"Check Electron Loose Et Cut"<<endl;
  return (CalcElecSc(elec) > minElecLooseEt);
}//--- PassElecLooseEtCut

bool WZAnalyzer::PassElecTightEtCut(const pat::Electron* elec){
  if(debugme) cout<<"Check Electron Tight Et Cut"<<endl;
  return (CalcElecSc(elec) > minElecTightEt);
}//--- PassElecTightEtCut

bool WZAnalyzer::PassElecLooseWPCut(const pat::Electron* elec){
//-----------------------------------------------------------
  if(debugme) cout<<"Check Electron WP Loose Cut"<<endl;
  return ((int)elec->electronID("simpleEleId95relIso") & cutElecWPLooseMask) 
    == cutElecWPLooseMask;
}//--- PassElecLooseWPCut

bool WZAnalyzer::PassElecWPRelIsoCut(const pat::Electron* elec){
//-----------------------------------------------------------
  if(debugme) cout<<"Check Electron WP RelIso Cut"<<endl;
  return ((int)elec->electronID("simpleEleId80relIso") & cutWenuWPRelIsoMask) 
    == cutWenuWPRelIsoMask;
}//--- PassElecWPRelIsoElecCut

////////////////////////////////
/////////Check Muon Properties/////
////////////////////////////////
bool WZAnalyzer::PassMuonLooseCut(const pat::Muon* mu){
  for(uint i=0; i<LooseMuonCutFns_.size(); ++i){
    if(!(this->*LooseMuonCutFns_[i])(mu)) return false;
  }
  return true;
}

bool WZAnalyzer::PassMuonTightCut(const pat::Muon* mu){
  for(uint i=0; i<TightMuonCutFns_.size(); ++i){
    if(!(this->*TightMuonCutFns_[i])(mu)) return false;
  }
  return true;
}

bool WZAnalyzer::PassMuonLoosePtCut(const pat::Muon* mu){
  if(debugme) cout<<"Check Muon Loose Pt Cut"<<endl;
  return (mu->pt() > minMuonLoosePt);
}

bool WZAnalyzer::PassMuonTightPtCut(const pat::Muon* mu){
  if(debugme) cout<<"Check Muon Tight Pt Cut"<<endl;
  return (mu->pt() > minMuonTightPt);
}//--- PassMuonPtCut

bool WZAnalyzer::PassMuonGlobalCut(const pat::Muon* mu){
  if(debugme) cout<<"Check Muon Global Cut"<<endl;
  return (mu->isGlobalMuon()); 
}//--- PassMuonGlobalCut

bool WZAnalyzer::PassMuonNpixhitCut(const pat::Muon* mu){
  if(debugme) cout<<"Check Muon NpixhitCut"<<endl;
  return (mu->globalTrack()->hitPattern().numberOfValidPixelHits() > minMuonNPixHit);
}//--- PassMuonNpixhitCut

bool WZAnalyzer::PassMuonNtrkhitCut(const pat::Muon* mu){
  if(debugme) cout<<"Check Muon NtrkhitCut"<<endl;
  return (mu->globalTrack()->hitPattern().numberOfValidTrackerHits() > minMuonNTrkHit);
}//--- PassMuonNtrkhitCut

bool WZAnalyzer::PassMuonNormChi2Cut(const pat::Muon* mu){
  if(debugme) cout<<"Check Muon Chi2 Cut"<<endl;
  return (mu->globalTrack()->normalizedChi2() < maxMuonNormChi2);
}//--- PassMuonChi2Cut

bool WZAnalyzer::PassMuonHitsUsedCut(const pat::Muon* mu){
  //Num Valid Muon Hits
  if(debugme) cout<<"Check Muon Hits Used Cut"<<endl;
  return (mu->globalTrack()->hitPattern().numberOfValidMuonHits() > minMuonHitsUsed);
}//--- PassMuonHits Used Cut

bool WZAnalyzer::PassMuonStationsCut(const pat::Muon* mu){
  if(debugme) cout<<"Check Muon Stations Cut"<<endl;
  return (mu->numberOfMatches() > minMuonStations);
}//--- PassMuonStationsCut

bool WZAnalyzer::PassMuonEtaCut(const pat::Muon* mu){
  if(debugme) cout<<"Check Muon Eta Cut"<<endl;
  return (fabs(mu->eta()) < maxMuonEta);
}//--- PassMuonEta Cut

bool WZAnalyzer::PassMuonCombRelIsoCut(const pat::Muon* mu){
  if(debugme) cout<<"Check Muon CombRelIso Cut"<<endl;
  return (Calc_MuonRelIso(mu) < maxWmunuCombRelIso);
}//--- PassMuonCombRelIsoCut

bool WZAnalyzer::PassMuonDxyCut(const pat::Muon* mu){
  if(debugme) cout<<"Check Muon Dxy Cut"<<endl;
  return (fabs(mu->userFloat("d0")) < maxMuonDxy);
}//--- PassMuonDxyCut

////////////////////////////////
/////Check Other Properties/////
////////////////////////////////

//Check Ht Properties
//-----------------------------------------------------------
bool WZAnalyzer::PassHtCut(){
//-----------------------------------------------------------
  if(debugme) cout<<"Check Ht Cuts"<<endl;
   
  return Ht > minHt;
    
}//--- PassHtCut

//Calc Ht
//-----------------------------------------------------------
float WZAnalyzer::Calc_Ht(){
  return wCand.daughter(0)->pt() +
    zCand.daughter(0)->pt() +
    zCand.daughter(1)->pt();
}//--- CalcHt

float WZAnalyzer::Calc_Q(){
  return wzCand.mass("minPz") - zCand.mass() - W_mass;
}

int WZAnalyzer::Calc_EvtType(){
  return (zCand && wCand) ?  2 * (zCand.flavor() != 11) + (wCand.flavor() != 11) : -999;

}

float
WZAnalyzer::CalcElecSc(const pat::Electron* elec){
  float scEt=-999.;
  if (elec->core().isNonnull()) {
    const reco::SuperClusterRef sc = elec->superCluster();
    scEt = sc.get()->energy() / cosh(sc.get()->eta());
  }
  return scEt;
}

float
WZAnalyzer::Calc_MuonRelIso(const pat::Muon* mu){
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

//--------------------------------------------------------------
double WZAnalyzer::deltaPhi(double phi1, double phi2)
{
//--------------------------------------------------------------

  double phi=fabs(phi1-phi2);
  return (phi <= PI) ? phi : TWOPI - phi;
}

//--------------------------------------------------------------
double WZAnalyzer::deltaEta(double eta1, double eta2)
{
//--------------------------------------------------------------
  double eta = fabs(eta1 - eta2);
  return eta;
}

//Just a function to calculate DeltaR
//--------------------------------------------------------------
double WZAnalyzer::deltaR(double eta1, double phi1, double eta2, double phi2)
{
//--------------------------------------------------------------
  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return sqrt(deta * deta + dphi * dphi);
}

////////////
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
  ScaleHistos();
  printSummary(fi->samplename);  
  listOfHists.clear();
}

void WZAnalyzer::endAnalysis(ofstream & out){
}
