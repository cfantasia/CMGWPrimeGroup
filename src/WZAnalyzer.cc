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

  muonAlgo_ = cfg.getParameter<int>("muonAlgo");
  
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
  if(wCand && zCand){
    hWZInvMass[index]->Fill(wzCand.mass("minPz"), weight);
    if     (evtType == 0) hWZ3e0muInvMass[index]->Fill(wzCand.mass("minPz"), weight);
    else if(evtType == 1) hWZ2e1muInvMass[index]->Fill(wzCand.mass("minPz"), weight);
    else if(evtType == 2) hWZ1e2muInvMass[index]->Fill(wzCand.mass("minPz"), weight);
    else if(evtType == 3) hWZ0e3muInvMass[index]->Fill(wzCand.mass("minPz"), weight);
    hEvtType[index]->Fill(evtType, weight);
    hQ[index]->Fill(Q, weight); 
    hWZTransMass[index]->Fill(wzCand.transMass(), weight);
    hHt[index]->Fill(Ht, weight);
  }
  if(zCand){
    hZpt[index]->Fill(zCand.pt(), weight);
    hZMass[index]->Fill(zCand.mass(), weight);
    if      (zCand.flavor() == PDGELEC){
      hZeeMass[index]->Fill(zCand.mass(), weight);
    }else if (zCand.flavor() == PDGMUON){
      hZmumuMass[index]->Fill(zCand.mass(), weight);
    }
    if     (evtType == 0) hZ3e0muMass[index]->Fill(zCand.mass(), weight);
    else if(evtType == 1) hZ2e1muMass[index]->Fill(zCand.mass(), weight);
    else if(evtType == 2) hZ1e2muMass[index]->Fill(zCand.mass(), weight);
    else if(evtType == 3) hZ0e3muMass[index]->Fill(zCand.mass(), weight);
  }
  if(wCand){
    hWpt[index]->Fill(wCand.pt(), weight);
    hWTransMass[index]->Fill(wCand.mt(), weight);
    if      (wCand.flavor() == PDGELEC) hWenuTransMass[index]->Fill(wCand.mt(), weight);
    else if (wCand.flavor() == PDGMUON) hWmunuTransMass[index]->Fill(wCand.mt(), weight);
    if(evtType == 0) hW3e0muTransMass[index]->Fill(wCand.mt(), weight);
    if(evtType == 1) hW2e1muTransMass[index]->Fill(wCand.mt(), weight);
    if(evtType == 2) hW1e2muTransMass[index]->Fill(wCand.mt(), weight);
    if(evtType == 3) hW0e3muTransMass[index]->Fill(wCand.mt(), weight);
  }  
  hMET[index]->Fill(met.et(), weight);
}//Fill_Histos

void
WZAnalyzer::CalcEventVariables(){
  if (debugme) cout<<"In CalcEventVariables\n";
////Calculate Important Quantities for each event
  Ht = (zCand && wCand) ? Calc_Ht() : -999.;
  Q = (zCand && wCand) ? Calc_Q() : -999.;
  evtType = (zCand && wCand) ? Calc_EvtType() : -999;
  verbose("evt Type: %i, Z Flav: %i, W Flav: %i\n", evtType, zCand.flavor(), wCand.flavor());

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
  const ElectronV patElectrons = getProduct<ElectronV>(event,electronsLabel_);
  const vector<pat::Muon> patMuons = getProduct<vector<pat::Muon> >(event, muonsLabel_);
  verbose("    Contains: %i electron(s), %i muon(s)",
          patElectrons.size(), patMuons.size());

  // Get various types of MET
//  METV mets;
//  VInputTag metTags = pset_.getParameter<VInputTag>("mets");
//  for (size_t i = 0; i < metTags.size(); i++)
//    mets.push_back(getProduct<METV>(event, string(metTags[i].label()))[0]);
//  met = mets[2];
  met = getProduct<METV>(event, metLabel_)[0];

  // Make vectors of leptons passing various criteria
  ElectronV looseElectrons, tightElectrons;
  for (size_t i = 0; i < patElectrons.size(); i++) {
    if (PassElecLooseCut(&patElectrons[i])){
      looseElectrons.push_back(patElectrons[i]);
      if (PassElecTightCut(&patElectrons[i]))
        tightElectrons.push_back(patElectrons[i]);
    }
  }

  MuonV muons, looseMuons, tightMuons;
  for (size_t i = 0; i < patMuons.size(); i++) {
    muons.push_back(TeVMuon(patMuons[i],muonAlgo_));   
    if (PassMuonLooseCut(&muons[i])){
      looseMuons.push_back(muons[i]);
      if (PassMuonTightCut(&muons[i]))
        tightMuons.push_back(muons[i]);
    }
  }

  verbose("    Contains: %i loose electron(s), %i loose muon(s)",
          looseElectrons.size(), looseMuons.size());

  verbose("    Contains: %i tight electron(s), %i tightmuon(s)",
          tightElectrons.size(), tightMuons.size());

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
  wCand = getWCand(tightElectrons, tightMuons, met, zCand, minDeltaR_);
  verbose("    Contains: %i tight W candidate(s)", (bool)wCand);

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

      <<" zCand.flavor(): "<<zCand.flavor()
      <<" zCand.mass(): "<<zCand.mass()
      <<" wCand.flavor(): "<<wCand.flavor()
      <<" wCand.mt(): "<<wCand.mt()
      <<endl;

  cout<<" Z lep1 pt "<<zCand.daughter(0)->pt()
      <<" Z lep2 pt "<<zCand.daughter(1)->pt()
      <<endl;

  cout<<" W lep pt "<<wCand.daughter(0)->pt()
      <<" pfMet et: "<<met.et()
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
    PrintMuon((const TeVMuon*)zCand.daughter(0), PDGZ);
    PrintMuon((const TeVMuon*)zCand.daughter(1), PDGZ);
  }

  if     (wCand.flavor() == PDGELEC) PrintElectron((const pat::Electron*)wCand.daughter(0), PDGW);
  else if(wCand.flavor() == PDGMUON) PrintMuon    ((const TeVMuon*    )wCand.daughter(0), PDGW);
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
        if(debugme) cout<<"Match A: "<<acceptedPaths[i]->name()<< " B: "<<triggersToUse_[j]<<endl;
        if(acceptedPaths[i]->prescale() == 1  && acceptedPaths[i]->wasAccept()) return true;
        break;
      }
    }
  }
  return false;
}//--- PassTriggersCut()

bool
WZAnalyzer::PassValidWandZCut(){
  return PassValidZCut() && PassValidWCut();
}

bool
WZAnalyzer::PassValidWCut(){
  return wCand && wCand.mt()>0;
}

bool
WZAnalyzer::PassValidZCut(){
  return zCand && zCand.mass()>0.;
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
  for(uint i=0; i<TightElecCutFns_.size(); ++i){
    if(!(this->*TightElecCutFns_[i])(elec)) return false;
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
  return (mu->pt() > minMuonLoosePt);
}

bool WZAnalyzer::PassMuonTightPtCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Tight Pt Cut"<<endl;
  return (mu->pt() > minMuonTightPt);
}//--- PassMuonPtCut

bool WZAnalyzer::PassMuonGlobalCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Global Cut"<<endl;
  return (mu->isGlobalMuon()); 
}//--- PassMuonGlobalCut

bool WZAnalyzer::PassMuonNpixhitCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon NpixhitCut"<<endl;
  return (mu->globalTrack()->hitPattern().numberOfValidPixelHits() > minMuonNPixHit);
}//--- PassMuonNpixhitCut

bool WZAnalyzer::PassMuonNtrkhitCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon NtrkhitCut"<<endl;
  return (mu->globalTrack()->hitPattern().numberOfValidTrackerHits() > minMuonNTrkHit);
}//--- PassMuonNtrkhitCut

bool WZAnalyzer::PassMuonNormChi2Cut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Chi2 Cut"<<endl;
  return (mu->globalTrack()->normalizedChi2() < maxMuonNormChi2);
}//--- PassMuonChi2Cut

bool WZAnalyzer::PassMuonHitsUsedCut(const TeVMuon* mu){
  //Num Valid Muon Hits
  if(debugme) cout<<"Check Muon Hits Used Cut"<<endl;
  return (mu->globalTrack()->hitPattern().numberOfValidMuonHits() > minMuonHitsUsed);
}//--- PassMuonHits Used Cut

bool WZAnalyzer::PassMuonStationsCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Stations Cut"<<endl;
  return (mu->numberOfMatches() > minMuonStations);
}//--- PassMuonStationsCut

bool WZAnalyzer::PassMuonEtaCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Eta Cut"<<endl;
  return (fabs(mu->eta()) < maxMuonEta);
}//--- PassMuonEta Cut

bool WZAnalyzer::PassMuonCombRelIsoCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon CombRelIso Cut"<<endl;
  return (Calc_MuonRelIso(mu) < maxWmunuCombRelIso);
}//--- PassMuonCombRelIsoCut

bool WZAnalyzer::PassMuonDxyCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Dxy Cut"<<endl;
  return (fabs(mu->userFloat("d0")) < maxMuonDxy);
}//--- PassMuonDxyCut

////////////////////////////////
/////Check TeV Properties/////
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
  return wzCand.mass("minPz") - zCand.mass() - WMASS;
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
