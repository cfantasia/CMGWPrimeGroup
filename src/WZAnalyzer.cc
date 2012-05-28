#include "UserCode/CMGWPrimeGroup/interface/WZAnalyzer.h"

using namespace std;

WZAnalyzer::WZAnalyzer(){}
WZAnalyzer::WZAnalyzer(const edm::ParameterSet & cfg, int fileToRun) :
  AnalyzerBase(cfg, fileToRun){
  setupCutOrder();
  if(debug_) printf("Using %i cuts\n",NCuts_);

  wzAlgo_ = (NuAlgos)cfg.getUntrackedParameter<int>("NuAlgo", kMinPz);
  doSystematics_ = cfg.getUntrackedParameter<bool>("doSystematics", false);


  rhoFastJetLabel_ = cfg.getParameter<edm::InputTag>("rhoFastJet");
  effectiveElecArea_ = cfg.getParameter<vector<double> >("effectiveElecArea");
  effectiveMuonArea_ = cfg.getParameter<vector<double> >("effectiveMuonArea");
  
// +++++++++++++++++++General Cut values
  maxNumZs_ = cfg.getParameter<uint>("maxNumZs");
  minLeadPt_ = cfg.getParameter<double>("minLeadPt");
  minNLeptons_ = cfg.getUntrackedParameter<uint>("minNLeptons", 0);
  maxNVLLeptons_ = cfg.getUntrackedParameter<uint>("maxNVLLeptons", 999);
  minNTightLeptons_ = cfg.getUntrackedParameter<uint>("minNTightLeptons", 0);
  maxNJets_ = cfg.getUntrackedParameter<uint>("maxNJets", 999);
 
  minMET_ = cfg.getUntrackedParameter<double>("minMET", 0.);

// +++++++++++++++++++Lt Cuts
  minLt_ = cfg.getUntrackedParameter<double>("minLt", -1);

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
  string vlElectronType = cfg.getUntrackedParameter<string>("VLElectronType", "wp95");
  if(debug_) cout<<"Using "<<vlElectronType<<" for vl electrons "<<endl;
  vlElectron_ = ElectronSelector(eSelectorPset, vlElectronType);
  string looseElectronType = cfg.getUntrackedParameter<string>("LooseElectronType", "wp95");
  if(debug_) cout<<"Using "<<looseElectronType<<" for loose electrons "<<endl;
  looseElectron_ = ElectronSelector(eSelectorPset, looseElectronType);
  string tightElectronType = cfg.getUntrackedParameter<string>("TightElectronType", "wp95");
  if(debug_) cout<<"Using "<<tightElectronType<<" for tight electrons "<<endl;
  tightElectron_ = ElectronSelector(eSelectorPset, tightElectronType);

  Pset mSelectorPset = cfg.getParameter<Pset>("muonSelectors");
  string vlMuonType = cfg.getUntrackedParameter<string>("VLMuonType", "exotica");
  if(debug_) cout<<"Using "<<vlMuonType<<" for vl muons"<<endl;
  vlMuon_ = MuonSelector(mSelectorPset, vlMuonType);
  string looseMuonType = cfg.getUntrackedParameter<string>("LooseMuonType", "exotica");
  if(debug_) cout<<"Using "<<looseMuonType<<" for loose muons"<<endl;
  looseMuon_ = MuonSelector(mSelectorPset, looseMuonType);
  string tightMuonType = cfg.getUntrackedParameter<string>("TightMuonType", "exotica");
  if(debug_) cout<<"Using "<<tightMuonType<<" for tight muons"<<endl;
  tightMuon_ = MuonSelector(mSelectorPset, tightMuonType);

  Pset jSelectorPset = cfg.getParameter<Pset>("jetSelectors");
  string looseJetType = cfg.getUntrackedParameter<string>("LooseJetType", "Base");
  if(debug_) cout<<"Using "<<looseJetType<<" for jets\n";
  looseJet_ = JetSelector(jSelectorPset, looseJetType);
}

WZAnalyzer::~WZAnalyzer(){
}

void WZAnalyzer::setupCutOrder(){
  map<string,fnCut > mFnPtrs;
  mFnPtrs["NoCuts"] = boost::bind(&WZAnalyzer::passNoCut, this);
  mFnPtrs["HLT"] = boost::bind(&WZAnalyzer::passTriggersCut, this);
  mFnPtrs["MinNLeptons"] = boost::bind(&WZAnalyzer::passMinNLeptonsCut, this, boost::cref(looseElectrons_), boost::cref(looseMuons_), boost::cref(minNLeptons_));
  mFnPtrs["MaxNJets"] = boost::bind(&WZAnalyzer::passMaxNJetsCut, this, boost::cref(looseJets_), boost::cref(maxNJets_));
  mFnPtrs["MaxNVLLeptons"] = boost::bind(&WZAnalyzer::passMaxNLeptonsCut, this, boost::cref(vlElectrons_), boost::cref(vlMuons_), boost::cref(maxNVLLeptons_));
  mFnPtrs["MinNTightLeptons"] = boost::bind(&WZAnalyzer::passMinNLeptonsCut, this, boost::cref(tightElectrons_), boost::cref(tightMuons_), boost::cref(minNTightLeptons_));
  mFnPtrs["ValidW"] = boost::bind(&WZAnalyzer::passValidWCut, this, boost::ref(wCand_));
  mFnPtrs["ValidZ"] = boost::bind(&WZAnalyzer::passValidZCut, this, boost::ref(zCand_));
  mFnPtrs["ValidWZCand"] = boost::bind(&WZAnalyzer::passValidWZCut, this, boost::ref(wzCand_));
  mFnPtrs["LeadLepPt"] = boost::bind(&WZAnalyzer::passLeadingLeptonPtCut, this);
  mFnPtrs["NumZs"] = boost::bind(&WZAnalyzer::passNumberOfZsCut, this);
  mFnPtrs["ZMass"] = boost::bind(&WZAnalyzer::passZMassCut, this, boost::cref(zCand_), boost::cref(minZmass_), boost::cref(maxZmass_));
  mFnPtrs["WTransMass"] = boost::bind(&WZAnalyzer::passWtransMassCut, this, boost::cref(wCand_), boost::cref(minWtransMass_));
  mFnPtrs["MET"] = boost::bind(&WZAnalyzer::passMinMETCut, this, boost::cref(met_), boost::cref(minMET_));
  mFnPtrs["Lt"] = boost::bind(&WZAnalyzer::passLtCut, this);
  mFnPtrs["Zpt"] = boost::bind(&WZAnalyzer::passZptCut, this, boost::cref(zCand_), boost::cref(minZpt_));
  mFnPtrs["Wpt"] = boost::bind(&WZAnalyzer::passWptCut, this, boost::cref(wCand_), boost::cref(minWpt_));
  mFnPtrs["AllCuts"] = boost::bind(&WZAnalyzer::passNoCut, this);

  mFnPtrs["WLepTight"] = boost::bind(&WZAnalyzer::passWLepTightCut, this);
  mFnPtrs["ZLep1Tight"] = boost::bind(&WZAnalyzer::passZLepTightCut, this, true);
  mFnPtrs["ZLep2Tight"] = boost::bind(&WZAnalyzer::passZLepTightCut, this, false);
  mFnPtrs["WFlavorElec"] = boost::bind(&WZAnalyzer::passWFlavorElecCut, this);
  mFnPtrs["WFlavorMuon"] = boost::bind(&WZAnalyzer::passWFlavorMuonCut, this);
  mFnPtrs["FakeEvt"] = boost::bind(&WZAnalyzer::passFakeEvtCut, this);
  mFnPtrs["FakeLepProbe"] = boost::bind(&WZAnalyzer::passFakeLeptonProbeCut, this);

  fillCuts(mFnPtrs);

}

void WZAnalyzer::defineHistos(const TFileDirectory & dir){
  if(debug_) printf("Declare histos\n");

  if(!doSystematics_){
    defineHistoSet("hWZMass", "Reconstructed WZ Invariant Mass",
                   "M_{WZ} (GeV)", 250, 0, 2500, "GeV", hWZMass,dir);
    defineHistoSet("hWZ3e0mMass", "Reconstructed WZ(3e0#mu) Invariant Mass",
                   "M_{WZ}^{3e0#mu} (GeV)", 250, 0, 2500, "GeV", hWZ3e0mMass,dir);
    defineHistoSet("hWZ2e1mMass", "Reconstructed WZ(2e1#mu) Invariant Mass",
                   "M_{WZ}^{2e1#mu} (GeV)", 250, 0, 2500, "GeV", hWZ2e1mMass,dir);
    defineHistoSet("hWZ1e2mMass", "Reconstructed WZ(1e2#mu) Invariant Mass",
                   "M_{WZ}^{1e2#mu} (GeV)", 250, 0, 2500, "GeV", hWZ1e2mMass,dir);
    defineHistoSet("hWZ0e3mMass", "Reconstructed WZ(0e3#mu) Invariant Mass",
                   "M_{WZ}^{0e3#mu} (GeV)", 250, 0, 2500, "GeV", hWZ0e3mMass,dir);
  
//Q=M_{WZ} - M_W - M_Z
    defineHistoSet("hQ", "Q=M_{WZ} - M_{W} - M_{Z}",
                   "Q #equiv M_{WZ} - M_{W} - M_{Z} (GeV)", 50, 0, 500, "GeV", hQ,dir);
    defineHistoSet("hWZTransMass", "Reconstructed WZ Transverse Mass",
                   "M_{WZ}^{T} (GeV)", 250, 0, 2500, "GeV", hWZTransMass,dir);
//WZpt Histos
    defineHistoSet("hWZpt", "Reconstructed WZ Transverse Momentum",
                   "p_{WZ}^{T} (GeV)", 50, 0, 500, "GeV", hWZpt,dir);
    
//Lt Histos
    defineHistoSet("hLt", "L_{T}", 
                   "L_{T} #equiv #Sigma p_{T}^{Lep} (GeV)", 80, 0, 800, "GeV", hLt,dir);
    defineHistoSet("hTriLepMass", "hTriLepMass",
                   "Trilepton Invariant Mass", 100, 0., 1000., "GeV", hTriLepMass, dir);
  
    defineHistoSet("hEvtType", "Event Type",
                   "N_{#mu}", 4, 0, 4, "NONE", hEvtType,dir);
    defineHistoSet("hEvtTypeP", "Event Type for Q=+1",
                   "N_{#mu},W^{+}", 4, 0, 4, "NONE", hEvtTypeP,dir);
    defineHistoSet("hEvtTypeM", "Event Type for Q=-1",
                   "N_{#mu},W^{-}", 4, 0, 4, "NONE", hEvtTypeM,dir);
//Lead Lepton Pt
    defineHistoSet("hLeadPt", "Leading Lepton Pt",
                   "p_{T}^{Max}", 50, 0, 1000., "GeV", hLeadPt,dir);
    defineHistoSet("hLeadPtZee", "Leading Lepton Pt Zee",
                   "p_{T}^{Max, ee}", 50, 0, 1000., "GeV", hLeadPtZee,dir);
    defineHistoSet("hLeadPtZmm", "Leading Lepton Pt Zmm",
                   "p_{T}^{Max #mu#mu}", 50, 0, 1000., "GeV", hLeadPtZmm,dir);
    defineHistoSet("hLeadElecPt", "Leading Electron Pt",
                   "p_{T}^{Max e}", 50, 0, 1000., "GeV", hLeadElecPt,dir);
    defineHistoSet("hLeadMuonPt", "Leading Muon Pt",
                   "p_{T}^{Max #mu}", 50, 0, 1000., "GeV", hLeadMuonPt,dir);
    
///////////////////////////
//Z Mass Histos
    defineHistoSet("hZMass" , "Reconstructed Mass of Z",
                   "M_{Z} (GeV)", 30, 60, 120, "GeV", hZMass,dir);
    defineHistoSet("hZeeMass","Reconstructed Mass of Zee",
                   "M_{Z}^{ee} (GeV)", 30, 60, 120, "GeV", hZeeMass,dir);
    defineHistoSet("hZmmMass","Reconstructed Mass of Z#mu#mu",
                   "M_{Z}^{#mu#mu} (GeV)", 30, 60, 120, "GeV", hZmmMass,dir);
    defineHistoSet("hZ3e0mMass","Reconstructed Mass of Z (3e0#mu)",
                   "M_{Z}^{3e0#mu} (GeV)", 30, 60, 120, "GeV", hZ3e0mMass,dir);
    defineHistoSet("hZ2e1mMass","Reconstructed Mass of Z (2e1#mu)",
                   "M_{Z}^{2e1#mu} (GeV)", 30, 60, 120, "GeV", hZ2e1mMass,dir);
    defineHistoSet("hZ1e2mMass","Reconstructed Mass of Z (1e2#mu)",
                   "M_{Z}^{1e2#mu} (GeV)", 30, 60, 120, "GeV", hZ1e2mMass,dir);
    defineHistoSet("hZ0e3mMass","Reconstructed Mass of Z (0e3#mu)",
                   "M_{Z}^{0e3#mu} (GeV)", 30, 60, 120, "GeV", hZ0e3mMass,dir);
    
//Zpt Histos
    defineHistoSet("hZpt", "p_{T}^{Z}", 
                   "p_{T}^{Z} (GeV)", 50, 0, 1000, "GeV", hZpt,dir);
    defineHistoSet("hZeept", "p_{T}^{Z#rightarrowee}", 
                   "p_{T}^{Z#rightarrowee} (GeV)", 50, 0, 1000, "GeV", hZeept,dir);
    defineHistoSet("hZmmpt", "p_{T}^{Z#rightarrow#mu#mu}", 
                   "p_{T}^{Z#rightarrow#mu#mu} (GeV)", 50, 0, 1000, "GeV", hZmmpt,dir);

//MET Histos
    defineHistoSet("hMET", "MET",
                   "#slash{E}_{T} (GeV)", 50, 0, 500, "GeV", hMET,dir);
    defineHistoSet("hMETee", "MET",
                   "#slash{E}_{T}^{ee} (GeV)", 50, 0, 500, "GeV", hMETee,dir);
    defineHistoSet("hMETmm", "MET",
                   "#slash{E}_{T}^{#mu#mu} (GeV)", 50, 0, 500, "GeV", hMETmm,dir);
    defineHistoSet("hMET3e0m", "MET (3e0#mu)",
                   "#slash{E}_{T}^{3e0#mu} (GeV)", 50, 0, 500, "GeV", hMET3e0m,dir);
    defineHistoSet("hMET2e1m", "MET (2e1#mu)",
                   "#slash{E}_{T}^{2e1#mu} (GeV)", 50, 0, 500, "GeV", hMET2e1m,dir);
    defineHistoSet("hMET1e2m", "MET (1e2#mu)",
                   "#slash{E}_{T}^{1e2#mu} (GeV)", 50, 0, 500, "GeV", hMET1e2m,dir);
    defineHistoSet("hMET0e3m", "MET (0e3#mu)",
                   "#slash{E}_{T}^{0e3#mu} (GeV)", 50, 0, 500, "GeV", hMET0e3m,dir);
    defineHistoSet("hMETSig", "MET Significance",
                   "#slash{E}_{T}^{Signif}", 50, 0, 500, "NONE", hMETSig,dir);
    
//W Trans Mass Histos
    defineHistoSet("hWTransMass", "Reconstructed Transverse Mass of W",
                   "M_{T}^{W} (GeV)", 40, 0, 200, "GeV", hWTransMass,dir);
    defineHistoSet("hWenuTransMass", "Reconstructed Transverse Mass of We\\nu",
                   "M_{T}^{W#rightarrowe#nu} (GeV)", 40, 0, 200, "GeV", hWenuTransMass,dir);
    defineHistoSet("hWmnuTransMass", "Reconstructed TransverseMass of W#mu\\nu",
                   "M_{T}^{W#rightarrow#mu#nu} (GeV)", 40, 0, 200, "GeV", hWmnuTransMass,dir);
    defineHistoSet("hW3e0mTransMass", "Reconstructed TransverseMass of W (3e0#mu)",
                   "M_{T}^{W (3e0#mu)} (GeV)", 40, 0, 200, "GeV", hW3e0mTransMass,dir);
    defineHistoSet("hW2e1mTransMass", "Reconstructed TransverseMass of W (2e1#mu)",
                   "M_{T}^{W (2e1#mu)} (GeV)", 40, 0, 200, "GeV", hW2e1mTransMass,dir);
    defineHistoSet("hW1e2mTransMass", "Reconstructed TransverseMass of W (1e2#mu)",
                   "M_{T}^{W (1e2#mu)} (GeV)", 40, 0, 200, "GeV", hW1e2mTransMass,dir);
    defineHistoSet("hW0e3mTransMass", "Reconstructed TransverseMass of W (0e3#mu)",
                   "M_{T}^{W (0e3#mu)} (GeV)", 40, 0, 200, "GeV", hW0e3mTransMass,dir);
    
//Wpt Histos
    defineHistoSet("hWpt", "p_{T}^{W}", 
                   "p_{T}^{W} (GeV)", 50, 0, 1000, "GeV", hWpt,dir);
    defineHistoSet("hWenupt", "p_{T}^{W#rightarrowe#nu}", 
                   "p_{T}^{W#rightarrowe#nu} (GeV)", 50, 0, 1000, "GeV", hWenupt,dir);
    defineHistoSet("hWmnupt", "p_{T}^{W#rightarrow#mu#nu}", 
                   "p_{T}^{W#rightarrow#mu#nu} (GeV)", 50, 0, 1000, "GeV", hWmnupt,dir);
    
//W Charge Histos
    defineHistoSet("hWQ", "Reconstructed Charge of W",
                   "q_{W}", 2, -1.5, 1.5, "", hWQ,dir);
    defineHistoSet("hWenuQ", "Reconstructed Charge of W#rightarrowe#nu",
                   "q_{W#rightarrowe#nu}", 2, -1.5, 1.5, "", hWenuQ,dir);
    defineHistoSet("hWmnuQ", "Reconstructed Charge of W#rightarrow#mu#nu",
                   "q_{W#rightarrow#mu#nu}", 2, -1.5, 1.5, "", hWmnuQ,dir);

    defineHistoSet("hWenuCombRelIso", "Comb Rel Iso of W Electron",
                   "Electron Combined Relative Isolation", 20, 0, 0.2, "NONE", hWenuCombRelIso,dir);
    defineHistoSet("hWmnuCombRelIso", "Comb Rel Iso of W Muon",
                   "Muon Combined Relative Isolation", 20, 0, 0.2, "NONE", hWmnuCombRelIso,dir);
    
    defineHistoSet("hNLElec", "Number of Loose Electrons in Event",
                   "N_{e}^{Loose}", 10, 0, 10, "NONE", hNLElec,dir);
    defineHistoSet("hNLMuon", "Number of Loose Muons in Event",
                   "N_{#mu}^{Loose}", 10, 0, 10, "NONE", hNLMuon,dir);
    defineHistoSet("hNLLeps", "Number of Loose Leptons in Event",
                   "N_{Lep}^{Loose}", 10, 0, 10, "NONE", hNLLeps,dir);
    
    defineHistoSet("hNTElec", "Number of Tight Electrons in Event",
                   "N_{e}^{Tight}", 10, 0, 10, "NONE", hNTElec,dir);
    defineHistoSet("hNTMuon", "Number of Tight Muons in Event",
                   "N_{#mu}^{Tight}", 10, 0, 10, "NONE", hNTMuon,dir);
    defineHistoSet("hNTLeps", "Number of Tight Leptons in Event",
                   "N_{Lep}^{Tight}", 10, 0, 10, "NONE", hNTLeps,dir);
    
  
    
    hDiscriminant = dir.make<TH1F>("hDiscriminant","disc", 1e4, -1e9, 1e9);
    hDiscriminantFrac = dir.make<TH1F>("hDiscriminantFrac","Disc Frac", 1e4, -1e6, 1e6);
    hDiscriminantAngle = dir.make<TH1F>("hDiscriminantAngle","Disc Angle", 40, -2, 2);
    hDiscriminantReal = dir.make<TH1F>("hDiscriminantReal","Disc Real Part", 1e4, -1e7, 1e7);
    hDiscriminantImag = dir.make<TH1F>("hDiscriminantImag","Disc Imag Part", 1e4, -1e7, 1e7);
    hVtxMatch = dir.make<TH1F>("hVtxMatch","Mask of leptons in PV", 10, 0, 10);
    
    tEvts.assign(NCuts_,NULL);
    for(int i=0; i<NCuts_; ++i){
      string name = "tEvts_" + CutNames_[i];
      string title = "Analysis Variables (After " + CutDescs_[i] + " Cut)";
      tEvts[i] = dir.make<TTree>(name.c_str(), title.c_str());
      tEvts[i]->Branch("Run", &runNumber_);
      tEvts[i]->Branch("Lumi", &lumiNumber_);
      tEvts[i]->Branch("Event", &evtNumber_);
      tEvts[i]->Branch("WZMass", &WZMass_);
      tEvts[i]->Branch("EvtType", &evtType_);
      tEvts[i]->Branch("Lt", &Lt_);
      tEvts[i]->Branch("Zpt", &Zpt_);
      tEvts[i]->Branch("ZMass", &ZMass_);
      tEvts[i]->Branch("Wpt", &Wpt_);
      tEvts[i]->Branch("WTransMass", &WTransMass_);
      tEvts[i]->Branch("MET", &MET_);
      tEvts[i]->Branch("METSig", &METSig_);
      tEvts[i]->Branch("Discriminant", &Discriminant_);
      tEvts[i]->Branch("Q", &Q_);
      tEvts[i]->Branch("TriLepMass", &TriLepMass_);
      tEvts[i]->Branch("NVtxs", &NVtxs_);
      tEvts[i]->Branch("weight", &weight_);
      tEvts[i]->Branch("ZLep1Pt", &ZLep1Pt_);
      tEvts[i]->Branch("ZLep1Eta", &ZLep1Eta_);
      tEvts[i]->Branch("ZLep2Pt", &ZLep2Pt_);
      tEvts[i]->Branch("ZLep2Eta", &ZLep2Eta_);
      tEvts[i]->Branch("WLepPt", &WLepPt_);
      tEvts[i]->Branch("WLepEta", &WLepEta_);
    }
    
  }else{//Systematics
    //Z Mass for Tag and Probe Eff Calc
    defineHistoSet("hZeeMassTT","Reconstructed MassTT of ZeeTT",
                   "M_{Z}^{ee,TT} (GeV)", 30, 60, 120, "GeV", hZeeMassTT,dir);
    defineHistoSet("hZeeMassTF","Reconstructed Mass of ZeeTF",
                   "M_{Z}^{ee,TF} (GeV)", 30, 60, 120, "GeV", hZeeMassTF,dir);
    defineHistoSet("hZmmMassTT","Reconstructed Mass of Z#mu#muTT",
                   "M_{Z}^{#mu#mu,TT} (GeV)", 30, 60, 120, "GeV", hZmmMassTT,dir);
    defineHistoSet("hZmmMassTF","Reconstructed Mass of Z#mu#muTF",
                   "M_{Z}^{#mu#mu,TF} (GeV)", 30, 60, 120, "GeV", hZmmMassTF,dir);

    //2D ZMass plots
    defineHistoSet("hPtVsZeeBarrelMassTT","Reconstructed Mass of p_{T}^{e Barrel} vs ZeeTT",
                   "M_{Z}^{ee,TT} (GeV)", 30, 60, 120, "p_{T}^{e}", 100, 0, 100., hPtVsZeeBarrelMassTT,dir);
    defineHistoSet("hPtVsZeeBarrelMassTF","Reconstructed Mass of p_{T}^{e Barrel} vs ZeeTF",
                   "M_{Z}^{ee,TF} (GeV)", 30, 60, 120, "p_{T}^{e}", 100, 0, 100., hPtVsZeeBarrelMassTF,dir);
    defineHistoSet("hPtVsZeeEndCapMassTT","Reconstructed Mass of p_{T}^{e EndCap} vs ZeeTT",
                   "M_{Z}^{ee,TT} (GeV)", 30, 60, 120, "p_{T}^{e}", 100, 0, 100., hPtVsZeeEndCapMassTT,dir);
    defineHistoSet("hPtVsZeeEndCapMassTF","Reconstructed Mass of p_{T}^{e EndCap} vs ZeeTF",
                   "M_{Z}^{ee,TF} (GeV)", 30, 60, 120, "p_{T}^{e}", 100, 0, 100., hPtVsZeeEndCapMassTF,dir);
    defineHistoSet("hPtVsZmmMassTT","Reconstructed Mass of p_{T}^{#mu} vs Z#mu#muTT",
                   "M_{Z}^{#mu#mu,TT} (GeV)", 30, 60, 120, "p_{T}^{#mu}", 100, 0, 100., hPtVsZmmMassTT,dir);
    defineHistoSet("hPtVsZmmMassTF","Reconstructed Mass of p_{T}^{#mu} vs Z#mu#muTF",
                   "M_{Z}^{#mu#mu,TF} (GeV)", 30, 60, 120, "p_{T}^{#mu}", 100, 0, 100., hPtVsZmmMassTF,dir);

    
    defineHistoSet("hWTransMass", "Reconstructed Transverse Mass of W",
                   "M_{T}^{W} (GeV)", 40, 0, 200, "GeV", hWTransMass,dir);
    //Lead Lepton Pt
    defineHistoSet("hLeadPt", "Leading Lepton Pt",
                   "p_{T}^{Max}", 50, 0, 1000., "GeV", hLeadPt,dir);
    defineHistoSet("hLeadPtZee", "Leading Lepton Pt Zee",
                 "p_{T}^{Max, ee}", 50, 0, 1000., "GeV", hLeadPtZee,dir);
    defineHistoSet("hLeadPtZmm", "Leading Lepton Pt Zmm",
                   "p_{T}^{Max #mu#mu}", 50, 0, 1000., "GeV", hLeadPtZmm,dir);
    defineHistoSet("hLeadElecPt", "Leading Electron Pt",
                 "p_{T}^{Max e}", 50, 0, 1000., "GeV", hLeadElecPt,dir);
    defineHistoSet("hLeadMuonPt", "Leading Muon Pt",
                   "p_{T}^{Max #mu}", 50, 0, 1000., "GeV", hLeadMuonPt,dir);

    
  }//end of systematic histos
  //Eta-Pt 2D Plot
  defineHistoSet("hEtaVsPt", "#eta Vs p_{T}",
                 "p_{T}", 100, 0, 100., "#eta", 50, -2.5, 2.5, hEtaVsPt,dir);
  defineHistoSet("hEtaVsPtElec", "#eta Vs p_{T}",
                 "p_{T}^{e}", 100, 0, 100., "#eta", 50, -2.5, 2.5, hEtaVsPtElec,dir);
  defineHistoSet("hEtaVsPtMuon", "#eta Vs p_{T}",
                 "p_{T}^{#mu}", 100, 0, 100., "#eta", 50, -2.5, 2.5, hEtaVsPtMuon,dir);
  defineHistoSet("hEtaVsPt3e0m", "#eta Vs p_{T}",
                 "p_{T}^{3e0#mu}", 100, 0, 100., "#eta", 50, -2.5, 2.5, hEtaVsPt3e0m,dir);
  defineHistoSet("hEtaVsPt2e1m", "#eta Vs p_{T}",
                 "p_{T}^{2e1#mu}", 100, 0, 100., "#eta", 50, -2.5, 2.5, hEtaVsPt2e1m,dir);
  defineHistoSet("hEtaVsPt1e2m", "#eta Vs p_{T}",
                 "p_{T}^{1e2#mu}", 100, 0, 100., "#eta", 50, -2.5, 2.5, hEtaVsPt1e2m,dir);
  defineHistoSet("hEtaVsPt0e3m", "#eta Vs p_{T}",
                 "p_{T}^{0e3#mu}", 100, 0, 100., "#eta", 50, -2.5, 2.5, hEtaVsPt0e3m,dir);

  defineHistoSet("hNJets", "Number of Jets in Event",
                 "N_{Jets}", 15, 0, 15, "NONE", hNJets,dir);
  
  defineHistoSet("hNVtxs", "Number of Vertexs in Event",
                 "N_{Vtx}", 50, 0, 50, "NONE", hNVtxs,dir);
  defineHistoSet("hNVtxsZee", "Number of Vertexs in Event Z#rightarrowee",
                 "N_{Vtx}^{Z#rightarrowee}", 50, 0, 50, "NONE", hNVtxsZee,dir);
  defineHistoSet("hNVtxsZmm", "Number of Vertexs in Event Z#rightarrow#mu#mu",
                 "N_{Vtx}^{Z#rightarrow#mu#mu}", 50, 0, 50, "NONE", hNVtxsZmm,dir);
  
  defineHistoSet("hWeight", "PU Weight",
                 "Weight", 40, 0, 2, "NONE", hWeight,dir);
  defineHistoSet("hL1FastJet", "L1 Fast Jet Correction",
                 "#rho", 50, 0, 25, "NONE", hL1FastJet,dir);

}//defineHistos

//fill Histograms
void WZAnalyzer::fillHistos(const int& index, const float& weight){
  if(debug_) printf("filling Histos\n");
  if(!doSystematics_){
    if(wCand_ && zCand_){
      hWZMass[index]->Fill(WZMass_, weight);
      if     (evtType_ == 0){
        hWZ3e0mMass[index]->Fill(WZMass_, weight);
        hZ3e0mMass[index]->Fill(ZMass_, weight);
        hW3e0mTransMass[index]->Fill(wCand_.mt(), weight);
        hMET3e0m[index]->Fill(MET_, weight);
        const heep::Ele& e = *wCand_.elec();
        hEtaVsPt3e0m[index]->Fill(e.patEle().pt(), e.patEle().eta(), weight);
      }else if(evtType_ == 1){ 
        hWZ2e1mMass[index]->Fill(WZMass_, weight);
        hZ2e1mMass[index]->Fill(ZMass_, weight);
        hW2e1mTransMass[index]->Fill(wCand_.mt(), weight);
        hMET2e1m[index]->Fill(MET_, weight);
        const TeVMuon& m = *wCand_.muon();
        hEtaVsPt2e1m[index]->Fill(m.pt(), m.eta(), weight);
      }else if(evtType_ == 2){ 
        hWZ1e2mMass[index]->Fill(WZMass_, weight);
        hZ1e2mMass[index]->Fill(ZMass_, weight);
        hW1e2mTransMass[index]->Fill(wCand_.mt(), weight);
        hMET1e2m[index]->Fill(MET_, weight);
        const heep::Ele& e = *wCand_.elec();
        hEtaVsPt1e2m[index]->Fill(e.patEle().pt(), e.patEle().eta(), weight);
      }else if(evtType_ == 3){ 
        hWZ0e3mMass[index]->Fill(WZMass_, weight);
        hZ0e3mMass[index]->Fill(ZMass_, weight);
        hW0e3mTransMass[index]->Fill(wCand_.mt(), weight);
        hMET0e3m[index]->Fill(MET_, weight);
        const TeVMuon& m = *wCand_.muon();
        hEtaVsPt0e3m[index]->Fill(m.pt(), m.eta(), weight);
      }
      hQ[index]->Fill(Q_, weight); 
      hWZTransMass[index]->Fill(wzCand_().mt(), weight);
      hWZpt[index]->Fill(wzCand_().pt(), weight);
      hLt[index]->Fill(Lt_, weight);
      hTriLepMass[index]->Fill(TriLepMass_, weight);
      hEvtType[index]->Fill(evtType_, weight);
      if     (wCand_.charge() > 0) hEvtTypeP[index]->Fill(evtType_, weight);
      else if(wCand_.charge() < 0) hEvtTypeM[index]->Fill(evtType_, weight);
    }
    if(zCand_){
      hZMass[index]->Fill(ZMass_, weight);
      hZpt[index]->Fill(Zpt_, weight);
      if      (zCand_.flavor() == PDG_ID_ELEC){
        hZeeMass[index]->Fill(ZMass_, weight);
        hZeept[index]->Fill(Zpt_, weight);
        hMETee[index]->Fill(MET_, weight);
        hLeadPtZee[index]->Fill(LeadPt_, weight);
        hNVtxsZee[index]->Fill((*verticesH_).size(), weight);
      }else if (zCand_.flavor() == PDG_ID_MUON){
        hZmmMass[index]->Fill(ZMass_, weight);
        hMETmm[index]->Fill(MET_, weight);
        hZmmpt[index]->Fill(Zpt_, weight);
        hLeadPtZmm[index]->Fill(LeadPt_, weight);
        hNVtxsZmm[index]->Fill((*verticesH_).size(), weight);
      }
    }
    if(wCand_){
      hWTransMass[index]->Fill(wCand_.mt(), weight);
      hWpt[index]->Fill(wCand_.pt(), weight);
      hWQ[index]->Fill(wCand_.charge(), weight);
      if      (wCand_.flavor() == PDG_ID_ELEC){
        hWenuTransMass[index]->Fill(wCand_.mt(), weight);
        hWenupt[index]->Fill(wCand_.pt(), weight);
        hWenuQ[index]->Fill(wCand_.charge(), weight);
        const heep::Ele& e = *wCand_.elec();
        hWenuCombRelIso[index]->Fill(calcCombRelIso(e.patEle(), ElecPU(e)), weight);
      }else if (wCand_.flavor() == PDG_ID_MUON){
        hWmnuTransMass[index]->Fill(wCand_.mt(), weight);
        hWmnupt[index]->Fill(wCand_.pt(), weight);
        hWmnuQ[index]->Fill(wCand_.charge(), weight);
        const TeVMuon& m = *wCand_.muon();
        hWmnuCombRelIso[index]->Fill(m.combRelIsolation03(MuonPU(m)), weight);
      }
    }  
    
    hLeadPt[index]->Fill(LeadPt_, weight);
    hLeadElecPt[index]->Fill(LeadElecPt_, weight);
    hLeadMuonPt[index]->Fill(LeadMuonPt_, weight);
    
    hMET[index]->Fill(MET_, weight);
    hMETSig[index]->Fill(METSig_, weight);
    
    hNLElec[index]->Fill(looseElectrons_.size(), weight);
    hNLMuon[index]->Fill(looseMuons_    .size(), weight);
    hNLLeps[index]->Fill(looseElectrons_.size()+looseMuons_.size(), weight);
    
    hNTElec[index]->Fill(tightElectrons_.size(), weight);
    hNTMuon[index]->Fill(tightMuons_    .size(), weight);
    hNTLeps[index]->Fill(tightElectrons_.size()+tightMuons_.size(), weight);
    
    if(wCand_.flavor() == PDG_ID_ELEC){
      const heep::Ele& e = *wCand_.elec();
      hEtaVsPt[index]->Fill(e.patEle().pt(), e.patEle().eta(), weight);
      hEtaVsPtElec[index]->Fill(e.patEle().pt(), e.patEle().eta(), weight);
    }else if(wCand_.flavor() == PDG_ID_MUON){
      const TeVMuon& m = *wCand_.muon();
      hEtaVsPt[index]->Fill(m.pt(), m.eta(), weight);
      hEtaVsPtMuon[index]->Fill(m.pt(), m.eta(), weight);
    }
    
    if(index > 4) tEvts[index]->Fill();//trying to keep the file size down
  }else{//Systematics plots
    if(zCand_){
      if(zCand_.flavor() == PDG_ID_ELEC){
        if(TT){//Probe passed
          hZeeMassTT[index]->Fill(ZMass_, weight);
          if(zCand_.elec1()->patEle().isEB()){
            hPtVsZeeBarrelMassTT[index]->Fill(ZMass_, zCand_.daughter(0)->pt(), weight);
          }else{
            hPtVsZeeEndCapMassTT[index]->Fill(ZMass_, zCand_.daughter(0)->pt(), weight);
          }
          if(zCand_.elec2()->patEle().isEB()){
            hPtVsZeeBarrelMassTT[index]->Fill(ZMass_, zCand_.daughter(1)->pt(), weight);
          }else{
            hPtVsZeeEndCapMassTT[index]->Fill(ZMass_, zCand_.daughter(1)->pt(), weight);
          }
        }//TT
        if(TF || FT){//Probe failed 
          hZeeMassTF[index]->Fill(ZMass_, weight);
          if(FT){//Probe is daughter 0
            if(zCand_.elec1()->patEle().isEB()){
              hPtVsZeeBarrelMassTF[index]->Fill(ZMass_, zCand_.daughter(0)->pt(), weight);
            }else{
              hPtVsZeeEndCapMassTF[index]->Fill(ZMass_, zCand_.daughter(0)->pt(), weight);
            }
          }
          if(TF){//Probe is daughter 1
            if(zCand_.elec2()->patEle().isEB()){
              hPtVsZeeBarrelMassTF[index]->Fill(ZMass_, zCand_.daughter(1)->pt(), weight);
            }else{
              hPtVsZeeEndCapMassTF[index]->Fill(ZMass_, zCand_.daughter(1)->pt(), weight);
            }
          }
        }//TF || FT
      }else if(zCand_.flavor() == PDG_ID_MUON){
        if(TT){//Probe passed
          hZmmMassTT[index]->Fill(ZMass_, weight);
          hPtVsZmmMassTT[index]->Fill(ZMass_, zCand_.daughter(0)->pt(), weight/2);
          hPtVsZmmMassTT[index]->Fill(ZMass_, zCand_.daughter(1)->pt(), weight/2);
        }
        if(TF || FT){//Probe failed
          hZmmMassTF[index]->Fill(ZMass_, weight);
          if(FT){//Probe is daughter 0 
            hPtVsZmmMassTF[index]->Fill(ZMass_, zCand_.daughter(0)->pt(), weight);
          }
          if(TF){//Probe is daughter 1
            hPtVsZmmMassTF[index]->Fill(ZMass_, zCand_.daughter(1)->pt(), weight);
          }
        }//TF || FT
      }//Z lep flavor
    }//valid Z
    hLeadPt[index]->Fill(LeadPt_, weight);
    hLeadElecPt[index]->Fill(LeadElecPt_, weight);
    hLeadMuonPt[index]->Fill(LeadMuonPt_, weight);

    hNJets[index]->Fill(looseJets_.size(), weight);
    hNVtxs[index]->Fill((*verticesH_).size(), weight);
    hWeight[index]->Fill(weight_/wprimeUtil_->getSampleWeight(), 1.);//Don't weight
    hL1FastJet[index]->Fill(*rhoFastJetH_, weight);

    hWTransMass[index]->Fill(wCand_.mt(), weight);
    if(wCand_.flavor() == PDG_ID_ELEC){//remember to take the opposite
      hEtaVsPt[index]->Fill(looseMuons_[0].pt(), looseMuons_[0].eta(), weight);
      hEtaVsPtMuon[index]->Fill(looseMuons_[0].pt(), looseMuons_[0].eta(), weight);
    }else if(wCand_.flavor() == PDG_ID_MUON){
      hEtaVsPt[index]->Fill(looseElectrons_[0].patEle().pt(), looseElectrons_[0].patEle().eta(), weight);
      hEtaVsPtElec[index]->Fill(looseElectrons_[0].patEle().pt(), looseElectrons_[0].patEle().eta(), weight);
    }
  }


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
  if (debug_) cout<<"In calc Z Variables\n";
  // Reconstruct the Z
  float matchptcut = 0.;
  
  matchptcut = 8.;
  vector<bool> zElectronsMask(looseElectrons_.size(), false);
  for (size_t i=0; i < looseElectrons_.size(); i++)
    if(WPrimeUtil::passTriggerMatch(looseElectrons_[i], matchptcut, triggersToUse_))
      zElectronsMask[i] = true;

  matchptcut = 8.;
  vector<bool> zMuonsMask(looseMuons_.size(), false);
  for (size_t i=0; i < looseMuons_.size(); i++)
    if(WPrimeUtil::passTriggerMatch(looseMuons_[i], matchptcut, triggersToUse_))
      zMuonsMask[i] = true;

//  if(debug_) printf("    Contains: %lu z electron(s), %lu z muon(s)\n",
//                     zElectronsMask.size(), zMuonsMask.size());

  matchptcut = 17.;
  ZCandV zeeCands = getZCands(looseElectrons_, maxZMassDiff_, false, zElectronsMask);
  removeLowLepPtCands(zeeCands, minZeePt1_, minZeePt2_);
  for (ZCandV::iterator i = zeeCands.begin(); i != zeeCands.end(); ++i){
    const heep::Ele& e1 = *i->elec1();
    const heep::Ele& e2 = *i->elec2();
    if(!WPrimeUtil::passTriggerMatch(e1, matchptcut, triggersToUse_) &&
       !WPrimeUtil::passTriggerMatch(e2, matchptcut, triggersToUse_)){
      zeeCands.erase(i);
      i--;
    }
  }

  ZCandV zmmCands = getZCands(looseMuons_, maxZMassDiff_, false, zMuonsMask);
  removeLowLepPtCands(zmmCands, minZmmPt1_, minZmmPt2_);
  for (ZCandV::iterator i = zmmCands.begin(); i != zmmCands.end(); ++i){
    const TeVMuon& m1 = *i->muon1();
    const TeVMuon& m2 = *i->muon2();
    if(!WPrimeUtil::passTriggerMatch(m1, matchptcut, triggersToUse_) &&
       !WPrimeUtil::passTriggerMatch(m2, matchptcut, triggersToUse_)){
      zmmCands.erase(i);
      i--;
    }
  }
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
    ZCandV zmmCandsAll = getZCands(looseMuons_    , maxZMassDiff_, false);

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
      TT =   tight1 &&  tight2;
      TF = ( tight1 && !tight2);
      FT = (!tight1 &&  tight2);
    }
  }
  zCandsAll.insert(zCandsAll.begin(), zCand_);
  removeOverlapping(zCandsAll);
  numZs_ = countZCands(zCandsAll); 
  Zpt_ = zCand_.pt();
  ZMass_ = zCand_.mass();
  if(zCand_){
    ZLep1Pt_ = zCand_.daughter(0)->pt();
    ZLep1Eta_ = zCand_.daughter(0)->eta();
    ZLep2Pt_ = zCand_.daughter(1)->pt();
    ZLep2Eta_ = zCand_.daughter(1)->eta();
  }
  if(debug_){
    printf("    Contains: %lu Z candidate(s)\n", zCandsAll.size());
    for (uint i=0; i<zCandsAll.size(); ++i){
      printf("Cand: %u Flavor: %i Mass: %.2f\n", i, zCandsAll[i].flavor(), zCandsAll[i].mass());
    }
    printEventLeptons();
    printEventDetails();
  }
}

inline void
WZAnalyzer::calcWVariables(){
  if (debug_) cout<<"In calc W Variables\n";
  wCand_ = getWCand(tightElectrons_, tightMuons_, met_, zCand_, minDeltaR_);
  Wpt_ = wCand_.pt();
  WTransMass_ = wCand_.mt();
  if(wCand_){
    WLepPt_ = wCand_.daughter(0)->pt();
    WLepEta_ = wCand_.daughter(0)->eta();
  }
  if(debug_){
    printf("    Contains: %i tight W candidate(s)\n", (bool)wCand_);
    printEventLeptons(); 
    printEventDetails();
  }
}

inline void
WZAnalyzer::calcWZVariables(){
  if (debug_) cout<<"In calc WZ Variables\n";
  wzCand_ = (zCand_ && wCand_) ? XWLeptonic(zCand_, wCand_) : XWLeptonic();
  Discriminant_ = wzCand_.discriminant();
  hDiscriminant->Fill(wzCand_.discriminant(), weight_);
  hDiscriminantFrac->Fill(wzCand_.discriminantFrac(), weight_);
  hDiscriminantAngle->Fill(wzCand_.discriminantAngle(), weight_);
  hDiscriminantReal->Fill(wzCand_.discriminantReal(), weight_);
  hDiscriminantImag->Fill(wzCand_.discriminantImag(), weight_);
  WZMass_ = wzCand_(wzAlgo_).mass();
  Q_ = (zCand_ && wCand_) ? calcQ() : -999.;
  if(debug_) printEventDetails();
}

void
WZAnalyzer::calcEventVariables(){
  if (debug_) cout<<"In calc Event Variables\n";
  evtType_ = (zCand_ && wCand_) ? calcEvtType() : -999;
  if(debug_) printf("evt Type: %i, Z Flav: %i, W Flav: %i\n", evtType_, zCand_.flavor(), wCand_.flavor());
  LeadPt_ = calcLeadPt(); 
  LeadElecPt_ = calcLeadPt(PDG_ID_ELEC);
  LeadMuonPt_ = calcLeadPt(PDG_ID_MUON);
  Lt_ = (zCand_ && wCand_) ? calcLt() : -999.;
  TriLepMass_ = (zCand_ && wCand_) ? calcTriLepMass() : -999.;
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
  if(debug_){
    WPrimeUtil::printEvent(event, cout);
  }

  // Preselection - skip events that don't look promising
  if (doPreselect_){
    if(debug_) cout<<"Testing Preselection...\n";
    /*
    if (getProduct<double>(event, 
                           "wzPreselectionProducer:ZMassDiff") > 30.0 ||
        getProduct<double>(event, 
                           "wzPreselectionProducer:highestLeptonPt") < 10 ||
        getProduct<vector<uint> >(event, 
                                  "wzPreselectionProducer:nLeptonsEid")[5] < 3)
      return;
    */
  }

  // get leptons
  event.getByLabel(electronsLabel_,patElectronsH_);
  event.getByLabel(muonsLabel_,patMuonsH_);
  event.getByLabel(metLabel_, metH_);

  if(useAdjustedMET_) event.getByLabel(pfCandsLabel_, pfCandidatesH_);
  WPrimeUtil::getLeptonsMET(patElectronsH_, allElectrons_,
                            patMuonsH_, muReconstructor_, allMuons_,
                            metH_, useAdjustedMET_, met_,
                            pfCandidatesH_);
  if(debug_) printf("    Contains: %lu electron(s), %lu muon(s)\n",
                          allElectrons_.size(), allMuons_.size());

  event.getByLabel(rhoFastJetLabel_, rhoFastJetH_);

  // Make vectors of leptons passing various criteria
  PatElectronV corElectrons(allElectrons_.size());
  for (size_t i = 0; i < allElectrons_.size(); i++) {
    //adjust electron with corrections
    corElectrons[i] = allElectrons_[i].patEle();
    pat::Electron & patEle = corElectrons[i];
    float cor = patEle.userFloat("corrEt") / patEle.et();
    patEle.setP4(LorentzVector(patEle.px()*cor, patEle.py()*cor, 
                               patEle.pz()*cor, patEle.energy()*cor));
    allElectrons_[i] = heep::Ele(patEle);

    if(Overlap(allElectrons_[i].patEle(), *patMuonsH_.product(), 0.01)) continue;
    const float pu = ElecPU(allElectrons_[i]);
    if (doSystematics_ && vlElectron_(allElectrons_[i].patEle(), pu))
      vlElectrons_.push_back(allElectrons_[i]);

    if (looseElectron_(allElectrons_[i].patEle(), pu))
      looseElectrons_.push_back(allElectrons_[i]);

    if (tightElectron_(allElectrons_[i].patEle(), pu))
      tightElectrons_.push_back(allElectrons_[i]);
  }

  for (size_t i = 0; i < allMuons_.size(); i++) {
    const float pu = MuonPU(allMuons_[i]);
    if (doSystematics_ && vlMuon_(allMuons_[i], pu))
      vlMuons_.push_back(allMuons_[i]);

    if (looseMuon_(allMuons_[i], pu))
      looseMuons_.push_back(allMuons_[i]);

    if (tightMuon_(allMuons_[i], pu))
      tightMuons_.push_back(allMuons_[i]);
  }
  if(looseElectrons_.size() + looseMuons_.size() == 0) return;

  if(debug_){
    print(allElectrons_);
    print(allMuons_);
    printf("    Contains: %lu loose electron(s), %lu loose muon(s)\n",
           looseElectrons_.size(), looseMuons_.size());
    if(doSystematics_)
      printf("    Contains: %lu vl electron(s), %lu vl muon(s)\n",
             vlElectrons_.size(), vlMuons_.size());
    printf("    Contains: %lu tight electron(s), %lu tightmuon(s)\n",
           tightElectrons_.size(), tightMuons_.size());
  }

  //get Jets
  event.getByLabel(jetsLabel_, patJetsH_);
  const JetV & allJets  = *patJetsH_;
  for (size_t i = 0; i < allJets.size(); i++) {
    if (looseJet_(allJets[i]) && !Overlap(allJets[i], looseMuons_, 1.0) && !Overlap(allJets[i], looseElectrons_, 1.0))
      looseJets_.push_back(allJets[i]);
  }

  //get Trigger 
  event.getByLabel(hltEventLabel_, triggerEventH_);
  if(debug_) printDebugEvent();

  //get Vertex
  event.getByLabel(vertexLabel_, verticesH_);

  MET_ = met_.et();
  METSig_ = met_.significance();
  NVtxs_ = (*verticesH_).size();

  /*
  if(!wprimeUtil_->runningOnData()){//Don't do this for data
    if(debug_){
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
  */

  weight_ = wprimeUtil_->getWeight();
  if(!passCuts(weight_)) return;
  if(wprimeUtil_->runningOnData() && !doSystematics_){
    cout<<" The following data event passed All Cuts!!!\n";
    printPassingEvent(event);
    if(debug_) printEventLeptons();
    cout<<" ------------------\n";
  }
  if(debug_) printEventLeptons();

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
  /*
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
x  */
}

void WZAnalyzer::printDebugEvent() const{
  WPrimeUtil::printPassingTriggers(*triggerEventH_,triggersToUse_);
  printEventDetails();
  printEventLeptons();
  print(allElectrons_);
  print(allMuons_);
}

void WZAnalyzer::printEventDetails() const{
  if(zCand_){
    cout<<" Z Flavor: "<<zCand_.flavor()
        <<" Z Mass: "<<ZMass_
        <<" Z lep1 pt "<<zCand_.daughter(0)->pt()
        <<" Z lep2 pt "<<zCand_.daughter(1)->pt()
        <<endl;
  }
  if(wCand_){
    cout<<" W Flavor: "<<wCand_.flavor()
        <<" W MT: "<<WTransMass_
        <<" W lep pt "<<wCand_.daughter(0)->pt()
        <<" pfMet et: "<<MET_
        <<" pfMet phi: "<<met_.phi()
        <<endl;
  }
  if(zCand_ && wCand_ && WZMass_>0.){
    cout<<" WZ Mass: "<<WZMass_
        <<" Neu Pz: "<<wzCand_.neutrinoPz(wzAlgo_)
        <<" Lt: "<<Lt_
        <<" Zpt: "<<Zpt_
        <<" Wpt: "<<Wpt_
        <<endl;
  }
  return;
}

void
WZAnalyzer::printEventLeptons() const{
  if     (zCand_.flavor() == PDG_ID_ELEC){
    cout<<"------- Electron 1 from Z -------\n";
    print(*zCand_.elec1());
    cout<<"------- Electron 2 from Z -------\n";
    print(*zCand_.elec2());
  }else if(zCand_.flavor() == PDG_ID_MUON){
    cout<<"------- Muon 1 from Z -------\n";
    print(*zCand_.muon1());
    cout<<"------- Muon 2 from Z -------\n";
    print(*zCand_.muon2());
  }

  if     (wCand_.flavor() == PDG_ID_ELEC){   
    cout<<"------- Electron from W -------\n";
    print(*wCand_.elec());
  }else if(wCand_.flavor() == PDG_ID_MUON){
    cout<<"------- Muon from W -------\n";
    print(*wCand_.muon());
  }
}

float
WZAnalyzer::calcLeadPt(int type) const{
  if(type){
    double leadpt = -999.;
    if(type == PDG_ID_ELEC)
      for (size_t i=0; i < looseElectrons_.size(); i++)
        leadpt = TMath::Max(leadpt, looseElectrons_[i].patEle().pt());
    if(type == PDG_ID_MUON)
      for (size_t i=0; i < looseMuons_.size(); i++)
        leadpt = TMath::Max(leadpt, looseMuons_[i].pt());

    return (float)leadpt;
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

inline bool WZAnalyzer::passZLepTightCut(bool firstDaughter) const{
  if(zCand_.flavor() == PDG_ID_ELEC){
    const heep::Ele & e = firstDaughter ? *zCand_.elec1() : *zCand_.elec2();
    return WPrimeUtil::Contains(e.patEle(), tightElectrons_);
  }else if(zCand_.flavor() == PDG_ID_MUON){
    const TeVMuon & m = firstDaughter ? *zCand_.muon1() : *zCand_.muon2();
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
    return ( max(zCand_.daughter(0)->pt(),zCand_.daughter(1)->pt()) > minZeePt1_ && 
             min(zCand_.daughter(0)->pt(),zCand_.daughter(1)->pt()) > minZeePt2_ );
  }else if(zCand_.flavor() == PDG_ID_MUON){

    return ( max(zCand_.daughter(0)->pt(),zCand_.daughter(1)->pt()) > minZmmPt1_ && 
             min(zCand_.daughter(0)->pt(),zCand_.daughter(1)->pt()) > minZmmPt2_ );
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
  return wCand_.daughter(0)->pt() > minWlepPt_;
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

//Check Lt Properties
//-----------------------------------------------------------
inline bool WZAnalyzer::passLtCut() const{
//-----------------------------------------------------------
  if(debug_) cout<<"Check Lt Cuts"<<endl;
  return Lt_ > minLt_;   
}//--- passLtCut

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

//calc Lt
//-----------------------------------------------------------
inline float WZAnalyzer::calcLt() const{
  return wCand_.daughter(0)->pt() + zCand_.daughter(0)->pt() + zCand_.daughter(1)->pt();
}//--- calcLt

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
  looseJets_.clear();
  allElectrons_.clear();
  looseElectrons_.clear(); vlElectrons_.clear();
  tightElectrons_.clear();
  allMuons_.clear();
  looseMuons_.clear(); vlMuons_.clear();
  tightMuons_.clear();
  met_ = pat::MET();
  zCand_ = ZCandidate();
  wCand_ = WCandidate();
  wzCand_ = XWLeptonic();
  evtType_ = -999;
  LeadPt_ = -999;
  LeadElecPt_ = -999;
  LeadMuonPt_ = -999;
  WZMass_ = -999;
  Lt_= -999;
  TriLepMass_ = -999;
  Zpt_ = -999;
  ZMass_ = -999;
  Wpt_ = -999;
  WTransMass_ = -999;
  Q_ = -999;
  Discriminant_ = 0;
  TT = TF = FT = false;
  runNumber_ = 0;
  lumiNumber_ = 0;
  evtNumber_ = 0;
  MET_ = 0;
  METSig_ = 0;
  NVtxs_ = 0;
  ZLep1Pt_ = ZLep1Eta_ = ZLep2Pt_ = ZLep2Eta_ = WLepPt_ = WLepEta_ = 0;
  weight_ = 0;
}

inline float
  WZAnalyzer::ElecPU(const heep::Ele & e) const{
  return (*rhoFastJetH_)*effectiveElecArea_[e.patEle().isEE()];
}

inline float
  WZAnalyzer::MuonPU(const TeVMuon & m) const{
  return (*rhoFastJetH_)*effectiveMuonArea_[inEE(m)];
}

