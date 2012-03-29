#include "UserCode/CMGWPrimeGroup/interface/HadronicVWAnalyzer.h"

using namespace std;

HadronicVWAnalyzer::HadronicVWAnalyzer(){}
HadronicVWAnalyzer::HadronicVWAnalyzer(const edm::ParameterSet & cfg, int fileToRun) :
  AnalyzerBase(cfg, fileToRun){
  setupCutOrder();
  if(debug_) printf("Using %i cuts\n",NCuts_);

// +++++++++++++++++++General Cut values
  minNLeptons_ = cfg.getUntrackedParameter<uint>("minNLeptons", 0);
  minNJets_ = cfg.getUntrackedParameter<uint>("minNJets", 0);
  
  minMET_ = cfg.getUntrackedParameter<double>("minMET", 0.);

// +++++++++++++++++++V Cuts
  minVmass_ = cfg.getUntrackedParameter<double>("minVmass", 0.);
  maxVmass_ = cfg.getUntrackedParameter<double>("maxVmass", 9e9);
  minVpt_ = cfg.getUntrackedParameter<double>("minVpt", 0.);

// +++++++++++++++++++W Cuts
  minWtransMass_ = cfg.getUntrackedParameter<double>("minWtransMass", 0.);
  minWpt_ = cfg.getUntrackedParameter<double>("minWpt", 0.);
  nuAlgo_ = (NuAlgos)cfg.getUntrackedParameter<int>("NuAlgo", kMinPz);


  //Selectors
  Pset eSelectorPset = cfg.getParameter<Pset>("electronSelectors");
  string looseElectronType = cfg.getUntrackedParameter<string>("LooseElectronType", "wp95");
  looseElectron_ = ElectronSelector(eSelectorPset, looseElectronType);
  if(debug_) cout<<"Using "<<looseElectronType<<" for loose electrons\n";

  Pset mSelectorPset = cfg.getParameter<Pset>("muonSelectors");
  string looseMuonType = cfg.getUntrackedParameter<string>("LooseMuonType", "exotica");
  looseMuon_ = MuonSelector(mSelectorPset, looseMuonType);
  if(debug_) cout<<"Using "<<looseMuonType<<" for loose muons\n";

  Pset jSelectorPset = cfg.getParameter<Pset>("jetSelectors");
  string looseJetType = cfg.getUntrackedParameter<string>("LooseJetType", "Base");
  looseJet_ = JetSelector(jSelectorPset, looseJetType);
  if(debug_) cout<<"Using "<<looseJetType<<" for jets\n";

}

HadronicVWAnalyzer::~HadronicVWAnalyzer(){
}

void HadronicVWAnalyzer::setupCutOrder(){
  map<string,fnCut > mFnPtrs;

  mFnPtrs["NoCuts"] = boost::bind(&HadronicVWAnalyzer::passNoCut, this);
  mFnPtrs["HLT"] = boost::bind(&HadronicVWAnalyzer::passTriggersCut, this);
  mFnPtrs["MinNLeptons"] = boost::bind(&HadronicVWAnalyzer::passMinNLeptonsCut, this, boost::cref(looseElectrons_), boost::cref(looseMuons_), boost::cref(minNLeptons_));
  mFnPtrs["MaxNLeptons"] = boost::bind(&HadronicVWAnalyzer::passMaxNLeptonsCut, this, boost::cref(looseElectrons_), boost::cref(looseMuons_), boost::cref(minNLeptons_));
  mFnPtrs["MinNJets"] = boost::bind(&HadronicVWAnalyzer::passMinNJetsCut, this, boost::cref(looseJets_), boost::cref(minNJets_));
  mFnPtrs["ValidW"] = boost::bind(&HadronicVWAnalyzer::makeAndPassValidWCut, this);
  mFnPtrs["ValidV"] = boost::bind(&HadronicVWAnalyzer::makeAndPassValidVCut, this);
  mFnPtrs["ValidVWCand"] = boost::bind(&HadronicVWAnalyzer::makeAndPassValidVWCut, this);
  mFnPtrs["VMass"] = boost::bind(&HadronicVWAnalyzer::passZMassCut, this, boost::cref(vCand_), boost::cref(minVmass_), boost::cref(maxVmass_));
  mFnPtrs["WTransMass"] = boost::bind(&HadronicVWAnalyzer::passWtransMassCut, this, boost::cref(wCand_), boost::cref(minWtransMass_));
  mFnPtrs["MET"] = boost::bind(&HadronicVWAnalyzer::passMinMETCut, this, boost::cref(met_), boost::cref(minMET_));
  mFnPtrs["Vpt"] = boost::bind(&HadronicVWAnalyzer::passZptCut, this, boost::cref(vCand_), boost::cref(minVpt_));
  mFnPtrs["Wpt"] = boost::bind(&HadronicVWAnalyzer::passWptCut, this, boost::cref(wCand_), boost::cref(minWpt_));
  mFnPtrs["AllCuts"] = boost::bind(&HadronicVWAnalyzer::passNoCut, this);

  fillCuts(mFnPtrs);
}

//--------------------------------------------------------------
void HadronicVWAnalyzer::defineHistos(const TFileDirectory & dir){
  if(debug_) printf("Declare histos\n");

  defineHistoSet("hVWMass", "Reconstructed VW Invariant Mass",
                  "M_{VW} (GeV)", 250, 0, 2500, "GeV", hVWMass,dir);
  defineHistoSet("hVWenuMass", "Reconstructed VWe#nu Invariant Mass",
                  "M_{VW}^{e#nu} (GeV)", 250, 0, 2500, "GeV", hVWenuMass,dir);
  defineHistoSet("hVWmnuMass", "Reconstructed VW#mu#nu Invariant Mass",
                  "M_{VW}^{#mu#nu} (GeV)", 250, 0, 2500, "GeV", hVWmnuMass,dir);

//Q=M_{VW} - M_W - M_V
  defineHistoSet("hQ", "Q=M_{VW} - M_{W} - M_{V}",
                  "Q (GeV)", 50, 0, 500, "GeV", hQ,dir);
  defineHistoSet("hVWTransMass", "Reconstructed VW Transverse Mass",
                  "M_{VW}^{T} (GeV)", 100, 0, 1000, "GeV", hVWTransMass,dir);
//VWpt Histos
  defineHistoSet("hVWpt", "Reconstructed VW Transverse Momentum",
                  "p_{VW}^{T} (GeV)", 50, 0, 500, "GeV", hVWpt,dir);

  defineHistoSet("hEvtType", "Event Type",
                  "N_{#mu}", 2, 0, 2, "NONE", hEvtType,dir);
  defineHistoSet("hEvtTypeP", "Event Type for Q=+1",
                  "N_{#mu},W^{+}", 2, 0, 2, "NONE", hEvtTypeP,dir);
  defineHistoSet("hEvtTypeM", "Event Type for Q=-1",
                  "N_{#mu},W^{-}", 2, 0, 2, "NONE", hEvtTypeM,dir);

///////////////////////////
//V Mass Histos
  defineHistoSet("hVMass" , "Reconstructed Mass of V",
                  "M_{V} (GeV)", 30, 60, 120, "GeV", hVMass,dir);

//Vpt Histos
  defineHistoSet("hVpt", "p_{T}^{V}", 
                  "p_{T}^{V} (GeV)", 40, 0, 400, "GeV", hVpt,dir);
//MET Histos
  defineHistoSet("hMET", "MET",
                  "#slash{E}_{T} (GeV)", 30, 0, 300, "GeV", hMET,dir);
  defineHistoSet("hMETenu", "MET",
                  "#slash{E}_{T}^{e#nu} (GeV)", 30, 0, 300, "GeV", hMETenu,dir);
  defineHistoSet("hMETmnu", "MET",
                  "#slash{E}_{T}^{#mu#nu} (GeV)", 30, 0, 300, "GeV", hMETmnu,dir);
  defineHistoSet("hMETSig", "MET Sig",
                  "#slash{E}_{T}^{Sig} (GeV)", 30, 0, 30, "GeV", hMETSig,dir);

//W Trans Mass Histos
  defineHistoSet("hWTransMass", "Reconstructed Transverse Mass of W",
                  "M_{T} (GeV)", 20, 0, 100, "GeV", hWTransMass,dir);
  defineHistoSet("hWenuTransMass", "Reconstructed Transverse Mass of We\\nu",
                  "M_{T}^{e#nu} (GeV)", 20, 0, 100, "GeV", hWenuTransMass,dir);
  defineHistoSet("hWmnuTransMass", "Reconstructed TransverseMass of W#mu\\nu",
                  "M_{T}^{#mu#nu} (GeV)", 20, 0, 100, "GeV", hWmnuTransMass,dir);

//Wpt Histos
  defineHistoSet("hWpt", "p_{T}^{W}", 
                  "p_{T}^{W} (GeV)", 40, 0, 400, "GeV", hWpt,dir);
  defineHistoSet("hWenupt", "p_{T}^{W#rightarrowe#nu}", 
                  "p_{T}^{W#rightarrowe#nu} (GeV)", 40, 0, 400, "GeV", hWenupt,dir);
  defineHistoSet("hWmnupt", "p_{T}^{W#rightarrow#mu#nu}", 
                  "p_{T}^{W#rightarrow#mu#nu} (GeV)", 40, 0, 400, "GeV", hWmnupt,dir);

//W Charge Histos
  defineHistoSet("hWQ", "Reconstructed Charge of W",
                  "q_{W}", 3, -1, 1, "", hWQ,dir);
  defineHistoSet("hWenuQ", "Reconstructed Charge of We\\nu",
                  "q_{W}^{e#nu}", 3, -1.5, 1.5, "", hWenuQ,dir);
  defineHistoSet("hWmnuQ", "Reconstructed TransverseMass of W#mu\\nu",
                  "q_{W}^{#mu#nu}", 3, -1.5, 1.5, "", hWmnuQ,dir);

  defineHistoSet("hNLElec", "Number of Loose Electrons in Event",
                  "N_{e}^{Loose}", 10, 0, 10, "NONE", hNLElec,dir);
  defineHistoSet("hNLMuon", "Number of Loose Muons in Event",
                  "N_{#mu}^{Loose}", 10, 0, 10, "NONE", hNLMuon,dir);
  defineHistoSet("hNLLeps", "Number of Loose Leptons in Event",
                  "N_{Lep}^{Loose}", 10, 0, 10, "NONE", hNLLeps,dir);

  defineHistoSet("hNJets", "Number of Jets in Event",
                  "N_{Jets}", 20, 0, 20, "NONE", hNJets,dir);

  defineHistoSet("hNVtxs", "Number of Vertexs in Event",
                  "N_{Vtx}", 50, 0, 50, "NONE", hNVtxs,dir);

  defineHistoSet("hWeight", "PU Weight",
                 "Weight", 40, 0, 2, "NONE", hWeight,dir);

  tVWCand = dir.make<TTree>("tVWCand", "Analysis Variables after VWCand");//Only 1 for now;
  tVWCand->Branch("Run", &runNumber_);
  tVWCand->Branch("Lumi", &lumiNumber_);
  tVWCand->Branch("Event", &evtNumber_);
  tVWCand->Branch("VWMass", &VWMass_);
  tVWCand->Branch("EvtType", &evtType_);
  tVWCand->Branch("VMass", &VMass_);
  tVWCand->Branch("Vpt", &Vpt_);
  tVWCand->Branch("WTransMass", &WTransMass_);
  tVWCand->Branch("Wpt", &Wpt_);
  tVWCand->Branch("Q", &Q_);
  tVWCand->Branch("MET", &MET_);
  tVWCand->Branch("METSig", &METSig_);
  tVWCand->Branch("NVtxs", &NVtxs_);
  tVWCand->Branch("weight", &weight_);

}//defineHistos

//fill Histograms
void HadronicVWAnalyzer::fillHistos(const int& index, const float& weight){
  if(debug_) printf("filling Histos\n");
  if(wCand_ && vCand_){
    hVWMass[index]->Fill(vwCand_(nuAlgo_).mass(), weight);
    hQ[index]->Fill(Q_, weight); 
    hVWTransMass[index]->Fill(vwCand_().mt(), weight);
    hVWpt[index]->Fill(vwCand_().pt(), weight);
    hEvtType[index]->Fill(evtType_, weight);
    if     (wCand_.charge() > 0) hEvtTypeP[index]->Fill(evtType_, weight);
    else if(wCand_.charge() < 0) hEvtTypeM[index]->Fill(evtType_, weight);
    if     (wCand_.flavor() == PDG_ID_ELEC){
      hVWenuMass[index]->Fill(vwCand_(nuAlgo_).mass(), weight);      
      hWenupt[index]->Fill(wCand_.pt(), weight);
    }else if(wCand_.flavor() == PDG_ID_MUON){ 
      hVWmnuMass[index]->Fill(vwCand_(nuAlgo_).mass(), weight);      
      hWmnupt[index]->Fill(wCand_.pt(), weight);
    }
    if(CutNames_[index] == "ValidVWCand"){
      tVWCand->Fill();
    }
  }
  if(vCand_){
    hVMass[index]->Fill(vCand_.mass(), weight);
    hVpt[index]->Fill(vCand_.pt(), weight);
    if      (vCand_.flavor() == PDG_ID_ELEC){
      hMETenu[index]->Fill(met_.et(), weight);
    }else if (vCand_.flavor() == PDG_ID_MUON){
      hMETmnu[index]->Fill(met_.et(), weight);
    }

  }
  if(wCand_){
    hWTransMass[index]->Fill(wCand_.mt(), weight);
    hWpt[index]->Fill(wCand_.pt(), weight);
    hWQ[index]->Fill(wCand_.charge(), weight);
    if      (wCand_.flavor() == PDG_ID_ELEC){
      hWenuTransMass[index]->Fill(wCand_.mt(), weight);
      hWenuQ[index]->Fill(wCand_.charge(), weight);
    }else if (wCand_.flavor() == PDG_ID_MUON){
      hWmnuTransMass[index]->Fill(wCand_.mt(), weight);
      hWmnuQ[index]->Fill(wCand_.charge(), weight);
    }
  }  
  hMET[index]->Fill(met_.et(), weight);
  hMETSig[index]->Fill(met_.significance(), weight);

  hNLElec[index]->Fill(looseElectrons_.size(), weight);
  hNLMuon[index]->Fill(looseMuons_    .size(), weight);
  hNLLeps[index]->Fill(looseElectrons_.size()+looseMuons_.size(), weight);

  hNJets[index]->Fill(looseJets_.size(), weight);
  //cory: hNVtxs[index]->Fill(verticesH_->size(), weight);
  hWeight[index]->Fill(weight_/wprimeUtil_->getSampleWeight(), 1.);//Don't weight

}//fillHistos

inline void
HadronicVWAnalyzer::calcVVariables(){
  if (debug_) cout<<"In calc V Variables\n";
  vCand_ = getVCand2(looseJets_);
  VMass_ = vCand_.mass();
  Vpt_ = vCand_.pt();
  if(debug_){
    printf("    Contains: %i V candidate(s)\n", (bool)vCand_);
    printEventLeptons(); 
    printEventDetails();
  }
}

inline void
HadronicVWAnalyzer::calcWVariables(){
  if (debug_) cout<<"In calc W Variables\n";
  wCand_ = getWCand(looseElectrons_, looseMuons_, met_);
  WTransMass_ = wCand_.mt();
  Wpt_ = wCand_.pt();
  if(debug_) printf("    Contains: %i W candidate(s)\n", (bool)wCand_);
  if(debug_){
    printEventLeptons(); 
    printEventDetails();
  }
}

inline void
HadronicVWAnalyzer::calcVWVariables(){
  if (debug_) cout<<"In calc VW Variables\n";
  vwCand_ = (vCand_ && wCand_) ? XWLeptonic(vCand_, wCand_) : XWLeptonic();
  VWMass_ = vwCand_(nuAlgo_).mass();
  Q_ = (vCand_ && wCand_) ? calcQ() : -999.;
  if(debug_) printEventDetails();
}

void
HadronicVWAnalyzer::calcEventVariables(){
  if (debug_) cout<<"In calc Event Variables\n";
  evtType_ = (vCand_ && wCand_) ? calcEvtType() : -999;
  if(debug_) printf("evt Type: %i, V Flav: %i, W Flav: %i\n", evtType_, (int)vCand_.flavor(), (int)wCand_.flavor());
}

void 
HadronicVWAnalyzer::eventLoop(edm::EventBase const & event){
  clearEvtVariables();
  runNumber_ = event.id().run();
  lumiNumber_ = event.id().luminosityBlock();
  evtNumber_ = event.id().event();
  if(debug_) WPrimeUtil::printEvent(event, cout);

  // Preselection - skip events that don't look promising
  if (doPreselect_){
    if(debug_) cout<<"Testing Preselection...\n";
  }

  //get Jets
  event.getByLabel(jetsLabel_,patJetsH_);
  if(debug_) printf("    Contains: %i pat jets(s)\n",
                     (int)patJetsH_->size());
  for (size_t i = 0; i < patJetsH_->size(); i++) {
    const pat::Jet & jet = (*patJetsH_.product())[i];
    if (looseJet_(jet))
      looseJets_.push_back(jet);
  }
  if(looseJets_.size() == 0) return;


  // get leptons
  event.getByLabel(electronsLabel_,patElectronsH_);
  event.getByLabel(muonsLabel_,patMuonsH_);
  event.getByLabel(metLabel_, metH_);
  if(useAdjustedMET_) event.getByLabel(pfCandsLabel_, pfCandidatesH_);
/*
  WPrimeUtil::getMuonsMET(patMuonsH_, muReconstructor_, allMuons_,
                          metH_, useAdjustedMET_, met_,
                          pfCandidatesH_);
  
*/
  WPrimeUtil::getLeptonsMET(patElectronsH_, allElectrons_,
                            patMuonsH_, muReconstructor_, allMuons_,
                            metH_, useAdjustedMET_, met_,
                            pfCandidatesH_);
  if(debug_) printf("    Contains: %i electron(s), %i muon(s)\n",
                          (int)allElectrons_.size(), (int)allMuons_.size());
  MET_ = met_.et();
  METSig_ = met_.significance();

  // Make vectors of leptons passing various criteria
  for (size_t i = 0; i < allElectrons_.size(); i++) {
    if(Overlap(allElectrons_[i].patEle(), allMuons_)) continue;
    if (looseElectron_(allElectrons_[i].patEle()))
      looseElectrons_.push_back(allElectrons_[i]);
  }

  for (size_t i = 0; i < allMuons_.size(); i++) {
    if (looseMuon_(allMuons_[i]))
      looseMuons_.push_back(allMuons_[i]);
  }
  if(looseElectrons_.size() + looseMuons_.size() == 0) return;

  if(debug_){
    print(allElectrons_);
    print(allMuons_);
    printf("    Contains: %i loose electron(s), %i loose muon(s), %i loose jet(s)\n",
           (int)looseElectrons_.size(), (int)looseMuons_.size(), (int)looseJets_.size());
  }

  ///////////////////
  //Cory: Remove overlap btw jets and leptons! 
  //(maybe I don't have to for pfJets)
  ///////////////////


  //get Trigger 
  event.getByLabel(hltEventLabel_, triggerEventH_);

  //get Vertex
  //Cory: not yet in files
  //event.getByLabel(vertexLabel_, verticesH_);
  //NVtxs_ = (*verticesH_).size();

  if(debug_){
    cout<<"This is a debug event\n";
    printPassingEvent(event);
    printDebugEvent();
  }

  weight_ = wprimeUtil_->getWeight();
  if(!passCuts(weight_)) return;
  if(wprimeUtil_->runningOnData()){
    cout<<" The following data events passed All Cuts!!!\n";
    printPassingEvent(event);
    if(debug_) printEventLeptons();
    cout<<" ------------------\n";
  }
  if(debug_) printEventLeptons();
}

void HadronicVWAnalyzer::printDebugEvent() const{
  WPrimeUtil::printPassingTriggers(*triggerEventH_,triggersToUse_);
  printEventDetails();
  printEventLeptons();
  print(allElectrons_);
  print(allMuons_);
  print(allJets_);
}

void HadronicVWAnalyzer::printEventDetails() const{
  if(vCand_){
    cout<<" V Mass: "<<vCand_.mass()
        <<endl;
  }
  if(wCand_){
    cout<<" W Flavor: "<<wCand_.flavor()
        <<" W MT: "<<wCand_.mt()
        <<" W lep pt "<<wCand_.daughter(0)->pt()
        <<" pfMet et: "<<met_.et()
        <<" pfMet phi: "<<met_.phi()
        <<endl;
  }
  if(vCand_ && wCand_ && vwCand_(nuAlgo_).mass()>0.){
    cout<<" VW Mass: "<<vwCand_(nuAlgo_).mass()
        <<" Neu Pz: "<<vwCand_.neutrinoPz(nuAlgo_)
        <<" Vpt: "<<vCand_.pt()
        <<" Wpt: "<<wCand_.pt()
        <<endl;
  }
  return;
}

void
HadronicVWAnalyzer::printEventLeptons() const{
  print(*(pat::Jet*)vCand_.daughter(0));

  if      (wCand_.flavor() == PDG_ID_ELEC){   
    print(*wCand_.elec());
  }else if(wCand_.flavor() == PDG_ID_MUON){
    print(*wCand_.muon());
  }
}

/////////////////Accessors///////////////////////

/////////////////Modifiers///////////////////////

/////////////////Cuts///////////////////////

inline bool HadronicVWAnalyzer::makeAndPassValidVCut(){
  calcVVariables();
  return AnalyzerBase::passValidZCut(vCand_);
}

inline bool HadronicVWAnalyzer::makeAndPassValidWCut(){
  calcWVariables();
  calcEventVariables();
  return AnalyzerBase::passValidWCut(wCand_);
}

inline bool HadronicVWAnalyzer::makeAndPassValidVWCut(){
  calcVWVariables();
  return AnalyzerBase::passValidXWCut(vwCand_);
}

///////////////////////////////////

inline float HadronicVWAnalyzer::calcQ() const{
  return vwCand_(nuAlgo_).mass() - vCand_.mass() - WMASS;
}

inline int HadronicVWAnalyzer::calcEvtType() const{
  return (wCand_) ? (wCand_.flavor() != PDG_ID_ELEC) : -999;
}

///////////////Utilities//////////////////
//--------------------------------------------------------------

inline void
HadronicVWAnalyzer::clearEvtVariables(){
  allJets_.clear();
  looseJets_.clear();
  allElectrons_.clear();
  looseElectrons_.clear();
  looseMuons_.clear();
  met_ = pat::MET();
  wCand_ = WCandidate();
  vCand_ = ZCandidate();
  vwCand_ = XWLeptonic();
  runNumber_ = lumiNumber_ = evtNumber_ = 0;
  evtType_ = -999;
  VWMass_ = -999;
  VMass_ = Vpt_ = -999;
  WTransMass_  = Wpt_ = -999;
  Q_ = -999;
  MET_ = METSig_ = -999;
  NVtxs_ = -999;
  weight_ = 0;
}
