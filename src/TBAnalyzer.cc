#include "UserCode/CMGWPrimeGroup/interface/TBAnalyzer.h"

using namespace std;

TBAnalyzer::TBAnalyzer(){}
TBAnalyzer::TBAnalyzer(const edm::ParameterSet & cfg, int fileToRun) :
  AnalyzerBase(cfg, fileToRun){
  setupCutOrder();
  if(debug_) printf("Using %i cuts\n",NCuts_);

  
// +++++++++++++++++++General Cut values
  minNLeptons_ = cfg.getUntrackedParameter<uint>("minNLeptons", 0);
  minNJets_ = cfg.getUntrackedParameter<uint>("minNJets", 0);
  
  minMET_ = cfg.getUntrackedParameter<double>("minMET", 0.);

  minWtransMass_ = cfg.getUntrackedParameter<double>("minWtransMass", 0.);
  minWpt_ = cfg.getUntrackedParameter<double>("minWpt", 0.);

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

TBAnalyzer::~TBAnalyzer(){
}

void TBAnalyzer::setupCutOrder(){
  map<string,fnCut > mFnPtrs;

  mFnPtrs["NoCuts"] = boost::bind(&TBAnalyzer::passNoCut, this);
  mFnPtrs["HLT"] = boost::bind(&TBAnalyzer::passTriggersCut, this);
  mFnPtrs["MinNLeptons"] = boost::bind(&TBAnalyzer::passMinNLeptonsCut, this, boost::cref(looseElectrons_), boost::cref(looseMuons_), boost::cref(minNLeptons_));
  mFnPtrs["MinNJets"] = boost::bind(&TBAnalyzer::passMinNJetsCut, this, boost::cref(looseJets_), boost::cref(minNJets_));
  mFnPtrs["ValidW"] = boost::bind(&TBAnalyzer::passValidWCut, this, boost::cref(wCand_));
  mFnPtrs["ValidB"] = boost::bind(&TBAnalyzer::passValidBCut, this, boost::cref(bCand1_));
  mFnPtrs["ValidT"] = boost::bind(&TBAnalyzer::passValidTCut, this, boost::ref(tCand_));
  mFnPtrs["ValidTBCand"] = boost::bind(&TBAnalyzer::passValidTBCut, this, boost::ref(tbCand_));
  mFnPtrs["WTransMass"] = boost::bind(&TBAnalyzer::passWtransMassCut, this, boost::cref(wCand_), boost::cref(minWtransMass_));
  mFnPtrs["MET"] = boost::bind(&TBAnalyzer::passMinMETCut, this, boost::cref(met_), boost::cref(minMET_));
  mFnPtrs["Wpt"] = boost::bind(&TBAnalyzer::passWptCut, this, boost::cref(wCand_), boost::cref(minWpt_));
  mFnPtrs["AllCuts"] = boost::bind(&TBAnalyzer::passNoCut, this);

  fillCuts(mFnPtrs); 
}

void TBAnalyzer::defineHistos(const TFileDirectory & dir){
}//defineHistos


void TBAnalyzer::fillHistos(const int& index, const float& weight){
}//fillHistos


void 
TBAnalyzer::eventLoop(edm::EventBase const & event){
  clearEvtVariables();
  if(debug_) WPrimeUtil::printEvent(event);
  if (doPreselect_){
    if(debug_) cout<<"Testing Preselection...\n";
  }

  event.getByLabel(electronsLabel_,patElectronsH_);
  event.getByLabel(muonsLabel_,patMuonsH_);
  if(patElectronsH_->size() + patMuonsH_->size() == 0){
    cout << "Not enough leptons. Bad bad event, returning now..." << endl;
    return;
  }
  for (size_t i = 0; i < patElectronsH_->size(); i++) {
    allElectrons_.push_back(heep::Ele((*patElectronsH_)[i]));   
    if(Overlap(allElectrons_[i].patEle(), *patMuonsH_.product(), 0.01)) continue;
    if (looseElectron_(allElectrons_[i].patEle()))
      looseElectrons_.push_back(allElectrons_[i]);
  }

  // Loop over muons, and see if they pass the TeVMuon criteria  
  for (size_t i = 0; i < patMuonsH_->size(); i++) {
    allMuons_.push_back(TeVMuon((*patMuonsH_)[i],muReconstructor_));   
    
    if (looseMuon_(allMuons_[i]) )
      looseMuons_.push_back(allMuons_[i]);
  }
  if(debug_){
    printLeptons();
    printf("    Contains: %i electron(s), %i muon(s)\n",
           (int)allElectrons_.size(), (int)allMuons_.size());
    printf("    Contains: %i loose electron(s), %i loose muon(s)\n",
           (int)looseElectrons_.size(), (int)looseMuons_.size());
  }

  //get Jets
  allJets_      = getProduct<vector<pat::Jet     > >(event, jetsLabel_);
  if(allJets_.size() < 1){
    cout << "Not enough jets. Bad bad event, returning now..." << endl;
    return;
  }
  if(debug_)
    printf("    Contains: %i pat jet(s)\n",
           (int)allJets_.size());

  // Loop over jets, and see if they pass the jet criteria
  for (size_t i = 0; i < allJets_.size(); ++i) {
    if (looseJet_(allJets_[i]) && !Overlap(allJets_[i], looseMuons_, 0.5, 2))
      looseJets_.push_back(allJets_[i]);
  }
  if(debug_)
    printf("    Contains: %i loose jet(s)\n",
           (int)looseJets_.size());
  
  //////////////////////////////////////////////
  /// Start Applying Cuts///////////////////////
  //////////////////////////////////////////////
  weight_ = wprimeUtil_->getWeight();

  int iCut=0;
  if( !passNoCut() ) return;
  tabulateEvent(iCut++, weight_);

  if( !passMinNLeptonsCut(looseElectrons_, looseMuons_, minNLeptons_) ) return;
  tabulateEvent(iCut++, weight_);

  if( !passMinNJetsCut(looseJets_, minNJets_) ) return;
  tabulateEvent(iCut++, weight_);

  //////////////////////////////
  ///////  Make B Jet  /////////
  //////////////////////////////

  //TODO: Implement here
  
  if( !passValidBCut(bCand1_) ) return;
  tabulateEvent(iCut++, weight_);

  //////////////////////////////
  ///////  Make W  /////////////
  //////////////////////////////

  //TODO: Implement here
  wCand_ = getWCand(looseElectrons_, looseMuons_, met_);

  if( !passValidWCut(wCand_) ) return;
  tabulateEvent(iCut++, weight_);

  //if( !passValidWPtCut() ) return;
  //tabulateEvent(iCut++, weight_);

  //////////////////////////////
  /////  Make Top Quark  ///////
  //////////////////////////////

  //TODO: Implement here
  //tCand_ = XWLeptonic(bCand1_, wCand_);

  if( !passValidTCut(tCand_) ) return;
  tabulateEvent(iCut++, weight_);
  
  //if( !passTMassCut() ) return;
  //tabulateEvent(iCut++, weight_);

  //if( !passTPtCut() ) return;
  //tabulateEvent(iCut++, weight_);

  //////////////////////////////
  ///////  Make 2nd B  /////////
  //////////////////////////////

  //TODO: Implement here
  //How to distingush two B's??

  if( !passValidBCut(bCand2_) ) return;
  tabulateEvent(iCut++, weight_);

  //////////////////////////////
  ///////  Make TB Cand  ///////
  //////////////////////////////

  //tbCand_ = XWLeptonic(tCand_, bCand2_);

  if( !passValidTBCut(tbCand_) ) return;
  tabulateEvent(iCut++, weight_);

  //AllCuts
  tabulateEvent(iCut++, weight_);

  if(wprimeUtil_->runningOnData()){
    cout<<" The following data event passed All Cuts!!!\n";
    printPassingEvent(event);
    cout<<" ------------------\n";
  }

}

inline void
TBAnalyzer::clearEvtVariables(){
  AnalyzerBase::clearEvtVariables();
  bCand1_ = pat::Jet();
  bCand2_ = pat::Jet();
  wCand_ = WCandidate();
  tCand_ = XWLeptonic();
  tbCand_ = XWLeptonic();
}

/////////////////Cuts///////////////////////
inline bool TBAnalyzer::passValidBCut(const pat::Jet& b) const{
  //calcBVariables();
  return b.mass()>0.;
}

inline bool TBAnalyzer::passValidTCut(XWLeptonic & t) const{
  //calcWZVariables();
  return AnalyzerBase::passValidXWCut(t);
}

/////////Check TB Properties/////
inline bool TBAnalyzer::passValidTBCut(XWLeptonic & tb) const{
  //calcWZVariables();
  return AnalyzerBase::passValidXWCut(tb);
}
