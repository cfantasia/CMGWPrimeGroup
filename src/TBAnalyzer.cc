#include "UserCode/CMGWPrimeGroup/interface/TBAnalyzer.h"

using namespace std;

TBAnalyzer::TBAnalyzer(){}
TBAnalyzer::TBAnalyzer(const edm::ParameterSet & cfg, int fileToRun) :
  AnalyzerBase(cfg, fileToRun){
  setupCutOrder();
  if(debugme) printf("Using %i cuts\n",NCuts_);

  
// +++++++++++++++++++General Cut values

}

TBAnalyzer::~TBAnalyzer(){
}

void TBAnalyzer::setupCutOrder(){
  mFnPtrs_["NoCuts"] = boost::bind(&TBAnalyzer::passNoCut, this);
  mFnPtrs_["HLT"] = boost::bind(&TBAnalyzer::passTriggersCut, this);
  mFnPtrs_["MinNLeptons"] = boost::bind(&TBAnalyzer::passMinNLeptonsCut, this, boost::cref(looseElectrons_), boost::cref(looseMuons_), boost::cref(minNLeptons_));
  mFnPtrs_["MinNJets"] = boost::bind(&TBAnalyzer::passMinNJetsCut, this, boost::cref(looseJets_), boost::cref(minNJets_));
  mFnPtrs_["ValidW"] = boost::bind(&TBAnalyzer::passValidWCut, this, boost::cref(wCand_));
  mFnPtrs_["ValidB"] = boost::bind(&TBAnalyzer::passValidBCut, this);
  mFnPtrs_["ValidT"] = boost::bind(&TBAnalyzer::passValidTCut, this);
  mFnPtrs_["ValidTBCand"] = boost::bind(&TBAnalyzer::passValidTBCut, this);
  mFnPtrs_["WTransMass"] = boost::bind(&TBAnalyzer::passWtransMassCut, this, boost::cref(wCand_), boost::cref(minWtransMass_));
  mFnPtrs_["MET"] = boost::bind(&TBAnalyzer::passMinMETCut, this, boost::cref(met_), boost::cref(minMET_));
  mFnPtrs_["Wpt"] = boost::bind(&TBAnalyzer::passWptCut, this, boost::cref(wCand_), boost::cref(minWpt_));
  mFnPtrs_["AllCuts"] = boost::bind(&TBAnalyzer::passNoCut, this);

  fillCuts(); 
}

void TBAnalyzer::defineHistos(const TFileDirectory & dir){
  AnalyzerBase::defineHistos(dir);
}//defineHistos


void TBAnalyzer::fillHistos(const int& index, const float& weight){
}//fillHistos


void 
TBAnalyzer::eventLoop(edm::EventBase const & event){
  clearEvtVariables();
  if(debugme) WPrimeUtil::printEvent(event);
  if (doPreselect_){
    if(debugme) cout<<"Testing Preselection...\n";
  }

  event.getByLabel(electronsLabel_,patElectronsH_);
  event.getByLabel(muonsLabel_,patMuonsH_);
  if(patElectronsH_->size() + patMuonsH_->size() == 0){
    cout << "Not enough leptons. Bad bad event, returning now..." << endl;
    return;
  }
  for (size_t i = 0; i < patElectronsH_->size(); i++) {
    allElectrons_.push_back(heep::Ele((*patElectronsH_)[i]));   
    /////Cory:??if(Overlap(allElectrons_[i].patEle(), *patMuonsH_.product(), 0.01)) continue;
    if (looseElectron_(allElectrons_[i].patEle(), electronLooseResult_))
      looseElectrons_.push_back(allElectrons_[i]);

    if (tightElectron_(allElectrons_[i].patEle(), electronTightResult_))
      tightElectrons_.push_back(allElectrons_[i]);
  }

  // Loop over muons, and see if they pass the TeVMuon criteria  
  for (size_t i = 0; i < patMuonsH_->size(); i++) {
    allMuons_.push_back(TeVMuon((*patMuonsH_)[i],muReconstructor_));   
    
    if (looseMuon_(allMuons_[i], muonLooseResult_) )
      looseMuons_.push_back(allMuons_[i]);
    
    if (tightMuon_(allMuons_[i], muonTightResult_) )
      tightMuons_.push_back(allMuons_[i]);
  }
  if(debugme){
    printLeptons();
    printf("    Contains: %i electron(s), %i muon(s)\n",
           (int)allElectrons_.size(), (int)allMuons_.size());
    printf("    Contains: %i loose electron(s), %i loose muon(s)\n",
           (int)looseElectrons_.size(), (int)looseMuons_.size());
    printf("    Contains: %i tight electron(s), %i tightmuon(s)\n",
           (int)tightElectrons_.size(), (int)tightMuons_.size());
  }

  //get Jets
  allJets_      = getProduct<vector<pat::Jet     > >(event, jetsLabel_);
  if(allJets_.size() < 1){
    cout << "Not enough jets. Bad bad event, returning now..." << endl;
    return;
  }
  if(debugme)
    printf("    Contains: %i pat jet(s)\n",
           (int)allJets_.size());

  // Loop over jets, and see if they pass the jet criteria
  for (size_t i = 0; i < allJets_.size(); ++i) {
    if (looseJet_(allJets_[i], jetLooseResult_) && !Overlap(allJets_[i], looseMuons_, 0.5, 2))
      looseJets_.push_back(allJets_[i]);
  }
  if(debugme)
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
  
  if( !passValidBCut() ) return;
  tabulateEvent(iCut++, weight_);

  //////////////////////////////
  ///////  Make W  /////////////
  //////////////////////////////

  //TODO: Implement here
  wCand_ = getWCand(tightElectrons_, tightMuons_, met_);

  if( !passValidWCut(wCand_) ) return;
  tabulateEvent(iCut++, weight_);

  //if( !passValidWPtCut() ) return;
  //tabulateEvent(iCut++, weight_);

  //////////////////////////////
  /////  Make Top Quark  ///////
  //////////////////////////////

  //TODO: Implement here
  //tCand_ = XWLeptonic(bCand1_, wCand_);

  if( !passValidTCut() ) return;
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

  if( !passValidBCut() ) return;
  tabulateEvent(iCut++, weight_);

  //////////////////////////////
  ///////  Make TB Cand  ///////
  //////////////////////////////

  //tbCand_ = XWLeptonic(tCand_, bCand2_);

  if( !passValidTBCut() ) return;
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
bool
TBAnalyzer::passCuts(const float& weight){
  for(int i=0; i<NCuts_; ++i){
    if(!CutFns_[i]()) return false;
    tabulateEvent(i,weight); 
  }
  return true;
  
}

inline bool TBAnalyzer::passValidBCut() const{
  return bCand1_.mass()>0.;
}

inline bool TBAnalyzer::passValidTCut() const{
  return tCand_().mass()>0.;
}


/////////Check TB Properties/////
inline bool TBAnalyzer::passValidTBCut() const{
  return tbCand_().mass()>0.;
}
