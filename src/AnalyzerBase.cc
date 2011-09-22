#include "UserCode/CMGWPrimeGroup/interface/AnalyzerBase.h"

using namespace std;

AnalyzerBase::AnalyzerBase(){}
AnalyzerBase::AnalyzerBase(const edm::ParameterSet & cfg, WPrimeUtil * wprimeUtil){
  //Put here everything you want EVERY CONSTRUCTOR TO DO

  wprimeUtil_ = wprimeUtil;
  assert(wprimeUtil_);

   Cuts_          = cfg.getParameter<vstring>("Cuts");
  NCuts_          = Cuts_.size();

  debugme = cfg.getParameter<bool>("debugme");

  eSelectorPset_ = cfg.getParameter<Pset>("electronSelectors");
  looseElectronType_ = cfg.getUntrackedParameter<string>("LooseElectronType", "wp95");
  tightElectronType_ = cfg.getUntrackedParameter<string>("TightElectronType", "wp95");
  looseElectron_ = ElectronSelector(eSelectorPset_, looseElectronType_);
  tightElectron_ = ElectronSelector(eSelectorPset_, tightElectronType_);
  electronResult_ = looseElectron_.getBitTemplate();
  if(debugme) cout<<"Using "<<looseElectronType_<<" for loose electrons and "
                  <<tightElectronType_<<" for tight electrons\n";

  mSelectorPset_ = cfg.getParameter<Pset>("muonSelectors");
  looseMuonType_ = cfg.getUntrackedParameter<string>("LooseMuonType", "VBTF");
  tightMuonType_ = cfg.getUntrackedParameter<string>("TightMuonType", "VBTF");
  looseMuon_ = MuonSelector(mSelectorPset_, looseMuonType_);
  tightMuon_ = MuonSelector(mSelectorPset_, tightMuonType_);
  muonResult_ = looseMuon_.getBitTemplate();
  if(debugme) cout<<"Using "<<looseMuonType_<<" for loose muons and "
                  <<tightMuonType_<<" for tight muons\n";

  jSelectorPset_ = cfg.getParameter<Pset>("jetSelectors");
  looseJetType_ = cfg.getUntrackedParameter<string>("LooseJetType", "Base");
  looseJet_ = JetSelector(jSelectorPset_, looseJetType_);
  jetResult_ = looseJet_.getBitTemplate();
  if(debugme) cout<<"Using "<<looseJetType_<<" for jets\n";

  setCandEvtFile(cfg.getParameter<string>("candEvtFile"));

  doPreselect_ = cfg.getParameter<bool>("preselect");

  electronsLabel_ = cfg.getParameter<edm::InputTag>("electrons");
  muonsLabel_ = cfg.getParameter<edm::InputTag>("muons");
  jetsLabel_ = cfg.getParameter<edm::InputTag>("jets");
  metLabel_ = cfg.getParameter<edm::InputTag>("met");
  pfCandsLabel_ = cfg.getParameter<edm::InputTag>("particleFlow");
  vertexLabel_ = cfg.getParameter<edm::InputTag>("vertexTag");

  muReconstructor_ = cfg.getParameter<int>("muonReconstructor");
  if(debugme) cout<<"Using muon algo "<<algo_desc_long[muReconstructor_]<<endl;
  useAdjustedMET_ = cfg.getParameter<bool>("useAdjustedMET");
  
  hltEventLabel_ = cfg.getParameter<edm::InputTag>("hltEventTag");

  triggersToUse_ = cfg.getParameter<vstring>("triggersToUse");

  ////////////////////Default Cuts/////////////////
  minNLeptons_ = cfg.getUntrackedParameter<uint>("minNLeptons", 0);
  maxNLeptons_ = cfg.getUntrackedParameter<uint>("maxNLeptons", 99999);
  minNTightLeptons_ = cfg.getUntrackedParameter<uint>("minNTightLeptons", 0);
  minNJets_ = cfg.getUntrackedParameter<uint>("minNJets", 0);
  
  minMET_ = cfg.getUntrackedParameter<double>("minMET", 0.);

  minZmass_ = cfg.getUntrackedParameter<double>("minZmass", 0.);
  maxZmass_ = cfg.getUntrackedParameter<double>("maxZmass", 9e9);
  minZpt_ = cfg.getUntrackedParameter<double>("minZpt", 0.);

  minVmass_ = cfg.getUntrackedParameter<double>("minVmass", 0.);
  maxVmass_ = cfg.getUntrackedParameter<double>("maxVmass", 9e9);
  minVpt_ = cfg.getUntrackedParameter<double>("minVpt", 0.);

  minWtransMass_ = cfg.getUntrackedParameter<double>("minWtransMass", 0.);
  minWpt_ = cfg.getUntrackedParameter<double>("minWpt", 0.);

}

AnalyzerBase::~AnalyzerBase(){
  outCandEvt_.close(); 
}

///////////////Utilities//////////////////

//Tabulate results after the cut has been passed
void AnalyzerBase::tabulateEvent(const int& cut_index, const float& weight){
  if(debugme) cout<<"Tabulating results for cut_index = "
                  <<cut_index<<" = "<<Cuts_[cut_index]<<endl;

  //increase the number of events passing the cuts
  hNumEvts->Fill(cut_index,weight);

  results_[cut_index].Nsurv_evt_cut_w += weight;
  results_[cut_index].Nsurv_evt_cut++;

  //fill the histograms
  fillHistos(cut_index,weight);
}//tabulateEvent

//Writing results to a txt file
void AnalyzerBase::tabulateFile(wprime::EffV& results){
  for(uint i = 0; i < results.size(); ++i){
    //calculate efficiencies
    int idx = std::max((int)i-1,0);
    float num =   results[i].Nsurv_evt_cut_w;
    float denom = results[idx].Nsurv_evt_cut_w;
    WPrimeUtil::getEff(results[i].eff, results[i].deff, num, denom);
    WPrimeUtil::getEff(results[i].eff_abs, results[i].deff_abs, num, results[0].Nsurv_evt_cut_w);
  } // loop over different cuts
}//tabulateFile

void AnalyzerBase::printFileSummary(std::vector<wprime::InputFile>::const_iterator fi, ofstream& out){ 
  out << setiosflags(std::ios::fixed) << std::setprecision(2);
  out<<"$$$$$$$$$$$$$$$$$$$$$$$ Sample: "<<fi->samplename<<" ( "<<fi->description<<" )"<<endl;
  out << " Total # of produced events for " << wprimeUtil_->getLumi_ipb() 
      << " ipb = " << fi->Nprod_evt*fi->weight << endl;
  out << " Total # of events after pre-selection for " 
      << wprimeUtil_->getLumi_ipb() << " ipb = " << fi->Nact_evt*fi->weight 
      << endl;
  float eff, deff;
  WPrimeUtil::getEff(eff, deff, fi->Nact_evt, fi->Nprod_evt);
  out << " Preselection efficiency = " << eff*100 << " % +- " << deff*100 << " %\n";
  
  for(uint i = 0; i < Cuts_.size(); ++i){
    out<<std::right<<"Cut " << std::setw(2) << i << "("
       <<std::left<< std::setw(15) << Cuts_[i]
       <<std::right << "): " <<"Evts = " << std::setw(10) << results_[i].Nsurv_evt_cut_w
       <<" (" << std::right << std::setw(7) << results_[i].Nsurv_evt_cut << ")";
    
    out << std::setw(9) <<"\tRel eff: "<<std::setw(6)<<results_[i].eff*100
        << " +/- " << std::setw(6)<<results_[i].deff*100 << "%"
        << std::setw(9) <<"\tAbs eff: "<<std::setw(6)<<results_[i].eff_abs *100
        << " +/- " << std::setw(6)<<results_[i].deff_abs*100 << "%"
        << endl;
  } // loop over different cuts
}//printFileSummary

void AnalyzerBase::resetcounters(){
  results_.assign(NCuts_,wprime::FilterEff());
}

void AnalyzerBase::clearEvtVariables(){
  allJets_.clear();
  looseJets_.clear();
  tightJets_.clear();
  allElectrons_.clear();
  looseElectrons_.clear();
  tightElectrons_.clear();
  allMuons_.clear();
  looseMuons_.clear();
  tightMuons_.clear();
  zCand_ = ZCandidate();
  vCand_ = WCandidate();
  wCand_ = WCandidate();
}

/////printers//////////////
void
AnalyzerBase::printEventFull(edm::EventBase const & event) const{
  WPrimeUtil::printEvent(event);
  WPrimeUtil::printPassingTriggers(triggerEvent_,triggersToUse_);
  printEventDetails();
  printEventLeptons();
}

void AnalyzerBase::printPassingEvent(edm::EventBase const & event){
  printEventToFile(event);
  WPrimeUtil::printEvent(event);
  printEventDetails();
}

void AnalyzerBase::printDebugEvent() const{
  WPrimeUtil::printPassingTriggers(triggerEvent_,triggersToUse_);
  printEventDetails();
  printEventLeptons();
  printLeptons();
}

void AnalyzerBase::printEventToFile(edm::EventBase const & event){
  outCandEvt_<<event.id().run()<<":"
             <<event.id().luminosityBlock()<<":"
             <<event.id().event()<<endl;
}

void AnalyzerBase::printEventDetails() const{
  if(zCand_){
    cout<<" Z Flavor: "<<zCand_.flavor()
        <<" Z Mass: "<<zCand_.mass()
        <<" Z Eta: "<<zCand_.eta()
        <<" Z Phi: "<<zCand_.phi()
        <<endl;
  }
  if(wCand_){
    cout<<" W Flavor: "<<wCand_.flavor()
        <<" W MT: "<<wCand_.mt()
        <<" pfMet et: "<<met_.et()
        <<" pfMet phi: "<<met_.phi()
        <<endl;
  }
  if(vCand_){
    cout<<" V Flavor: "<<vCand_.flavor()
        <<" V Mass: "<<vCand_.mass()
        <<" V Eta: "<<vCand_.eta()
        <<" V Phi: "<<vCand_.phi()
        <<endl;
  }
}

void
AnalyzerBase::printEventLeptons() const{
  cout<<"You shouldn't be seeing this! Implement your own.\n";
}

void
AnalyzerBase::printLeptons() const{
  printElectrons();
  printMuons();
  cout<<"----------------------\n";
}

void
AnalyzerBase::printElectrons() const{
  cout<<"----All Electrons: ( "<<allElectrons_.size()<<" )------\n";
  for(uint i=0; i<allElectrons_.size(); ++i){
    cout<<" ---  ";
    if(WPrimeUtil::Contains(allElectrons_[i], looseElectrons_)) cout<<" Loose ";
    if(WPrimeUtil::Contains(allElectrons_[i], tightElectrons_)) cout<<" Tight ";
    cout<<"  ---\n";
    printElectron(allElectrons_[i]);
  }
}

void
AnalyzerBase::printMuons() const{
  cout<<"----All Muons: ( "<<allMuons_.size()<<" )------\n";
  for(uint i=0; i<allMuons_.size(); ++i){
    cout<<" ---  ";
    if(WPrimeUtil::Contains(allMuons_[i], looseMuons_)) cout<<" Loose ";
    if(WPrimeUtil::Contains(allMuons_[i], tightMuons_)) cout<<" Tight ";
    cout<<"  ---\n";
    printMuon(allMuons_[i]);
  }
}

void
AnalyzerBase::printJets() const{
  cout<<"----All Jets: ( "<<allJets_.size()<<" )------\n";
  for(uint i=0; i<allJets_.size(); ++i){
    cout<<" ---  ";
    if(WPrimeUtil::Contains(allJets_[i], looseJets_)) cout<<" Loose ";
    cout<<"  ---\n";
    printJet(allJets_[i]);
  }
}

void
AnalyzerBase::printElectron(const heep::Ele& elec) const{
  cout << setiosflags(ios::fixed) << setprecision(2);
  cout<<" Elec ScEt: "<<elec.et()<<endl; //ScEt
  if(!elec.isPatEle()){
    cout<<"Not a pat electron, why???\n";
    return;
  }
  printElectron(elec.patEle());
}

void
AnalyzerBase::printElectron(const pat::Electron& elec) const{
  cout << setiosflags(ios::fixed) << setprecision(2);
  cout<<" Elec Pt: "<<elec.pt()<<endl
      <<" Elec P4Pt: "<<elec.p4().Pt()<<endl
      <<" Elec energy: "<<elec.energy()<<endl
      <<" Elec Charge: "<<elec.charge()<<endl
      <<" Elec Eta: "<<elec.eta()<<", isEB="<<elec.isEB()<<endl //Eta
      <<" Elec Phi: "<<elec.phi()<<endl
      <<" Elec NMiss: "<<elec.gsfTrack().get()->trackerExpectedHitsInner().numberOfHits()<<endl
      <<" Elec Dist: "<<elec.convDist()<<endl
      <<" Elec dCotTheta: "<<elec.convDcot()<<endl
      <<" Elec SigmaNN: "<<elec.sigmaIetaIeta()<<endl //sigmaNN
      <<" Elec dPhi: "<<elec.deltaPhiSuperClusterTrackAtVtx()<<endl //DeltaPhi
      <<" Elec dEta: "<<elec.deltaEtaSuperClusterTrackAtVtx()<<endl //DeltaEta
      <<" Elec HoverE: "<<elec.hadronicOverEm()<<endl// H/E
      <<" Elec EoverP: "<<elec.eSuperClusterOverP()<<endl;// E/P
}


void
AnalyzerBase::printMuon(const TeVMuon& mu) const{
  cout << setiosflags(ios::fixed) << setprecision(2);
  reco::TrackRef gm = mu.globalTrack();
  cout<<" Muon Pt: "  <<mu.pt()<<endl
      <<" Muon Global Pt: "  <<mu.pt(kGLOBAL)<<endl
      <<" Muon Inner  Pt: "  <<mu.pt(kINNER)<<endl
      <<" Muon TPFMS  Pt: "  <<mu.pt(kTPFMS)<<endl
      <<" Muon COCKTA Pt: "  <<mu.pt(kCOCKTAIL)<<endl
      <<" Muon Picky  Pt: "  <<mu.pt(kPICKY)<<endl
      <<" Muon TeV    Pt: "  <<mu.pt(kTEV)<<endl
      <<" Muon DYT    Pt: "  <<mu.pt(kDYT)<<endl
      <<" Muon Pat    Pt: "  <<mu.pt(kPAT)<<endl
      <<" Muon Charge: "<<mu.charge()<<endl
      <<" Muon Eta: " <<mu.eta()<<endl
      <<" Muon Phi: " <<mu.phi()<<endl
      <<" Muon Dxy: " <<mu.userFloat("d0")<<endl //Dxy
      <<" Muon NormX2: "<<gm->normalizedChi2()<<endl //NormX2
      <<" Muon NPix: "  <<gm->hitPattern().numberOfValidPixelHits()<<endl //Npixhit
      <<" Muon NTrk: "  <<gm->hitPattern().numberOfValidTrackerHits()<<endl //Ntrk hit
      <<" Muon NMatches: "<<mu.numberOfMatches()<<endl //MuonStations
      <<" Muon Hits: "  <<gm->hitPattern().numberOfValidMuonHits()<<endl; //Muon Hits
}

void
AnalyzerBase::printJet(const pat::Jet& jet) const{
  cout << setiosflags(ios::fixed) << setprecision(2);
  cout<<" Jet Pt: "  <<jet.pt()<<endl
      <<" Jet Mass: " <<jet.mass()<<endl
      <<" Jet Eta: " <<jet.eta()<<endl
      <<" Jet Phi: " <<jet.phi()<<endl
      <<" Jet NHF: " <<jet.neutralHadronEnergyFraction()<<endl
      <<" Jet NEF: " <<jet.neutralEmEnergyFraction()<<endl
      <<" Jet NDau: " <<jet.numberOfDaughters()<<endl
      <<" Jet CHF: " <<jet.chargedHadronEnergyFraction()<<endl
      <<" Jet CEF: " <<jet.chargedEmEnergyFraction()<<endl
      <<" Jet CMult: " <<jet.chargedMultiplicity()<<endl;
}

////////////////////////////
//////////setters///////////
////////////////////////////
void AnalyzerBase::setCandEvtFile(const std::string&  s){
  outCandEvt_.open(s.c_str());
  WPrimeUtil::CheckStream(outCandEvt_, s);
}

////////////////////
//////Cuts//////////
////////////////////
bool
AnalyzerBase::passCuts(const float& weight){
  if (debugme) cout<<"In pass Cuts\n";
  
  int iCut=0;
  if(passNoCut()) return false;
  tabulateEvent(iCut++,weight);
  
  return true;
}

inline bool AnalyzerBase::passNoCut() const{
  return true;
}

inline bool AnalyzerBase::passTriggersCut() const{
  return WPrimeUtil::passTriggersCut(triggerEvent_,triggersToUse_);
}//--- passTriggersCut

inline bool
AnalyzerBase::passMinNLeptonsCut() const{
  return (looseElectrons_.size() + looseMuons_.size()) >= minNLeptons_;
}

inline bool
AnalyzerBase::passMinNTightLeptonsCut() const{
  return (tightElectrons_.size() + tightMuons_.size()) >= minNTightLeptons_;
}
inline bool
AnalyzerBase::passMaxNLeptonsCut() const{
  return (looseElectrons_.size() + looseMuons_.size()) <= maxNLeptons_;
}

inline bool
AnalyzerBase::passMinNJetsCut() const{
  return looseJets_.size() >= minNJets_;
}

inline bool
AnalyzerBase::passMinMETCut() const{
  return met_.et() > minMET_;
}

inline bool
AnalyzerBase::passMinPtCut(const reco::Candidate& cand, const float& cut) const{
  return cand.pt() > cut;
}

////////////////////////////////
/////////Check Z Properties/////
////////////////////////////////
inline bool
AnalyzerBase::passValidZCut() const{
  return zCand_ && zCand_.mass()>0.;
}

inline bool
AnalyzerBase::passZMassCut() const{
  return (zCand_.mass() > minZmass_) && (zCand_.mass() < maxZmass_);  
}

inline bool
AnalyzerBase::passZptCut() const{
  return zCand_.pt() > minZpt_;
}

////////////////////////////////
///////Check Had. V Properties//
////////////////////////////////
inline bool
AnalyzerBase::passValidVCut() const{
  return vCand_ && vCand_.mass()>0.;
}

inline bool AnalyzerBase::passVMassCut() const{
  return (vCand_.mass() > minVmass_) && (vCand_.mass() < maxVmass_);
}

inline bool
AnalyzerBase::passVptCut() const{
  return vCand_.pt() > minVpt_;
}

/////////Check W Properties/////
inline bool
AnalyzerBase::passValidWCut() const{
  return wCand_ && wCand_.mass()>0.;
}

inline bool
AnalyzerBase::passWtransMassCut() const{
  return wCand_.mt() > minWtransMass_;
}

inline bool
AnalyzerBase::passWptCut() const{
  return wCand_.pt() > minWpt_;
}

//////////////////
//file stuff//////
//////////////////
void AnalyzerBase::beginFile(std::vector<wprime::InputFile>::const_iterator fi){
  TFileDirectory dir = wprimeUtil_->getFileService()->mkdir(fi->samplename); 
  defineHistos(dir);
  resetcounters();
}

void AnalyzerBase::defineHistos(const TFileDirectory & dir){

  string title = Form("Expected # of Events / %.0f pb^{-1}",  wprimeUtil_->getLumi_ipb());
  title = title + ";;" + title;
  hNumEvts = NULL; hNumEvts = dir.make<TH1F>("hNumEvts",title.c_str(),NCuts_,0,NCuts_);
  for(int i=0; i<NCuts_; ++i) hNumEvts->GetXaxis()->SetBinLabel(i+1,Cuts_[i].c_str());

}//defineHistos

void
AnalyzerBase::defineHistoset(const string& n, const string& t, const string& xtitle,
                              const int& nbins, const float& min, const float& max, const string& units,
                              vector<TH1F*>& h, const TFileDirectory& d){
  h.assign(NCuts_,NULL);
  
  float binWidth = (max-min)/nbins;
  for(int i=0; i<NCuts_; ++i){
    
    string name = n + "_" + Cuts_[i];
    string title = t + " (After " + Cuts_[i] + " Cut);"; 
    title += xtitle + ";Events";
    if(units.compare("NONE"))
      title += Form(" / %.0f ",binWidth) + units;
    //title += Form(" / %.0f pb^{-1}", wprimeUtil_->getLumi_ipb());
    h[i] = d.make<TH1F>(name.c_str(),title.c_str(),nbins,min,max);
  }
}

void 
AnalyzerBase::eventLoop(edm::EventBase const & event){
  clearEvtVariables();
  if(debugme) WPrimeUtil::printEvent(event);

  if (doPreselect_){
  }

  event.getByLabel(electronsLabel_,patElectronsH_);
  event.getByLabel(muonsLabel_,patMuonsH_);
  event.getByLabel(metLabel_, metH_);

  if(wprimeUtil_->DebugEvent(event)){
    cout<<"This is a debug event\n";
    printPassingEvent(event);
    printDebugEvent();
  }

  if(!passCuts(wprimeUtil_->getWeight())) return;
  if(wprimeUtil_->runningOnData()){
    cout<<" The following data event passed All Cuts!!!\n";
    printPassingEvent(event);
    cout<<" ------------------\n";
  }
  
}

// operations to be done when closing input file 
// (e.g. print summary)
void AnalyzerBase::endFile(std::vector<wprime::InputFile>::const_iterator fi,
                                 ofstream & out){
  tabulateFile(results_);
  printFileSummary(fi, out);  
}

void AnalyzerBase::endAnalysis(ofstream & out){
}

