#include "UserCode/CMGWPrimeGroup/interface/AnalyzerBase.h"

using namespace std;

AnalyzerBase::AnalyzerBase(){}
AnalyzerBase::AnalyzerBase(const edm::ParameterSet & cfg, int fileToRun){
  //Put here everything you want EVERY CONSTRUCTOR TO DO

  reportAfter_ = cfg.getParameter<unsigned int>("reportAfter");
  maxEvents_   = cfg.getParameter<int>("maxEvents");
  useJSON_   = cfg.getParameter<bool>("useJSON") ;
  countGenEvts_ = cfg.getParameter<bool>("countGenEvts");
  genLabel_ = cfg.getParameter<edm::InputTag>("genParticles" );
  pfLabel_ = cfg.getParameter<edm::InputTag>("particleFlow" );
  pileupLabel_ = cfg.getParameter<edm::InputTag>("pileupTag" );
  doRecoilCorrectionForW_ = cfg.getParameter<bool>("doRecoilCorrectionForW");

  string sample_cross_sections = cfg.getParameter<string>("sample_cross_sections");
  edm::ParameterSet const& inputs = cfg.getParameter<edm::ParameterSet>("inputs");
  if ( inputs.exists("lumisToProcess") ) 
    {
      vector<edm::LuminosityBlockRange> const & lumisTemp =
	inputs.getUntrackedParameter<vector<edm::LuminosityBlockRange> > ("lumisToProcess");
      jsonVector.resize( lumisTemp.size() );
      copy( lumisTemp.begin(), lumisTemp.end(), jsonVector.begin() );
    }

  ctrNames_ = (cfg.getParameter<vstring>("eventCounters"));

  MCPUDistFile_   = cfg.getParameter<string>("MCPUDistFile" );
  MCPUDistHist_   = cfg.getParameter<string>("MCPUDistHist" );
  DataPUDistFile_ = cfg.getParameter<string>("DataPUDistFile" );
  DataPUDistHist_ = cfg.getParameter<string>("DataPUDistHist" );
  
  std::vector<edm::EventID> vEventsToDebug = cfg.getParameter<std::vector<edm::EventID> >("vEventsToDebug");
  
  wprimeUtil_ = new WPrimeUtil(outputFile_.c_str(), genLabel_, pfLabel_, sample_cross_sections);
  wprimeUtil_->setLumiWeights(MCPUDistFile_, DataPUDistFile_, MCPUDistHist_, DataPUDistHist_);
  wprimeUtil_->setEventsToDebug(vEventsToDebug);
  assert(wprimeUtil_);

  wprimeUtil_->getInputFiles(inputFiles_);
  if(fileToRun != -1){
    if(fileToRun < (int)inputFiles_.size()){
      inputFiles_.assign(1,inputFiles_[fileToRun]);
    }else{
      cerr<<"You asked for sample "<<fileToRun
          <<" but only "<<inputFiles_.size()
          <<" are listed!\n";
      inputFiles_.clear();
      abort();
    }
  }

  //////////////////////

   CutNames_ = cfg.getParameter<vstring>("Cuts");
   CutDescs_ = CutNames_;
  NCuts_     = CutNames_.size();

  debugme = cfg.getParameter<bool>("debugme");

  eSelectorPset_ = cfg.getParameter<Pset>("electronSelectors");
  looseElectronType_ = cfg.getUntrackedParameter<string>("LooseElectronType", "wp95");
  tightElectronType_ = cfg.getUntrackedParameter<string>("TightElectronType", "wp95");
  looseElectron_ = ElectronSelector(eSelectorPset_, looseElectronType_);
  tightElectron_ = ElectronSelector(eSelectorPset_, tightElectronType_);
  electronLooseResult_ = looseElectron_.getBitTemplate();
  electronTightResult_ = tightElectron_.getBitTemplate();
  if(debugme) cout<<"Using "<<looseElectronType_<<" for loose electrons and "
                  <<tightElectronType_<<" for tight electrons\n";

  mSelectorPset_ = cfg.getParameter<Pset>("muonSelectors");
  looseMuonType_ = cfg.getUntrackedParameter<string>("LooseMuonType", "VBTF");
  tightMuonType_ = cfg.getUntrackedParameter<string>("TightMuonType", "VBTF");
  looseMuon_ = MuonSelector(mSelectorPset_, looseMuonType_);
  tightMuon_ = MuonSelector(mSelectorPset_, tightMuonType_);
  muonLooseResult_ = looseMuon_.getBitTemplate();
  muonTightResult_ = tightMuon_.getBitTemplate();
  if(debugme) cout<<"Using "<<looseMuonType_<<" for loose muons and "
                  <<tightMuonType_<<" for tight muons\n";

  jSelectorPset_ = cfg.getParameter<Pset>("jetSelectors");
  looseJetType_ = cfg.getUntrackedParameter<string>("LooseJetType", "Base");
  looseJet_ = JetSelector(jSelectorPset_, looseJetType_);
  jetLooseResult_ = looseJet_.getBitTemplate();
  if(debugme) cout<<"Using "<<looseJetType_<<" for jets\n";

  //////Output Files
  outputFile_  = cfg.getParameter<string  >("outputFile" );
  logFile_  = cfg.getParameter<string  >("logFile" );
  candEvtFile_ = cfg.getParameter<string>("candEvtFile");

  if(fileToRun != -1){
    outputFile_ = Form("Sample%i_%s",fileToRun,outputFile_.c_str()); 
    logFile_ = Form("Sample%i_%s",fileToRun,logFile_.c_str()); 
    candEvtFile_ = Form("Sample%i_%s",fileToRun,candEvtFile_.c_str()); 
  }

  fs = new fwlite::TFileService(outputFile_);
  outLogFile_.open(logFile_.c_str());
  WPrimeUtil::CheckStream(outLogFile_, logFile_);
  outCandEvtFile_.open(candEvtFile_.c_str());
  WPrimeUtil::CheckStream(outCandEvtFile_, candEvtFile_);

  //////////////
  doPreselect_ = cfg.getParameter<bool>("preselect");

  electronsLabel_ = cfg.getParameter<edm::InputTag>("electrons");
  muonsLabel_ = cfg.getParameter<edm::InputTag>("muons");
  jetsLabel_ = cfg.getParameter<edm::InputTag>("jets");
  metLabel_ = cfg.getParameter<edm::InputTag>("met");
  pfCandsLabel_ = cfg.getParameter<edm::InputTag>("particleFlow");
  vertexLabel_ = cfg.getParameter<edm::InputTag>("vertexTag");

  muReconstructor_ = cfg.getParameter<int>("muonReconstructor");
  if(debugme) cout<<"Using muon algo "<<algo_desc_long[muReconstructor_]<<endl;
  assert(muReconstructor_ < Num_MuTeVtrkAlgos);

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

  muon_reconstructors.push_back(kGLOBAL);
  muon_reconstructors.push_back(kINNER);
  muon_reconstructors.push_back(kSTANDALONE);
  muon_reconstructors.push_back(kCOCKTAIL);
  muon_reconstructors.push_back(kDYT);
  muon_reconstructors.push_back(kTPFMS);
  muon_reconstructors.push_back(kPICKY);
  //  muon_reconstructors.push_back(kPAT); <- this fails for non-Cory's pat-tuples
}

AnalyzerBase::~AnalyzerBase(){
  if(wprimeUtil_) delete wprimeUtil_;
  delete fs; 
  outLogFile_.close(); 
  outCandEvtFile_.close(); 
}

///////////////Utilities//////////////////

//Fill Vector of Cuts based on map
void AnalyzerBase::fillCuts(){
  CutFns_.resize(NCuts_);
  for(int i=0; i<NCuts_; ++i){
    if(mFnPtrs_.find(CutNames_[i]) == mFnPtrs_.end()){
      cout<<"Didn't find cut named "<<CutNames_[i]<<endl;
      abort();
    }
    CutFns_[i] = mFnPtrs_.find(CutNames_[i])->second;
  } 
}

//Tabulate results after the cut has been passed
void AnalyzerBase::tabulateEvent(const int& cut_index, const float& weight){
  if(debugme) cout<<"Tabulating results for cut_index = "
                  <<cut_index<<" = "<<CutNames_[cut_index]<<endl;

  //increase the number of events passing the cuts
  hNumEvts->Fill(cut_index,weight);

  std::vector<wprime::InputFile>::iterator file = wprimeUtil_->getCurrentSample();
  file->results[cut_index].Nsurv_evt_cut_w += weight;
  file->results[cut_index].Nsurv_evt_cut++;

  //fill the histograms
  fillHistos(cut_index,weight);
}//tabulateEvent

//Writing results to a txt file
void AnalyzerBase::tabulateFile(std::vector<wprime::InputFile>::const_iterator fi, wprime::EffV& results){
  float nProd = fi->Nprod_evt*fi->weight;
  for(uint i = 0; i < results.size(); ++i){
    //calculate efficiencies
    float rel_denom = i==0 ? fi->Nact_evt*fi->weight : results[i-1].Nsurv_evt_cut_w;
    //Calc Relative Eff
    WPrimeUtil::getEff(results[i].eff, results[i].deff, 
                       results[i].Nsurv_evt_cut_w, rel_denom);
    //Calc Absolute Eff
    WPrimeUtil::getEff(results[i].eff_abs, results[i].deff_abs, 
                       results[i].Nsurv_evt_cut_w, nProd);
  } // loop over different cuts
}//tabulateFile

// print summary of efficiencies
void AnalyzerBase::printFileSummary(std::vector<wprime::InputFile>::const_iterator fi, ofstream& out){ 
  out<<"$$$$$$$$$$$$$$$$ Sample: "<<fi->samplename<<" ( "<<fi->description<<" )"
     <<". Int. Lumi = "<<wprimeUtil_->getLumi_ipb()<< ". Weight = "<<fi->weight<<endl;
  out << setiosflags(std::ios::fixed) << std::setprecision(2);
  float eff, deff;
  out<<"      (Produced       ): " <<"Evts = " << std::setw(10) << fi->Nprod_evt*fi->weight
     <<" (" << std::right << std::setw(7) << fi->Nprod_evt << ")";
  WPrimeUtil::getEff(eff, deff, fi->Nprod_evt, fi->Nprod_evt);
  out << std::setw(9) <<"\tRel eff: "<<std::setw(6)<<eff*100
      << " +/- " << std::setw(6)<<deff*100 << "%";
  out << std::setw(9) <<"\tAbs eff: "<<std::setw(6)<<eff*100
      << " +/- " << std::setw(6)<<deff*100 << "%"
      << endl;

  out<<"      (Pat Skim       ): " <<"Evts = " << std::setw(10) << fi->Nact_evt*fi->weight
     <<" (" << std::right << std::setw(7) << fi->Nact_evt << ")";
  WPrimeUtil::getEff(eff, deff, fi->Nact_evt, fi->Nprod_evt);
  out << std::setw(9) <<"\tRel eff: "<<std::setw(6)<<eff*100
      << " +/- " << std::setw(6)<<deff*100 << "%";
  out << std::setw(9) <<"\tAbs eff: "<<std::setw(6)<<eff*100
      << " +/- " << std::setw(6)<<deff*100 << "%"
      << endl;

  for(uint i = 0; i < CutNames_.size(); ++i){
    const wprime::FilterEff & curCut = fi->results[i];
    out<<std::right<<"Cut " << std::setw(2) << i << "("
       <<std::left<< std::setw(15) << CutNames_[i]
       <<std::right << "): " <<"Evts = " << std::setw(10) << curCut.Nsurv_evt_cut_w
       <<" (" << std::right << std::setw(7) << curCut.Nsurv_evt_cut << ")";
    
    out << std::setw(9) <<"\tRel eff: "<<std::setw(6)<<curCut.eff*100
        << " +/- " << std::setw(6)<<curCut.deff*100 << "%"
        << std::setw(9) <<"\tAbs eff: "<<std::setw(6)<<curCut.eff_abs *100
        << " +/- " << std::setw(6)<<curCut.deff_abs*100 << "%"
        << endl;
  } // loop over different cuts
}//printFileSummary

void AnalyzerBase::resetCounters(){
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
  outCandEvtFile_<<event.id().run()<<":"
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

void AnalyzerBase::printMuon(const TeVMuon& mu) const{
  cout << setiosflags(ios::fixed) << setprecision(2);
  
  cout << " Muon eta = " << mu.eta() << "  phi = " << mu.phi() << endl;

  typedef std::vector<unsigned>::const_iterator It;
  for(It it = muon_reconstructors.begin(); it != muon_reconstructors.end(); ++it)
    mu.printPtInfo(*it);
  
  mu.printTrackerInfo();

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

////////////////////
//////Cuts//////////
////////////////////
bool
AnalyzerBase::passCuts(const float& weight){
  if (debugme) cout<<"In pass Cuts\n";
  
  for(int i=0; i<NCuts_; ++i){
    if(!CutFns_[i]()) return false;
    tabulateEvent(i,weight); 
  }
  return true;
}

inline bool AnalyzerBase::passNoCut() const{
  return true;
}

inline bool AnalyzerBase::passTriggersCut() const{
  return WPrimeUtil::passTriggersCut(triggerEvent_,triggersToUse_);
}//--- passTriggersCut

inline bool
AnalyzerBase::passMinNLeptonsCut(const ElectronV& electrons, const MuonV& muons, const float & cut) const{
  return (electrons.size() + muons.size()) >= cut;
}

inline bool
AnalyzerBase::passMaxNLeptonsCut(const ElectronV& electrons, const MuonV& muons, const float & cut) const{
  return (electrons.size() + muons.size()) <= cut;
}

inline bool
AnalyzerBase::passMinNJetsCut(const JetV& jets, const float & cut) const{
  return jets.size() >= cut;
}

inline bool
AnalyzerBase::passMinMETCut(const pat::MET & met, const float& cut) const{
  return met.et() > cut;
}

inline bool
AnalyzerBase::passMinPtCut(const reco::Candidate& cand, const float& cut) const{
  return cand.pt() > cut;
}

////////////////////////////////
/////////Check Z Properties/////
////////////////////////////////
inline bool
AnalyzerBase::passValidZCut(const ZCandidate& z) const{
  return z && z.mass()>0.;
}

inline bool
AnalyzerBase::passZMassCut(const ZCandidate& z, const float& mincut, const float& maxcut) const{
  return (z.mass() > mincut) && (z.mass() < maxcut);  
}

inline bool
AnalyzerBase::passZptCut(const ZCandidate& z, const float& cut) const{
  return z.pt() > cut;
}

/////////Check W Properties/////
inline bool
AnalyzerBase::passValidWCut(const WCandidate& w) const{
  return w && w.mt()>0.;
}

inline bool
AnalyzerBase::passWtransMassCut(const WCandidate& w, const float& cut) const{
  return w.mt() > cut;
}

inline bool//For HadV which I dumbly made a W, to be removed
AnalyzerBase::passVMassCut(const WCandidate& w, const float& mincut, const float& maxcut) const{
  return (w.mass() > mincut) && (w.mass() < maxcut);  
}

inline bool
AnalyzerBase::passWptCut(const WCandidate& w, const float& cut) const{
  return w.pt() > cut;
}

//////////////////
//file stuff//////
//////////////////

// operations to be done when changing input file (e.g. create new histograms)
void AnalyzerBase::beginFile(std::vector<wprime::InputFile>::iterator fi){
  bool shouldCorrectMt = 
    ((fi->samplename=="W" || fi->samplename=="Wlowpt") 
     && doRecoilCorrectionForW_);
  wprimeUtil_->setApplyHadronicRecoilCorrection(shouldCorrectMt);

  wprimeUtil_->setSampleName(fi->samplename);
  wprimeUtil_->setSampleWeight(fi->weight);
  wprimeUtil_->setRunningOnData();
  wprimeUtil_->setIsSignalSample(fi->isSignal());
  wprimeUtil_->setCurrentSample(fi);
  wprimeUtil_->resetWarnings();

  if(wprimeUtil_->runningOnData())
    // Nprod_evt presumably contains the # of events before any filtering
    // that results in Nact_evt (< Nprod_evt) events contained in the file.
    // For data, we tend not to know how many events we started with,
    // so just assume pre-selection efficiency = 100%;
    // this affects only the efficiency calculations printed
    // at the end of the job - nothing else!
    fi->Nprod_evt = fi->Nact_evt;

  TFileDirectory dir = fs->mkdir(fi->samplename); 
  defineHistos(dir);
  defineResolutionHistos(dir, fi->signalMass);
  resetCounters();
  fi->results.assign(NCuts_,wprime::FilterEff());
}

void AnalyzerBase::defineHistos(const TFileDirectory & dir){

  string title = Form("Expected # of Events / %.0f pb^{-1}",  wprimeUtil_->getLumi_ipb());
  title = title + ";;" + title;
  hNumEvts = NULL; hNumEvts = dir.make<TH1F>("hNumEvts",title.c_str(),NCuts_,0,NCuts_);
  for(int i=0; i<NCuts_; ++i) hNumEvts->GetXaxis()->SetBinLabel(i+1,CutNames_[i].c_str());

}//defineHistos

void
AnalyzerBase::defineHistoset(const string& n, const string& t, 
			     const string& xtitle,
			     int nbins, float min, float max, 
			     const string& units,
			     vector<TH1F*>& h, const TFileDirectory& d){
  h.assign(NCuts_,NULL);
  for(int i=0; i<NCuts_; ++i)
    {
      string name = n + "_" + CutNames_[i];
      string title = t + " (After " + CutDescs_[i] + " Cut);"; 
      defineOneHisto(name, title, xtitle, nbins, min, max, units, h[i], d);
    }
}

void AnalyzerBase::defineOneHisto(const string & name, const string & title,
				  const string & xtitle, int nbins,
				  float min, float max,
				  const string & units, TH1F* & h, 
				  const TFileDirectory & d)
{
  float binWidth = (max-min)/nbins;
  string title2 = title + xtitle + ";Events";
  if(units.compare("NONE"))
    title2 += Form(" / %.0f ",binWidth) + units;
  h = d.make<TH1F>(name.c_str(),title2.c_str(),nbins,min,max);
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
void AnalyzerBase::endFile(std::vector<wprime::InputFile>::iterator fi,
                                 ofstream & out){
  tabulateFile(fi, fi->results);
  printFileSummary(fi, out);  
}

// e.g. print summmary of expected events for all samples
void AnalyzerBase::endAnalysis(ofstream & out){
  float N_SM = 0; 
  float N_DATA = 0; 
  std::vector<wprime::InputFile>::const_iterator it;

  out << endl;
  for(it = inputFiles_.begin(); it != inputFiles_.end(); ++it)
    { // loop over samples
      string sample = it->samplename;
      const wprime::FilterEff & stat = it->results.back();
      float N_evt = stat.Nsurv_evt_cut_w;
      out<< " "<< sample << ": " << N_evt
         << " evts (eff = " << 100.*stat.eff_abs
         << " +- " << 100.*stat.deff_abs
         << " %) " << endl;
      
      if     (wprimeUtil_->runningOnData()) 
        N_DATA += N_evt;
      else if(sample.find("wprime") == string::npos &&
              sample.find("Wprime") == string::npos &&
              sample.find("TC_WZ") == string::npos &&
              sample.find("RSZZ") == string::npos )
        N_SM += N_evt;
      
    } // loop over samples

  out << "  - - - - - - - - - - - - - " << endl;
  out << " Total # of SM events: " 
      << N_SM << endl;
  out << " Total # of DATA events: " 
      << N_DATA << endl;
}

bool AnalyzerBase::jsonContainsEvent (const vector<edm::LuminosityBlockRange>&jsonVec, const edm::EventBase &event)
{
  // if the jsonVec is empty, then no JSON file was provided so all
  // events should pass
  if (jsonVec.empty())
    {
      return true;
    }
  bool (* funcPtr) (edm::LuminosityBlockRange const &,
		    edm::LuminosityBlockID const &) = &edm::contains;
  edm::LuminosityBlockID lumiID (event.id().run(), 
				 event.id().luminosityBlock());
  vector< edm::LuminosityBlockRange >::const_iterator iter = 
    std::find_if (jsonVec.begin(), jsonVec.end(),
		  boost::bind(funcPtr, _1, lumiID) );
  return jsonVec.end() != iter;
}


void AnalyzerBase::run()
{
  int ievt_all=0;  int ievt_skipped = 0;
  unsigned i_sample = 1;
  vector<wprime::InputFile>::iterator it;

  for(it = inputFiles_.begin(); it != inputFiles_.end(); ++it, ++i_sample){
    int ievt=0;  
    cout << "\n Opening sample " << it->samplename 
         << " (" << it->description << ")... "<<std::flush;
    fwlite::ChainEvent ev(it->pathnames);
    it->Nact_evt = ev.size();
    cout<<" Done." << endl;
  
    cout << " Opened sample " << it->samplename << " with " << it->Nact_evt
         << " events (Input file #" << i_sample << " out of " << inputFiles_.size()
         << " samples) " << endl << endl;
    
    cout << std::fixed << std::setprecision(2);
    beginFile(it);

    unsigned runNumber = 0;
    unsigned lumiID = 0;
    nEvents_.assign(ctrNames_.size(), 0);
    for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){// loop over events
      edm::EventBase const & event = ev;

      if(countGenEvts_)
        updateEventcounts(ev, nEvents_, 
                          runNumber, lumiID, 
                          ctrNames_, false);
      
      // skip event if maximal number of events per input file is reached 
      if(maxEvents_>0 &&  ievt > maxEvents_) continue;
      
      // simple event counter
      if(reportAfter_!=0 ? (ievt>0 && ievt%reportAfter_==0) : false) 
        cout << " Processing event: " << ievt << " or " 
             << 100.*ievt/it->Nact_evt << "% of input file #" << i_sample
             << " (Total events processed: " << ievt_all 
             << ", non-certified/skipped: " << ievt_skipped << ") " << endl;
      
      if(useJSON_ && wprimeUtil_->runningOnData() &&
         !jsonContainsEvent (jsonVector, event))
      {
        ++ievt_skipped;
        continue;
      }
      else
        ++ievt_all;
      
      setEventWeight(event);
      eventLoop(event);
    } // loop over events
    if(countGenEvts_ && (int)nEvents_[0] != it->Nprod_evt) 
      cout<<"Weight Wrong: Found "<<nEvents_[0]<<" generated events and sample file lists "<<it->Nprod_evt<<endl;
    endFile(it, outLogFile_);
    
  } // loop over input files
  cout<<"Done with Input Samples\n";

  fs->cd(); 
  TH1F * h = new TH1F("lumi_ipb", "Integrated luminosity in pb^{-1}", 1, 0, 1);
  h->SetBinContent(1, wprimeUtil_->getLumi_ipb());
  h->SetBinContent(2, 1);//counter indicates number of files merged
  
  endAnalysis(outLogFile_);

}

void AnalyzerBase::setEventWeight(edm::EventBase const & event)
{
  wprimeUtil_->setHadronicMETcalculated(false);

  if(wprimeUtil_->runningOnData()){
    wprimeUtil_->setWeight(wprimeUtil_->getSampleWeight());
  }else{
    event.getByLabel(pileupLabel_, PupH_);
    //float PU_Weight = wprimeUtil_->getPUWeight1BX(*PupH_);//In time only
    //float PU_Weight = wprimeUtil_->getPUWeight3BX(*PupH_);//Average, in time and of of time
    float PU_Weight = wprimeUtil_->getPUWeight3D(*PupH_);//3D Matrix
    wprimeUtil_->setWeight(wprimeUtil_->getSampleWeight() * PU_Weight);
  }
}



