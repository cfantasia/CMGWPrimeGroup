#include "UserCode/CMGWPrimeGroup/interface/AnalyzerBase.h"

using namespace std;

AnalyzerBase::AnalyzerBase(){}
AnalyzerBase::AnalyzerBase(const edm::ParameterSet & cfg, int fileToRun){
  //Put here everything you want EVERY CONSTRUCTOR TO DO

  reportAfter_ = cfg.getParameter<int>("reportAfter");
  maxEvents_   = cfg.getParameter<int>("maxEvents");
  useJSON_   = cfg.getParameter<bool>("useJSON") ;
  doPreselect_ = cfg.getParameter<bool>("preselect");
  debug_ = cfg.getParameter<bool>("debug");

  genLabel_ = cfg.getParameter<edm::InputTag>("genParticles" );
  pfCandsLabel_ = cfg.getParameter<edm::InputTag>("particleFlow");
  pileupLabel_ = cfg.getParameter<edm::InputTag>("pileupTag" );
  hltEventLabel_ = cfg.getParameter<edm::InputTag>("hltEventTag");
  doRecoilCorrectionForW_ = cfg.getParameter<bool>("doRecoilCorrectionForW");

  // file with samples & cross-sections
  string sample_cross_sections = cfg.getParameter<string>("sample_cross_sections");
  edm::ParameterSet const& inputs = cfg.getParameter<edm::ParameterSet>("inputs");
  if ( inputs.exists("lumisToProcess") ) 
    {
      vector<edm::LuminosityBlockRange> const & lumisTemp =
	inputs.getUntrackedParameter<vector<edm::LuminosityBlockRange> > ("lumisToProcess");
      jsonVector.resize( lumisTemp.size() );
      copy( lumisTemp.begin(), lumisTemp.end(), jsonVector.begin() );
    }

  string MCPUDistFile   = cfg.getParameter<string>("MCPUDistFile" );
  string MCPUDistHist   = cfg.getParameter<string>("MCPUDistHist" );
  string DataPUDistFile = cfg.getParameter<string>("DataPUDistFile" );
  string DataPUDistHist = cfg.getParameter<string>("DataPUDistHist" );
  float  puScale        = cfg.getParameter<double>("puScale");
  
  std::vector<edm::EventID> vEventsToDebug = cfg.getParameter<std::vector<edm::EventID> >("vEventsToDebug");
  
  wprimeUtil_ = new WPrimeUtil(genLabel_, pfCandsLabel_, sample_cross_sections);
  wprimeUtil_->setLumiWeights(MCPUDistFile, DataPUDistFile, MCPUDistHist, DataPUDistHist, puScale);
  vEventsToDebug_ = vEventsToDebug;
  //wprimeUtil_->setEventsToDebug(vEventsToDebug);
  assert(wprimeUtil_);

  wprimeUtil_->getInputFiles(inputFiles_, fileToRun);
  if(fileToRun != -1){
    if(fileToRun < (int)inputFiles_.size()){
      inputFiles_.assign(1,inputFiles_[fileToRun]);
    }else{
      cerr<<"You asked for sample "<<fileToRun
          <<" but only "<<inputFiles_.size()
          <<" are listed!\nHere is the list\n";
      for(unsigned i=0; i<inputFiles_.size(); ++i) 
        cerr<<" Sample "<<i<<": "<<inputFiles_[i].samplename
            <<" (" << inputFiles_[i].description << ")\n";
      inputFiles_.clear();
      abort();
    }
  }

  //////Output Files
  string outputFile  = cfg.getParameter<string  >("outputFile" );
  string logFile  = cfg.getParameter<string  >("logFile" );
  string candEvtFile = cfg.getParameter<string>("candEvtFile");

  if(fileToRun != -1){
    string tail = Form("_Sample%i.",fileToRun);
    outputFile.replace(outputFile.find_last_of("."), 1, tail); 
    logFile.replace(logFile.find_last_of("."), 1, tail); 
    candEvtFile.replace(candEvtFile.find_last_of("."), 1, tail); 
  }

  fs = new fwlite::TFileService(outputFile);
  outLogFile_.open(logFile.c_str());
  WPrimeUtil::CheckStream(outLogFile_, logFile);
  outCandEvtFile_.open(candEvtFile.c_str());
  WPrimeUtil::CheckStream(outCandEvtFile_, candEvtFile);

  setNumSignalSamples();

  //////////////////////

   CutNames_ = cfg.getParameter<vstring>("Cuts");
   CutDescs_ = CutNames_;
  NCuts_     = CutNames_.size();

  //////////////

  electronsLabel_ = cfg.getParameter<edm::InputTag>("electrons");
  muonsLabel_ = cfg.getParameter<edm::InputTag>("muons");
  jetsLabel_ = cfg.getParameter<edm::InputTag>("jets");
  metLabel_ = cfg.getParameter<edm::InputTag>("met");
  vertexLabel_ = cfg.getParameter<edm::InputTag>("vertexTag");

  muReconstructor_ = cfg.getParameter<int>("muonReconstructor");
  if(debug_) cout<<"Using muon algo "<<algo_desc_long[muReconstructor_]<<endl;
  assert(muReconstructor_ < Num_MuTeVtrkAlgos);

  useAdjustedMET_ = cfg.getParameter<bool>("useAdjustedMET");
  
  triggersToUse_ = cfg.getParameter<vstring>("triggersToUse");

  //Muon Algos
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
void AnalyzerBase::fillCuts(const map<string,fnCut >& mFnPtrs){
  CutFns_.resize(NCuts_);
  for(int i=0; i<NCuts_; ++i){
    if(mFnPtrs.find(CutNames_[i]) == mFnPtrs.end()){
      cout<<"Didn't find cut named "<<CutNames_[i]<<endl;
      abort();
    }
    CutFns_[i] = mFnPtrs.find(CutNames_[i])->second;
  } 
}

//Tabulate results after the cut has been passed
void AnalyzerBase::tabulateEvent(const int& cut_index, const float& weight){
  if(debug_) cout<<"Tabulating results for cut_index = "
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
void AnalyzerBase::tabulateFile(std::vector<wprime::InputFile>::iterator fi){
  for(uint i = 0; i < fi->results.size(); ++i){
    //calculate efficiencies
    float rel_denom = i==0 ? fi->Nact_evt : fi->results[i-1].Nsurv_evt_cut;
    //Calc Relative Eff
    WPrimeUtil::getEff(fi->results[i].eff, fi->results[i].deff, 
                       fi->results[i].Nsurv_evt_cut, rel_denom);
    //Calc Absolute Eff
    WPrimeUtil::getEff(fi->results[i].eff_abs, fi->results[i].deff_abs, 
                       fi->results[i].Nsurv_evt_cut, fi->Nprod_evt);
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
}

const reco::Vertex& 
AnalyzerBase::findPV(const vector<reco::Vertex>& vtxs) const{
  //loop over vtx and find pv
  int pvidx = 0;
  float maxsumpt = 0;
  for(uint i=0; i<vtxs.size(); ++i){
    reco::Vertex::trackRef_iterator it;
    float sumpt = 0;
    for(it = vtxs[i].tracks_begin(); it!= vtxs[i].tracks_end(); ++it){
      sumpt += (*it)->pt();
    }
    if(sumpt > maxsumpt){
      maxsumpt = sumpt;
      pvidx = i;
    }
  }
  return vtxs[pvidx];
}

const float maxVtxDr = 0.01;
/*
const reco::Vertex& 
AnalyzerBase::findPV(const vector<reco::Vertex>& vtxs, const reco::Candidate& p) const{
  for(uint i=0; i<vtxs.size(); ++i){
    reco::Vertex::trackRef_iterator it;
    for(it = vtxs[i].tracks_begin(); it!= vtxs[i].tracks_end(); ++it){
      if(deltaR(**it, p) < maxVtxDr) return vtxs[i];
    }
  }
  cout<<"Didnt find a vertex!"<<endl;
  return vtxs[0];
}
*/

const reco::Vertex& 
AnalyzerBase::findPV(const vector<reco::Vertex>& vtxs, const reco::Candidate& p) const{
  float minDz = 9e9;
  int pvidx = 0;
  for(uint i=0; i<vtxs.size(); ++i){
    if(fabs(vtxs[i].z() - p.vz()) < minDz){
      minDz = fabs(vtxs[i].z() - p.vz());
      pvidx = i;
    }
  }
  cout<<"Best dist is "<<minDz<<endl;
  return vtxs[pvidx];
}

bool
AnalyzerBase::sameVertex(const reco::Vertex& vtx, const reco::Candidate& p) const{
  reco::Vertex::trackRef_iterator it;
  for(it = vtx.tracks_begin(); it!= vtx.tracks_end(); ++it){
    //if(deltaR(**it, p) < maxVtxDr) cout<<"track pt, p pt, dR"<<(*it)->pt()<<" "<<p.pt()<<" "<<deltaR(**it, p)<<endl;
    if(deltaR(**it, p) < maxVtxDr) return true;
  }
  return false;
}


/////printers//////////////
void
AnalyzerBase::printEventFull(edm::EventBase const & event) const{
  WPrimeUtil::printEvent(event, cout);
  WPrimeUtil::printPassingTriggers(*triggerEventH_,triggersToUse_);
  printEventDetails();
}

void AnalyzerBase::printPassingEvent(edm::EventBase const & event){
  WPrimeUtil::printEvent(event, outCandEvtFile_);
  WPrimeUtil::printEvent(event, cout);
  printEventDetails();
}

void AnalyzerBase::printDebugEvent() const{
  WPrimeUtil::printPassingTriggers(*triggerEventH_,triggersToUse_);
  printEventDetails();
}

void AnalyzerBase::printEventDetails() const{
}


void
AnalyzerBase::print(const heep::Ele& elec) const{
  cout << setiosflags(ios::fixed) << setprecision(2);
  cout<<" Elec ScEt: "<<elec.et()<<endl; //ScEt
  if(!elec.isPatEle()){
    cout<<"Not a pat electron, why???\n";
    return;
  }
  print(elec.patEle());
}

void
AnalyzerBase::print(const pat::Electron& elec) const{
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

void AnalyzerBase::print(const TeVMuon& mu) const{
  cout << setiosflags(ios::fixed) << setprecision(2);
  
  cout << " Muon eta = " << mu.eta() << "  phi = " << mu.phi() << endl;

  typedef std::vector<unsigned>::const_iterator It;
  for(It it = muon_reconstructors.begin(); it != muon_reconstructors.end(); ++it)
    mu.printPtInfo(*it);
  
  mu.printTrackerInfo();

}

void
AnalyzerBase::print(const pat::Jet& jet) const{
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
  if (debug_) cout<<"In pass Cuts\n";
  
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
  return WPrimeUtil::passTriggersCut(*triggerEventH_,triggersToUse_);
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

inline bool
AnalyzerBase::passWptCut(const WCandidate& w, const float& cut) const{
  return w.pt() > cut;
}

/////////Check XWLeptonic Properties/////
inline bool
AnalyzerBase::passValidXWCut(const XWLeptonic& xw) const{
  return xw && xw().mass()>0.;
}

inline bool
AnalyzerBase::passXWMassCut(const XWLeptonic& xw, const NuAlgos & algo, const float& min, const float& max) const{
  return xw(algo).mass() > min && xw(algo).mass() < max;
}

inline bool
AnalyzerBase::passXWptCut(const XWLeptonic& xw, const float& cut) const{
  return xw().pt() > cut;
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

  TFileDirectory dir = fs->mkdir(fi->samplename, fi->description); 

  //Save # of events after each cut
  string title = Form("Expected # of Events / %.0f pb^{-1}",  wprimeUtil_->getLumi_ipb());
  title = title + ";;" + title;
  hNumEvts = NULL; hNumEvts = dir.make<TH1F>("hNumEvts",title.c_str(),NCuts_,0,NCuts_);
  for(int i=0; i<NCuts_; ++i) hNumEvts->GetXaxis()->SetBinLabel(i+1,CutNames_[i].c_str());

  //Save Important File Info//////////
  title = "Info for Sample " + fi->samplename + " ( " + fi->description + " )";
  TH1F* hFileInfo = NULL; 
  hFileInfo = dir.make<TH1F>("hFileInfo",title.c_str(),8,0,8);
  hFileInfo->GetXaxis()->SetBinLabel  (1, "#intL dt (pb^{-1})");
  hFileInfo->SetBinContent(1, wprimeUtil_->getLumi_ipb());
  hFileInfo->GetXaxis()->SetBinLabel  (2, "#sigma (pb)");
  hFileInfo->SetBinContent(2, fi->x_sect);
  hFileInfo->GetXaxis()->SetBinLabel  (3, "Number of Events Produced");
  hFileInfo->SetBinContent(3, fi->Nprod_evt);
  hFileInfo->GetXaxis()->SetBinLabel  (4, "Number of Events in root Files");
  hFileInfo->SetBinContent(4, fi->Nact_evt);
  hFileInfo->GetXaxis()->SetBinLabel  (5, "Sample Weight");
  hFileInfo->SetBinContent(5, fi->weight);
  hFileInfo->GetXaxis()->SetBinLabel  (6, "Number of root files");
  hFileInfo->SetBinContent(6, fi->pathnames.size());
  hFileInfo->GetXaxis()->SetBinLabel  (7, "Signal Mass");
  hFileInfo->SetBinContent(7, fi->signalMass);
  hFileInfo->GetXaxis()->SetBinLabel  (8, "Number of Files Merged");
  hFileInfo->SetBinContent(8, 1);
  
  defineHistos(dir);
  if(wprimeUtil_->isSignalSample())
    {
      ++signalSample_index;
      defineResolutionHistos(dir, fi->signalMass);
    }
  resetCounters();
  fi->results.assign(NCuts_,wprime::FilterEff());
}

void
AnalyzerBase::defineHistoSet(const string& n, const string& t, 
                             const string& xtitle, 
                             int nxbins, float xmin, float xmax, 
                             const string& ytitle,
                             int nybins, float ymin, float ymax, 
                             vector<TH2F*>& h, const TFileDirectory& d){
  h.assign(NCuts_,NULL);
  for(int i=0; i<NCuts_; ++i)
    {
      string name = n + "_" + CutNames_[i];
      string title = t + " (After " + CutDescs_[i] + " Cut)"; 
      defineOneHisto(name, title, xtitle, nxbins, xmin, xmax, 
                     ytitle, nybins, ymin, ymax, h[i], d);
    }
}

void
AnalyzerBase::defineHistoSet(const string& n, const string& t, 
			     const string& xtitle,
			     int nbins, float xmin, float xmax, 
			     const string& units,
			     vector<TH1F*>& h, const TFileDirectory& d){
  h.assign(NCuts_,NULL);
  for(int i=0; i<NCuts_; ++i)
    {
      string name = n + "_" + CutNames_[i];
      string title = t + " (After " + CutDescs_[i] + " Cut)"; 
      defineOneHisto(name, title, xtitle, nbins, xmin, xmax, units, h[i], d);
    }
}

// mass format expected in <x>.<y> TeV, e.g. "1.2", corresponding to 1.2 TeV
// channel could be "e", "mu", "ee", "mumu", etc
void AnalyzerBase::createResolutionHist(const TFileDirectory & d, float Mass,
					const string & channel, 
					TH1F* & put_here)
{
  string name = string("Res_") + channel;
  string title = string("Resolution function for M = ") + Form("%.1f", Mass) + 
    string(" TeV (channel = ") + channel + string(")");
  string xtitle = "GeV/c^{2}";
  int nbins = 100;
  // need a smart way to adjust the resolution histogram range based on Mass
  float xmin = -800; float xmax = 800;
  defineOneHisto(name, title, xtitle, nbins, xmin, xmax, "GeV", put_here, d);
}

// mass format expected in <x>.<y> TeV, e.g. "1.2", corresponding to 1.2 TeV
// channel could be "e", "mu", "ee", "mumu", etc
void AnalyzerBase::createGenMtHist(const TFileDirectory & d, float Mass,
					const string & channel, 
					TH1F* & put_here)
{
  string name = string("GenMt_") + channel;
  string title = string("GEN Mt function for M = ") + Form("%.1f", Mass) + 
    string(" TeV (channel = ") + channel + string(")");
  string xtitle = "GeV/c^{2}";
  int nbins = 100;
  // need a smart way to adjust the resolution histogram range based on Mass
  float xmin = Mass*1000 - 1300; float xmax = Mass*1000 + 300;
  defineOneHisto(name, title, xtitle, nbins, xmin, xmax, "GeV", put_here, d);
}

void AnalyzerBase::defineOneHisto(const string & name, const string & title,
				  const string & xtitle, int nxbins,
				  float xmin, float xmax,
				  const string & ytitle, int nybins,
				  float ymin, float ymax,
          TH2F* & h, 
				  const TFileDirectory & d)
{
  string title2 = title + ";" + xtitle + ";" + ytitle;
  h = d.make<TH2F>(name.c_str(),title2.c_str(),nxbins,xmin,xmax,nybins,ymin,ymax);
}

void AnalyzerBase::defineOneHisto(const string & name, const string & title,
				  const string & xtitle, int nbins,
				  float xmin, float xmax,
				  const string & units, TH1F* & h, 
				  const TFileDirectory & d)
{
  float binWidth = (xmax-xmin)/nbins;
  string title2 = title + ";" + xtitle + ";Events";
  if(units.compare("NONE"))
    title2 += Form(" / %.0f ",binWidth) + units;
  h = d.make<TH1F>(name.c_str(),title2.c_str(),nbins,xmin,xmax);
}

// operations to be done when closing input file 
// (e.g. print summary)
void AnalyzerBase::endFile(std::vector<wprime::InputFile>::iterator fi,
                                 ofstream & out){
  tabulateFile(fi);
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

// return index of current signal sample (-1 if this is not a signal sample)
int AnalyzerBase::getSignalSampleIndex()
{
  if(!wprimeUtil_->isSignalSample())
    return -1;
  
  return signalSample_index;
}

void AnalyzerBase::setNumSignalSamples()
{
  NumSigSamples = 0;

  vector<wprime::InputFile>::const_iterator it;
  for(it = inputFiles_.begin(); it != inputFiles_.end(); ++it)
      if(it->isSignal())
	++NumSigSamples;
	
  // initialize index pointing to signal samples
  signalSample_index = -1;

}

void AnalyzerBase::run()
{
  int ievt_all=0;  int ievt_skipped = 0;
  unsigned i_sample = 1;
  float reportPercent = reportAfter_<0 ? reportAfter_/100. : 0;
  vector<wprime::InputFile>::iterator it;

  //Loop Over input files
  for(it = inputFiles_.begin(); it != inputFiles_.end(); ++it, ++i_sample){
    int ievt=0;//Evts checked in this input file  
    cout << "\n Opening sample " << it->samplename 
         << " (" << it->description << ")... "<<std::flush;
    fwlite::ChainEvent ev(it->pathnames);
    it->Nact_evt = ev.size();
    cout<<" Done." << endl;
  
    cout << " Opened sample " << it->samplename << " with " << it->Nact_evt
         << " events (Input file #" << i_sample << " out of " << inputFiles_.size()
         << " samples) " << endl << endl;   
    cout << std::fixed << std::setprecision(2);

    beginFile(it);//Set up for input file
    assert(it->Nact_evt <= it->Nprod_evt);

    if(reportPercent) reportAfter_ = fabs(it->Nact_evt * reportPercent);

    if(vEventsToDebug_.size()){//If events are in this vector, only process those
      cout<<" EventsToDebug has "<<vEventsToDebug_.size()<<" event(s), so only those events will be processed\n";
      debug_ = true;
      for(uint i=0; i<vEventsToDebug_.size(); ++i){
        if(!ev.to(vEventsToDebug_[i])){//Didn't find debug event
          cout<<" Didn't find debug event ( "<<vEventsToDebug_[i]<<" ) in this sample, skipping this event\n";
          continue;
        }
        edm::EventBase const & event = ev;
        eventLoop(event); //Analyze debug event
      }
    }else{
      for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){// loop over events
        edm::EventBase const & event = ev;

        // skip event if maximal number of events per input file is reached 
        if(maxEvents_>0 &&  ievt > maxEvents_) break;
      
        // simple event counter
        if(reportAfter_!=0 ? (ievt>0 && ievt%reportAfter_==0) : false) 
          cout << " Processing event: " << ievt << " or " 
               << 100.*ievt/it->Nact_evt << "% of input file #" << i_sample
               << " (Total events processed: " << ievt_all 
               << ", non-certified/skipped: " << ievt_skipped << ") " << endl;
      
        //skip event if not in JSON
        if(useJSON_ && wprimeUtil_->runningOnData() &&
           !jsonContainsEvent (jsonVector, event))
        {
          ++ievt_skipped;
          continue;
        }
        else
          ++ievt_all;
      
        setEventWeight(event); //PU Reweighting
        eventLoop(event); //Analyze each event
      } // loop over events
    }
    endFile(it, outLogFile_);//Finish up with this input sample
    
  } // loop over input files
  cout<<"Done with Input Samples\n";

  fs->cd(); 
  TH1F * h = new TH1F("lumi_ipb", "Integrated luminosity in pb^{-1}", 1, 0, 1);
  h->SetBinContent(1, wprimeUtil_->getLumi_ipb());
  h->SetBinContent(2, 1);//counter indicates number of files merged
  
  endAnalysis(outLogFile_);//Finish up with analsis

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



