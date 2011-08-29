#include "UserCode/CMGWPrimeGroup/interface/AnalyzerBase.h"

using namespace std;

AnalyzerBase::AnalyzerBase(){}
AnalyzerBase::AnalyzerBase(const edm::ParameterSet & cfg, WPrimeUtil * wprimeUtil){
  //Cory: Put here everything you want EVERY CONSTRUCTOR TO DO

  wprimeUtil_ = wprimeUtil;
  assert(wprimeUtil_);

   Cuts_          = cfg.getParameter<vstring>("Cuts");
  NCuts_          = Cuts_.size();
  //FillCutFns();

  debugme = cfg.getParameter<bool>("debugme");

  eSelectorPset_ = cfg.getParameter<PSet>("electronSelectors");
  looseElectronType_ = cfg.getParameter<string>("LooseElectronType");
  tightElectronType_ = cfg.getParameter<string>("TightElectronType");
  looseElectron_ = ElectronSelector(eSelectorPset_, looseElectronType_);
  tightElectron_ = ElectronSelector(eSelectorPset_, tightElectronType_);
  electronResult_ = looseElectron_.getBitTemplate();
  if(debugme) cout<<"Using "<<looseElectronType_<<" for loose electrons and "
                  <<tightElectronType_<<" for tight electrons\n";

  mSelectorPset_ = cfg.getParameter<PSet>("muonSelectors");
  looseMuonType_ = cfg.getParameter<string>("LooseMuonType");
  tightMuonType_ = cfg.getParameter<string>("TightMuonType");
  looseMuon_ = MuonSelector(mSelectorPset_, looseMuonType_);
  tightMuon_ = MuonSelector(mSelectorPset_, tightMuonType_);
  muonResult_ = looseMuon_.getBitTemplate();
  if(debugme) cout<<"Using "<<looseMuonType_<<" for loose muons and "
                  <<tightMuonType_<<" for tight muons\n";

  jSelectorPset_ = cfg.getParameter<PSet>("jetSelectors");
  looseJetType_ = cfg.getParameter<string>("LooseJetType");
  looseJet_ = JetSelector(jSelectorPset_, looseJetType_);
  jetResult_ = looseJet_.getBitTemplate();
  if(debugme) cout<<"Using "<<looseJetType_<<" for jets\n";

  SetCandEvtFile(cfg.getParameter<string>("candEvtFile"));
  results_.assign(NCuts_,wprime::FilterEff());

  doPreselect_ = cfg.getParameter<bool>("preselect");

  electronsLabel_ = cfg.getParameter<edm::InputTag>("electrons");
  muonsLabel_ = cfg.getParameter<edm::InputTag>("muons");
  jetsLabel_ = cfg.getParameter<edm::InputTag>("jets");
  pfCandsLabel_ = cfg.getParameter<edm::InputTag>("particleFlow");
  metLabel_ = cfg.getParameter<edm::InputTag>("met");

  muonAlgo_ = cfg.getParameter<uint>("muonReconstructor");
  if(debugme) cout<<"Using muon algo "<<algo_desc_long[muonAlgo_]<<endl;
  useAdjustedMET_ = cfg.getParameter<bool>("useAdjustedMET");
  
  hltEventLabel_ = cfg.getParameter<edm::InputTag>("hltEventTag");
  pileupLabel_ = cfg.getParameter<edm::InputTag>("pileupTag");

  triggersToUse_ = cfg.getParameter<vstring>("triggersToUse");

  ////////////////////Default Cuts/////////////////
  minNLeptons_ = cfg.getUntrackedParameter<uint>("minNLeptons", 0);
  maxNLeptons_ = cfg.getUntrackedParameter<uint>("maxNLeptons", 999);
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
void AnalyzerBase::Tabulate_Me(const int& cut_index, const float& weight){
  if(debugme) cout<<"Tabulating results for cut_index = "
                  <<cut_index<<" = "<<Cuts_[cut_index]<<endl;

//increase the number of events passing the cuts
  hNumEvts->Fill(cut_index,weight);

  results_[cut_index].Nsurv_evt_cut_w += weight;
  results_[cut_index].Nsurv_evt_cut++;

  //fill the histograms
  Fill_Histos(cut_index,weight);
}//Tabulate_Me

//void AnalyzerBase::Fill_Histos(const int& index, const float& weight){}
//void AnalyzerBase::Declare_Histos(const TFileDirectory & dir){}

void AnalyzerBase::ResetCounters(){
  results_.assign(NCuts_,wprime::FilterEff());
}

/////Printers//////////////
void
AnalyzerBase::PrintEventFull(edm::EventBase const & event) const{
  WPrimeUtil::PrintEvent(event);
  WPrimeUtil::PrintPassingTriggers(triggerEvent_,triggersToUse_);
  PrintEventDetails();
  PrintEventLeptons();
}

void AnalyzerBase::PrintPassingEvent(edm::EventBase const & event){
  PrintEventToFile(event);
  WPrimeUtil::PrintEvent(event);
  PrintEventDetails();
}

void AnalyzerBase::PrintDebugEvent() const{
  WPrimeUtil::PrintPassingTriggers(triggerEvent_,triggersToUse_);
  PrintEventDetails();
  PrintEventLeptons();
  PrintLeptons();
}

void AnalyzerBase::PrintEventToFile(edm::EventBase const & event){
  outCandEvt_<<event.id().run()<<":"
             <<event.id().luminosityBlock()<<":"
             <<event.id().event()<<endl;
}

void AnalyzerBase::PrintEventDetails() const{
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
        <<" V Mass: "<<vCand_.eta()
        <<" V Mass: "<<vCand_.phi()
        <<endl;
  }
}

void
AnalyzerBase::PrintEventLeptons() const{
  cout<<"You shouldn't be seeing this! Implement your own.\n";
}

void
AnalyzerBase::PrintLeptons() const{
  PrintElectrons();
  PrintMuons();
  cout<<"----------------------\n";
}

void
AnalyzerBase::PrintElectrons() const{
  cout<<"----All Electrons: ( "<<electrons_.size()<<" )------\n";
  for(uint i=0; i<electrons_.size(); ++i){
    if(WPrimeUtil::Contains(electrons_[i], looseElectrons_)) cout<<"( Loose ";
    if(WPrimeUtil::Contains(electrons_[i], tightElectrons_)) cout<<"& Tight ";
    cout<<" )\n";
    PrintElectron(electrons_[i]);
  }
}

void
AnalyzerBase::PrintMuons() const{
  cout<<"----All Muons: ( "<<muons_.size()<<" )------\n";
  for(uint i=0; i<muons_.size(); ++i){
    if(WPrimeUtil::Contains(muons_[i], looseMuons_)) cout<<"( Loose ";
    if(WPrimeUtil::Contains(muons_[i], tightMuons_)) cout<<"& Tight ";
    cout<<" )\n";
    PrintMuon(muons_[i]);
  }
}

void
AnalyzerBase::PrintJets() const{
  cout<<"----All Jets: ( "<<jets_.size()<<" )------\n";
  for(uint i=0; i<jets_.size(); ++i){
    if(WPrimeUtil::Contains(jets_[i], looseJets_)) cout<<"(Loose ";
    cout<<" )\n";
    PrintJet(jets_[i]);
  }
}

void
AnalyzerBase::PrintElectron(const heep::Ele& elec) const{
  cout << setiosflags(ios::fixed) << setprecision(3);
  cout<<" Elec ScEt: "<<elec.et()<<endl; //ScEt
  if(!elec.isPatEle()){
    cout<<"Not a pat electron, why???\n";
    return;
  }
  PrintElectron(elec.patEle());
}

void
AnalyzerBase::PrintElectron(const pat::Electron& elec) const{
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
AnalyzerBase::PrintMuon(const TeVMuon& mu) const{
  cout << setiosflags(ios::fixed) << setprecision(3);
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
AnalyzerBase::PrintJet(const pat::Jet& jet) const{
  cout << setiosflags(ios::fixed) << setprecision(3);
  cout<<" Jet Pt: "  <<jet.pt()<<endl
      <<" Jet Eta: " <<jet.eta()<<endl
      <<" Jet Phi: " <<jet.phi()<<endl;
}

////////////////////////////
//////////Setters///////////
////////////////////////////
void AnalyzerBase::SetCandEvtFile(const std::string&  s){
  outCandEvt_.open(s.c_str());
  WPrimeUtil::CheckStream(outCandEvt_, s);
}

////////////////////
//////Cuts//////////
////////////////////

inline bool AnalyzerBase::PassNoCut() const{
  return true;
}

inline bool AnalyzerBase::PassTriggersCut() const{
//-----------------------------------------------------------
  return WPrimeUtil::PassTriggersCut(triggerEvent_,triggersToUse_);
}//--- PassTriggersCut

inline bool
AnalyzerBase::PassMinNLeptonsCut() const{
  return (looseElectrons_.size() + looseMuons_.size()) >= minNLeptons_;
}

inline bool
AnalyzerBase::PassMinNTightLeptonsCut() const{
  return (tightElectrons_.size() + tightMuons_.size()) >= minNTightLeptons_;
}
inline bool
AnalyzerBase::PassMaxNLeptonsCut() const{
  return (looseElectrons_.size() + looseMuons_.size()) <= maxNLeptons_;
}

inline bool
AnalyzerBase::PassMinNJetsCut() const{
  return looseJets_.size() >= minNJets_;
}

inline bool
AnalyzerBase::PassMinMETCut() const{
  return met_.et() > minMET_;
}

////////////////////////////////
/////////Check Z Properties/////
////////////////////////////////
inline bool
AnalyzerBase::PassValidZCut() const{
  return zCand_ && zCand_.mass()>0.;
}

inline bool
AnalyzerBase::PassZMassCut() const{
  return (zCand_.mass() > minZmass_) && (zCand_.mass() < maxZmass_);  
}

inline bool
AnalyzerBase::PassZptCut() const{
  return zCand_.pt() > minZpt_;
}

////////////////////////////////
///////Check Had. V Properties//
////////////////////////////////
inline bool
AnalyzerBase::PassValidVCut() const{
  return vCand_ && vCand_.mass()>0.;
}

inline bool AnalyzerBase::PassVMassCut() const{
  return (vCand_.mass() > minVmass_) && (vCand_.mass() < maxVmass_);
}

inline bool
AnalyzerBase::PassVptCut() const{
  return vCand_.pt() > minVpt_;
}

/////////Check W Properties/////
inline bool
AnalyzerBase::PassValidWCut() const{
  return wCand_ && wCand_.mass()>0.;
}

inline bool
AnalyzerBase::PassWtransMassCut() const{
  return wCand_.mt() > minWtransMass_;
}

inline bool
AnalyzerBase::PassWptCut() const{
  return wCand_.pt() > minWpt_;
}

///////Check DiBoson Properties//
/*
inline bool AnalyzerBase::PassValidVZCandCut() const{
  return vzCand_ && vzCand_.mass()>0.;
}
*/
//////////////////
//file stuff//////
//////////////////
void AnalyzerBase::beginFile(std::vector<wprime::InputFile>::const_iterator fi){
  TFileDirectory dir = wprimeUtil_->getFileService()->mkdir(fi->samplename); 
  Declare_Histos(dir);
  ResetCounters();
}

void
AnalyzerBase::DeclareHistoSet(const string& n, const string& t, const string& xtitle,
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

// operations to be done when closing input file 
// (e.g. print summary)
void AnalyzerBase::endFile(std::vector<wprime::InputFile>::const_iterator fi,
                                 ofstream & out){
  WPrimeUtil::tabulateSummary(results_);
  WPrimeUtil::printSummary(fi->samplename, fi->description, Cuts_, results_, out);  
}

void AnalyzerBase::endAnalysis(ofstream & out){
}

