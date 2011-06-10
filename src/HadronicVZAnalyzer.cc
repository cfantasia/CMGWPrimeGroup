#include "UserCode/CMGWPrimeGroup/interface/HadronicVZAnalyzer.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember)) 

using namespace std;
HadronicVZAnalyzer::HadronicVZAnalyzer(){}
HadronicVZAnalyzer::HadronicVZAnalyzer(const edm::ParameterSet & cfg, WPrimeUtil * wprimeUtil){

  wprimeUtil_ = wprimeUtil;
  assert(wprimeUtil_);

  SetLogFile(cfg.getParameter<string>("logFile"));
  SetCandEvtFile(cfg.getParameter<string>("candEvtFile"));

  intOptions_["report"] = cfg.getParameter<uint>("reportAfter");
  intOptions_["verbose"] = cfg.getParameter<bool>("debugme");
  doPreselect_ = cfg.getParameter<bool>("preselect");
  intOptions_["events"] = 0;

  debugme = cfg.getParameter<bool>("debugme");

  muonsLabel_ = cfg.getParameter<string>("muons");
  jetsLabel_ = cfg.getParameter<string>("jets");

  muonAlgo_ = cfg.getParameter<uint>("muonAlgo");
  
  hltEventLabel_ = cfg.getParameter<string>("hltEventTag");
  pileupLabel_ = cfg.getParameter<string>("pileupTag");

  triggersToUse_          = cfg.getParameter<vstring>("triggersToUse");

  PDGMUON = 13;
  PDGELEC = 11;
  PDGW = 24;
  PDGZ = 23;
  PDGZPRIME = 33;
  PDGWPRIME = 34;

  PI    = TMath::Pi();
  TWOPI = TMath::TwoPi();
  NOCUT = 9e9;
 
// +++++++++++++++++++Event characteristics
  

// +++++++++++++++++++General Cut values
  maxNumZs = cfg.getParameter<uint>("maxNumZs");
  minNLeptons = cfg.getParameter<uint>("minNLeptons");
  minNJets = cfg.getParameter<uint>("minNumJets");
  maxNJets = cfg.getParameter<uint>("maxNumJets");
  minLeadPt = cfg.getParameter<double>("minLeadPt");
  maxAngleBetweenJets = cfg.getParameter<double>("maxAngleBetweenJets");

// +++++++++++++++++++Z Cuts
  minZpt = cfg.getParameter<double>("minZpt");
  minZmass = cfg.getParameter<double>("minZmass");
  maxZmass = cfg.getParameter<double>("maxZmass");

// +++++++++++++++++++Hadronic Boson Cuts
  minHadVpt = cfg.getParameter<double>("minHadVPt");
  minHadVmass = cfg.getParameter<double>("minHadVmass");
  maxHadVmass = cfg.getParameter<double>("maxHadVmass");

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

// +++++++++++++++++++Jet General Cuts
  minJetPt = cfg.getParameter<double>("minJetPt");
  maxJetEta = cfg.getParameter<double>("maxJetEta");
  maxJetNHF = cfg.getParameter<double>("maxJetNHF");
  maxJetNEF = cfg.getParameter<double>("maxJetNEF");
  minJetCHF = cfg.getParameter<double>("minJetCHF");
  maxJetCEF = cfg.getParameter<double>("maxJetCEF");
  minJetnumConst = cfg.getParameter<uint>("minJetnumConst");
  minJetcMult = cfg.getParameter<uint>("minJetcMult");
}

HadronicVZAnalyzer::~HadronicVZAnalyzer(){
  outCandEvt.close(); 
  outLogFile.close(); 
}

/// Declare Histograms
//--------------------------------------------------------------
void HadronicVZAnalyzer::Declare_Histos(TFileDirectory & dir)
{
  // For now we have only one histo for testing, may extend later.
  printf("Declare histos\n");
  h_HadVZMass = dir.make<TH1F>("h_HadVZMass","h_HadVZMass",200,0.0,2000.0);

}//Declare_Histos

//Fill Histograms
//-----------------------------------------------------------
void HadronicVZAnalyzer::Fill_Histos(int index, float weight)
{
  // For now we have only one histo, so no need for this function.
  printf("Filling Histos\n");
}//Fill_Histos

void 
HadronicVZAnalyzer::eventLoop(edm::EventBase const & event){

  // Preselection - skip events that don't look promising
  if (doPreselect_){
    if(debugme) cout<<"Testing Preselection...\n";
    // We could setup some preselection here. To be implemented.
  }

  // Get leptons
  const vector<pat::Muon    > patMuons     = getProduct<vector<pat::Muon    > >(event, muonsLabel_);
  // Get jets
  const vector<pat::Jet     > patJets      = getProduct<vector<pat::Jet     > >(event, jetsLabel_);

  printf("    Contains: %i muon(s)\n",
	 (int)patMuons.size());
  printf("    Contains: %i jet(s)\n",
	 (int)patJets.size());

  // Make vectors of leptons passing various criteria
  // Loop over muons, and see if they pass the TeVMuon criteria  
  for (size_t i = 0; i < patMuons.size(); i++) {
    muons_.push_back(TeVMuon(patMuons[i],muonAlgo_));   
    if (PassMuonCut(&muons_[i]))
      looseMuons_.push_back(muons_[i]);
  }

  // Loop over jets, and see if they pass the jet criteria
  for (size_t i = 0; i < patJets.size(); ++i) {
    if (PassJetCut(&patJets[i]))
      jets_.push_back(patJets[i]);
  }
  
  bool passedNLeptons = PassNLeptonsCut();
  bool passedNJets = PassNJetsCut();
  
  if(!(passedNLeptons && passedNJets))
    return;
  
  // Make a Z candidate out of the muons. 
  zCand = getZCands(looseMuons_).front();
  // Make a W candidate out of the jets.
  wCand = getWCand(jets_);

  bool validZ    = PassValidZCut();
  bool validHadV = PassValidHadVCut();
  bool goodZ     = PassZMassCut() && PassZptCut();
  bool goodHadV  = PassHadVMassCut() && PassHadVptCut();
  
  if(validZ && validHadV && goodZ && goodHadV) {
    reco::CompositeCandidate hadVZ;
    hadVZ.addDaughter(zCand);
    hadVZ.addDaughter(wCand);
    AddFourMomenta addP4;
    addP4.set(hadVZ);
    h_HadVZMass->Fill(hadVZ.mass());
  }
}

/////////////////Accessors///////////////////////

/////////////////Modifies///////////////////////
void HadronicVZAnalyzer::CheckStream(ofstream& stream, string s){
  if(!stream) { 
    cout << "Cannot open file " << s << endl; 
    abort();
  } 
}

void HadronicVZAnalyzer::SetCandEvtFile(string s){
  outCandEvt.open(s.c_str());
  CheckStream(outCandEvt, s);
}

void HadronicVZAnalyzer::SetLogFile(string s){
  outLogFile.open(s.c_str());      
  CheckStream(outLogFile, s); 
}

/////////////////Cuts///////////////////////

//Always true
bool HadronicVZAnalyzer::PassNoCut()
{
  return true;
}

//Trigger requirements
//-----------------------------------------------------------
bool HadronicVZAnalyzer::PassTriggersCut()
{
//-----------------------------------------------------------
  if(debugme) cout<<"Trigger requirements"<<endl;
  // To be implemented
  return true;
}//--- PassTriggersCut()

bool
HadronicVZAnalyzer::PassNLeptonsCut(){
  return (looseMuons_.size()) > minNLeptons;
}

bool
HadronicVZAnalyzer::PassNJetsCut(){
  return ((jets_.size() > minNJets) && (jets_.size() < maxNJets));
}

bool
HadronicVZAnalyzer::PassValidHadVCut(){
  return wCand && wCand.mass()>0.;
}

bool
HadronicVZAnalyzer::PassValidZCut(){
  return zCand && zCand.mass()>0.;
}

bool
HadronicVZAnalyzer::PassValidHadVZCandCut(){
  // To be implemented
  return true;
}

bool
HadronicVZAnalyzer::PassNumberOfZsCut(){
  return numZs < maxNumZs;
}

bool
HadronicVZAnalyzer::PassLeadingLeptonPtCut(){
  return LeadPt > minLeadPt;
}

////////////////////////////////
/////////Check Z Properties/////
////////////////////////////////
bool
HadronicVZAnalyzer::PassZMassCut(){
  return (zCand.mass() > minZmass) && (zCand.mass() < maxZmass);  
}

bool
HadronicVZAnalyzer::PassZptCut(){
  return zCand.pt() > minZpt;
}
////////////////////////////////
///////Check Had. V Properties//
////////////////////////////////

bool HadronicVZAnalyzer::PassHadVMassCut(){
  return (wCand.mass() > minHadVmass) && (wCand.mass() < maxHadVmass);
}

bool
HadronicVZAnalyzer::PassHadVptCut(){
  return wCand.pt() > minHadVpt;
}

////////////////////////////////
//////Check Muon Properties/////
////////////////////////////////
/// Big function to check all cuts.
bool HadronicVZAnalyzer::PassMuonCut(const TeVMuon* mu){
  for(uint i=0; i<MuonCutFns_.size(); ++i){
    bool result = CALL_MEMBER_FN(*this,MuonCutFns_[i])(mu);
    if(result==false) return false;
  }
  return true;
}

bool HadronicVZAnalyzer::PassMuonLoosePtCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Loose Pt Cut"<<endl;
  return (mu->pt() > minMuonLoosePt);
}

bool HadronicVZAnalyzer::PassMuonTightPtCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Tight Pt Cut"<<endl;
  return (mu->pt() > minMuonTightPt);
}//--- PassMuonPtCut

bool HadronicVZAnalyzer::PassMuonGlobalCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Global Cut"<<endl;
  return (mu->isGlobalMuon()); 
}//--- PassMuonGlobalCut

bool HadronicVZAnalyzer::PassMuonNpixhitCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon NpixhitCut"<<endl;
  return (mu->globalTrack()->hitPattern().numberOfValidPixelHits() > minMuonNPixHit);
}//--- PassMuonNpixhitCut

bool HadronicVZAnalyzer::PassMuonNtrkhitCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon NtrkhitCut"<<endl;
  return (mu->globalTrack()->hitPattern().numberOfValidTrackerHits() > minMuonNTrkHit);
}//--- PassMuonNtrkhitCut

bool HadronicVZAnalyzer::PassMuonNormChi2Cut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Chi2 Cut"<<endl;
  return (mu->globalTrack()->normalizedChi2() < maxMuonNormChi2);
}//--- PassMuonChi2Cut

bool HadronicVZAnalyzer::PassMuonHitsUsedCut(const TeVMuon* mu){
  //Num Valid Muon Hits
  if(debugme) cout<<"Check Muon Hits Used Cut"<<endl;
  return (mu->globalTrack()->hitPattern().numberOfValidMuonHits() > minMuonHitsUsed);
}//--- PassMuonHits Used Cut

bool HadronicVZAnalyzer::PassMuonStationsCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Stations Cut"<<endl;
  return (mu->numberOfMatches() > minMuonStations);
}//--- PassMuonStationsCut

bool HadronicVZAnalyzer::PassMuonEtaCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Eta Cut"<<endl;
  return (fabs(mu->eta()) < maxMuonEta);
}//--- PassMuonEta Cut

/*
bool HadronicVZAnalyzer::PassMuonCombRelIsoCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon CombRelIso Cut"<<endl;
  return (Calc_MuonRelIso(mu) < maxWmunuCombRelIso);
}//--- PassMuonCombRelIsoCut
*/

bool HadronicVZAnalyzer::PassMuonDxyCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Dxy Cut"<<endl;
  return (fabs(mu->userFloat("d0")) < maxMuonDxy);
}//--- PassMuonDxyCut

////////////////////////////////
/////Check jet Properties //////
////////////////////////////////
bool HadronicVZAnalyzer::PassJetCut(const pat::Jet* jet){
  bool jetKinStatus = PassJetPtCut(jet) && PassJetEtaCut(jet);
  bool jetIDStatus = PassJetIDCut(jet);
  return jetKinStatus && jetIDStatus;
}

bool HadronicVZAnalyzer::PassJetPtCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet Pt Cut"<<endl;
  return (jet->pt() > minJetPt);
}

bool HadronicVZAnalyzer::PassJetEtaCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet Eta Cut"<<endl;
  return (fabs(jet->eta()) < maxJetEta);
}

bool HadronicVZAnalyzer::PassJetNHFCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet NHF Cut"<<endl;
  return (jet->neutralHadronEnergyFraction() < maxJetNHF);
}

bool HadronicVZAnalyzer::PassJetNEFCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet NEF Cut"<<endl;
  return (jet->neutralEmEnergyFraction() < maxJetNEF);
}

bool HadronicVZAnalyzer::PassJetNConstCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet NConst Cut"<<endl;
  return (jet->numberOfDaughters() > minJetnumConst);
}

bool HadronicVZAnalyzer::PassJetCHFCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet CHF Cut"<<endl;
  return (jet->chargedHadronEnergyFraction() > minJetCHF);
}

bool HadronicVZAnalyzer::PassJetCMultCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet CMult Cut"<<endl;
  return ((unsigned int)jet->chargedMultiplicity() > minJetcMult);
}

bool HadronicVZAnalyzer::PassJetCEFCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet CEF Cut"<<endl;
  return (jet->chargedEmEnergyFraction() < maxJetCEF);
}

bool HadronicVZAnalyzer::PassJetIDCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Pass JetID Cut"<<endl;
  bool neutralStatus = 
    PassJetNHFCut(jet) &&
    PassJetNEFCut(jet) &&
    PassJetNConstCut(jet);
  bool chargedStatus = (fabs(jet->eta()) > 2.4) ||
    (PassJetCHFCut(jet) &&
     PassJetNHFCut(jet) &&
     PassJetCMultCut(jet)
     );
  return (neutralStatus && chargedStatus);
}

////////////////////////////////
/////Check TeV Properties/////
////////////////////////////////

/*float
HadronicVZAnalyzer::Calc_MuonRelIso(const TeVMuon* mu){
  return (mu->isolationR03().emEt + mu->isolationR03().hadEt + mu->isolationR03().sumPt)
    / mu->pt();
}
*/
///////////////Utilities//////////////////
/*
void
HadronicVZAnalyzer::ClearEvtVariables(){
  muons_.clear();
  looseMuons_.clear();
  tightMuons_.clear();
  zCand = ZCandidate();
  wCand = WCandidate();
}
*/
void 
HadronicVZAnalyzer::reportProgress(int eventNum) {
  if (eventNum % intOptions_["report"] == 0) {
    printf("\rWe've processed %i events so far...", eventNum);
    cout.flush();
    printf("\n");
  }
  printf("Event number: %i", ++eventNum);
}

/// Print to screen (like printf), but only if --verbose option is on
void HadronicVZAnalyzer::verbose(const char *string, ...)
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
void HadronicVZAnalyzer::beginFile(std::vector<wprime::InputFile>::const_iterator fi){
  TFileDirectory dir = wprimeUtil_->getFileService()->mkdir(fi->samplename); 
  Declare_Histos(dir);
  //ResetCounters();
}

// operations to be done when closing input file 
// (e.g. print summary)
void HadronicVZAnalyzer::endFile(std::vector<wprime::InputFile>::const_iterator fi,
                         ofstream & out){
  //ScaleHistos();//Already scaled
  //printSummary(fi->samplename);  
  //deleteHistos();
  //listOfHists.clear();
}

void HadronicVZAnalyzer::endAnalysis(ofstream & out){
}
