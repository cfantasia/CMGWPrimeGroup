#include "UserCode/CMGWPrimeGroup/interface/WPrimeUtil.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

using std::cout; using std::cerr; using std::endl;
using std::string; using std::vector;


WPrimeUtil::WPrimeUtil(const char * out_filename, edm::InputTag genParticles, string cross_sections)
{
  fs = new fwlite::TFileService(out_filename);

  hRecoilPerp = 0;
  hRecoilParalvsVBPt = 0;
  histRecoilParal = NULL;
  lumi_ipb = -1;
  genParticles_ = genParticles;
  sample_cross_sections = cross_sections;
  setupZMETcorrection();

}
WPrimeUtil::~WPrimeUtil()
{
  delete fs; delete [] histRecoilParal;
}

void WPrimeUtil::setupZMETcorrection()
{
  string filename = "ZMET_data.root";
  // open the Z data file with info about recoil
  TFile* File = new TFile(filename.c_str(), "READONLY");
  if(!File || File->IsZombie())
    {
      cerr << " **** Oops! Missing file " << filename << endl;
      cerr << " Exiting... " << endl;
      abort();
    }

  TDirectoryFile * demo=(TDirectoryFile*)File->Get("pflow");
  hRecoilParalvsVBPt = (TH2D*)demo->Get("hMETParalvsVBPt");
  hRecoilParalvsVBPt->SetName("hRecoilParalvsVBPt");
  hRecoilPerp = (TH1D*)demo->Get("hMETPerp");
  hRecoilPerp->SetName("hRecoilPerp");
  
  
  setRecoilProjections();

}

void WPrimeUtil::setRecoilProjections()
{
  assert(hRecoilParalvsVBPt != 0);
  assert(histRecoilParal == 0);
  int N = hRecoilParalvsVBPt->GetXaxis()->GetNbins();
  histRecoilParal = new TH1D *[N];

  for(int bin_no = 1; bin_no <= N; ++bin_no)
    {
      // Get projection in the W pt bin
      histRecoilParal[bin_no-1] = new TH1D(*(hRecoilParalvsVBPt->ProjectionY("_pbinWpt_paral", bin_no, bin_no)));

      if(histRecoilParal[bin_no-1]->Integral()==0) {
	cout << " *** Couldn't correct for recoil for bin_no = " << bin_no
	     << endl;

      }
    }

}

// get input files (to be retrieved from samples_cross_sections.txt)
void WPrimeUtil::getInputFiles(std::vector<wprime::InputFile> & inputFiles)
{
  string txt_file = "UserCode/CMGWPrimeGroup/config/" + sample_cross_sections;
  ifstream in(txt_file.c_str());
  string new_line; wprime::InputFile * new_file = 0;
  while(getline(in, new_line))
  {
    // if DONE, we are done!
    if(new_line == "DONE")break;

    if(new_line.find("samplename = ") != string::npos)
      // new file found! create structure to put in info
      new_file = new wprime::InputFile();
    if(!new_line.empty())
      parseLine(new_line, new_file);
    else{ //Empty line means the end of a sample
      if(new_file){
        if(new_file->samplename.find("data") == string::npos)
          new_file->weight = lumi_ipb*(new_file->x_sect)
            /(new_file->Nprod_evt);
        else
          new_file->weight = 1;
        cout << " Weight to be applied on " << new_file->samplename
             << " sample = " << new_file->weight << endl << endl;
        // all info should now be logged in; check!
        new_file->checkFile();
      // if we made it here, everything looks good: 
      // add to vector of input files
        inputFiles.push_back(*new_file);
        // release memory
        delete new_file;
      }
    }
  }
}

void WPrimeUtil::SetLumiWeights(const string & MCFile, const string & DataFile,
                                const string & MCHist, const string & DataHist){
  LumiWeights_ = edm::LumiReWeighting(MCFile, DataFile, MCHist, DataHist);
}

float WPrimeUtil::getTotalWeight3BX(edm::EventBase const & event, const std::string& label){
  return runningOnData() ? getWeight() : getWeight()*getPUWeight3BX(event, label);
}

int WPrimeUtil::getPU1BX(const std::vector< PileupSummaryInfo > & PupInfo){
  std::vector<PileupSummaryInfo>::const_iterator PVI;                       
  for(PVI = PupInfo.begin(); PVI != PupInfo.end(); ++PVI) {               
    if(PVI->getBunchCrossing() == 0){//Only care about in time PU for now 
      return PVI->getPU_NumInteractions();           
    }
  }//Looping over different Bunch Crossings
  return -1;
}

float WPrimeUtil::getPUWeight1BX(const std::vector< PileupSummaryInfo > & PupInfo){
  return LumiWeights_.weight(getPU1BX(PupInfo));
}

/////////////////////////////////////////////////
////Average Over In-time & out-of-time PU////////
/////////////////////////////////////////////////
float WPrimeUtil::getPUWeight3BX(edm::EventBase const & event, const std::string& label){
  std::vector< PileupSummaryInfo > PupInfo = getProduct<std::vector< PileupSummaryInfo > >(event, label);   
  return getPUWeight3BX(PupInfo);
}

float WPrimeUtil::getPU3BX(const std::vector< PileupSummaryInfo > & PupInfo){
  std::vector<PileupSummaryInfo>::const_iterator PVI;
  float sum_nvtx = 0;
  for(PVI = PupInfo.begin(); PVI != PupInfo.end(); ++PVI) {
    float npv = PVI->getPU_NumInteractions();
    sum_nvtx += float(npv);
  }
  return sum_nvtx/3.;//+1, 0, -1 BX
}

float WPrimeUtil::getPUWeight3BX(const std::vector< PileupSummaryInfo > & PupInfo){
  return LumiWeights_.weight3BX( getPU3BX(PupInfo) );
}

void WPrimeUtil::CheckStream(const ofstream& stream, const std::string & s){
  if(!stream) { 
    std::cout << "Cannot open file " << s << std::endl; 
    abort();
  } 
}

// Calculate efficiencies
//------------------------------------------------------------------------
void WPrimeUtil::getEff(float & eff, float & deff,float Num,float Denom)
{
  //------------------------------------------------------------------------
  if(Denom){
    eff = Num/Denom;
    deff = TMath::Sqrt(eff * (1-eff)/Denom);
  }else 
    eff = deff = 0.;
}//---------------getEff


// used for the parsing of samples_cross_sections.txt
void WPrimeUtil::parseLine(const string & new_line, wprime::InputFile * in_file)
{
  size_t i = 0;

  if(top_level_dir.empty())
  {
    string search = "Top Level Dir = ";
    i = new_line.find(search);
    if(i != string::npos)
    {
      top_level_dir = new_line.substr(search.length(), new_line.length() - search.length());
      cout << "\n Reading PAT-tuples from directory " << top_level_dir << endl;
      return;
    }
  }
  // this is where we extract the integrated luminosity
  if(lumi_ipb < 0)
	{
    i = new_line.find("integrated lumi = ");
    if(i != string::npos)
    {
      unsigned int begin = new_line.find("=") + 2;
      unsigned int end = new_line.find("pb^-1");
      string integ_lumi = new_line.substr(begin, end-begin);
      lumi_ipb = atof(integ_lumi.c_str());
      
      cout << "\n Will apply weights to MC samples to get distributions for "
           << lumi_ipb << " pb^-1" << endl << endl;
      if(lumi_ipb <= 0)
      {
        cout << " *** Oops, read in \"" << integ_lumi 
             << "\" which I failed to parse! Please seek help... " 
             << endl;
        abort();
      }
      return;
    }
  }

  i = new_line.find("samplename = ");
  if(i != string::npos)
    {
      in_file->samplename = new_line.substr(13, new_line.length() - 13);
      return;
    }

  i = new_line.find("description = ");
  if(i != string::npos)
    {
      in_file->description = new_line.substr(14, new_line.length() - 14);
      
      return;
    }

 i = new_line.find("subdir = ");
  if(i != string::npos)
    {
      in_file->subdir = new_line.substr(9, new_line.length() - 9);

      return;
    }


  i = new_line.find("pathname = ");
  if(i != string::npos)
    {
      string input = new_line.substr(11, new_line.length() - 11);
      size_t j = new_line.find(".txt");
      if(j != string::npos){
        cout<<"Using input file "<<input<<endl;
        ifstream infile;
        input = "UserCode/CMGWPrimeGroup/config/" + input;
        infile.open (input.c_str(), ifstream::in);
        
        while (infile.good()){
          string fname;
          infile>>fname;
          if(fname != ""){
            cout<<" filename: "<<fname.c_str()<<endl;
	    in_file->pathnames.push_back(top_level_dir + in_file->subdir+ "/" + fname);
          }
        }
        infile.close();
      }else{
        string pathname = top_level_dir + in_file->subdir + "/" + input;
	in_file->pathnames.push_back(pathname);
      }
      cout << " Input file: " << in_file->samplename;
      cout << " (" << in_file->description << ") " << endl;
      return;
    }

  i = new_line.find("x-section = ");
  if(i != string::npos)
    {
      in_file->x_sect = atof(new_line.substr(12, new_line.length() - 12).c_str());
      return;
    }

  i = new_line.find("Nprod_evt = ");
  if(i != string::npos)
    {
      in_file->Nprod_evt = atoi(new_line.substr(12, new_line.length() - 12).c_str());
      return;
    }
  
}

//////////////////
//Print Functions/
//////////////////
void WPrimeUtil::PrintEvent(edm::EventBase const & event){
  cout<<"run #: "<<event.id().run()
      <<" lumi: "<<event.id().luminosityBlock()
      <<" eventID: "<<event.id().event()<<endl;
}

//////////////////
//Trigger Fns/////
//////////////////
bool WPrimeUtil::PassTriggersCut(edm::EventBase const & event, std::string label, const std::vector<std::string>& triggerNames){
  pat::TriggerEvent triggerEvent = getProduct<pat::TriggerEvent>(event, label);
  return PassTriggersCut(triggerEvent,triggerNames);
}

bool WPrimeUtil::PassTriggersCut(const pat::TriggerEvent & triggerEvent,const std::vector<std::string>& triggerNames){
  const pat::TriggerPathRefVector acceptedPaths = triggerEvent.acceptedPaths();
  //cout<<"Using "<<acceptedPaths.size()<<" accepted paths from HLT"<<endl;
  for (size_t i = 0; i < acceptedPaths.size(); i++){
    if(FoundAndPassed(triggerEvent, acceptedPaths[i], triggerNames)) return true;
  }//acceptedPaths loop
  return false;
}

void
WPrimeUtil::PrintPassingTriggers(const pat::TriggerEvent & triggerEvent,const std::vector<std::string>& triggerNames){
  const pat::TriggerPathRefVector acceptedPaths = triggerEvent.acceptedPaths();
//  const pat::TriggerAlgorithmRefVector algoBits = triggerEvent.physAlgorithms();
//  for (size_t i = 0; i < algoBits.size(); i++){
//    cout<<" L1 algo: "<<algoBits[i]->name()<<" with prescale "<<algoBits[i]->prescale()<<endl;
//  }
  for (size_t i = 0; i < acceptedPaths.size(); i++){
    if(FoundAndPassed(triggerEvent, acceptedPaths[i], triggerNames))
      cout<<"Passed path: "<<acceptedPaths[i]->name()<<endl;
  }
}

inline bool
WPrimeUtil::FoundAndPassed(const pat::TriggerEvent & triggerEvent,const pat::TriggerPathRef path,const std::vector<std::string>& triggerNames){
  return FindTrigger(path, triggerNames) && Passed(triggerEvent,path);
}

inline bool
WPrimeUtil::Passed(const pat::TriggerEvent & triggerEvent, const pat::TriggerPathRef path){
  return (path->wasAccept() && path->prescale() == 1 && (1 || MaxL1Prescale(triggerEvent,path)==1) );
}

inline unsigned
WPrimeUtil::L1Prescale(const pat::TriggerEvent & triggerEvent, const pat::TriggerPathRef path){
  assert(path->l1Seeds().size() == 1);
  return triggerEvent.algorithm(path->l1Seeds()[0].second)->prescale();
}

unsigned
WPrimeUtil::MaxL1Prescale(const pat::TriggerEvent & triggerEvent,const pat::TriggerPathRef path){
  unsigned ps = 1;
  pat::L1SeedCollection l1algos = path->l1Seeds();
  //cout<<" Looking at "<<path->name()<<" and its "<<l1algos.size()<<" l1 seeds"<<endl;
  for (size_t j = 0; j < l1algos.size(); j++){
    //cout<<" and its seed "<<l1algos[j].second<<" and decision "<<l1algos[j].first<<endl;
    const pat::TriggerAlgorithm * l1algo = triggerEvent.algorithm(l1algos[j].second);
    ////////////////////
    //Looks like the format is l1algos[0] && l1algos[1] ....
    //but l1algos can be 'L1seedA OR L1seedB', with the word 'OR' actually in the string!
    ////////////////////
    //if(l1algo) cout<<" and name "<<l1algo->name()<<" and decision "<<l1algo->decision()<<" prescale "<<l1algo->prescale()<<endl;
    //else       cout<<" and name "<<l1algos[j].second<<" failed!!!!"<<endl;
    if(!l1algo) return 0;
    ps = std::max(ps, l1algo->prescale());
  }
  return ps;
}

inline bool
WPrimeUtil::FindTrigger(const pat::TriggerPathRef path, const std::vector<std::string>& triggerNames){
  for (size_t j = 0; j < triggerNames.size(); j++){
    if(SameTrigger(path->name(), triggerNames[j])) return true;
  }
  return false;
}

////////////////////
//Trigger Matching//
////////////////////
/*
bool 
HadronicVWAnalyzer::PassTriggerMatch(const pat::Electron & p, const float cut, const vstring& triggers) const{
  for (size_t i=0; i < triggers.size(); ++i){
    if (p.triggerObjectMatchesByPath(triggers[i], true, false).size() > 0){
      const pat::TriggerObjectStandAlone * trigRef = p.triggerObjectMatchByPath(triggers[i], true, false);
      if(trigRef->et() > cut) return true;
    }
  }
  return false;
}

bool 
HadronicVWAnalyzer::PassTriggerMatch(const TeVMuon & p, const float cut, const vstring& triggers) const{
  for(uint i=0; i<p.triggerObjectMatches().size(); ++i){
    vector<string> names = p.triggerObjectMatches()[i].pathNames(true, false);
    for(uint j=0; j<names.size(); ++j){
      for (size_t k=0; k < triggers.size(); ++k){
        if(WPrimeUtil::SameTrigger(names[j], triggers[k])){
          if (p.triggerObjectMatchesByPath(names[j], true, false).size() > 0){
            if(p.triggerObjectMatchByPath(names[j], true, false)->pt() > cut) return true;
          }
        }
      }
    }
  }
  return false;
}
*/


// get hadronic MET component (that needs to be corrected 
// if applyMETCorrection=true)from Z data; this will be done according to hadronic 
// activity from Z->mumu reconstructed events
TVector2 WPrimeUtil::getHadronicMET(edm::EventBase const & event)
{
  if(hadronicMETcalculated_)
    return hadronicMETcached;

  assert(event.isRealData() == false);

  int W_index = -1;
  event.getByLabel(genParticles_, genParticles);
  for(size_t i = 0; i != genParticles->size(); ++i) {
    const reco::GenParticle & W_p = (*genParticles)[i];
    if (W_p.pdgId() == 24 && W_p.status() == 3)
      {
	W_index = i;
	break;
      }
  } // loop over genParticles

  assert(W_index >= 0);
  const reco::GenParticle & W_p = (*genParticles)[W_index];

  // correct hadronic MET by taking into account recoil for given W pt
  
  // Find bin corresponding to W pt
  int binWpt = (hRecoilParalvsVBPt->GetXaxis())->FindBin(W_p.pt());

  // protect against outliers: EITHER use the last bin OR return MET w/o correction
  if(binWpt > hRecoilParalvsVBPt->GetXaxis()->GetNbins())
    binWpt = hRecoilParalvsVBPt->GetXaxis()->GetNbins();

  // Shoot random MET (parallel) from the Z histograms
  double dataSampledMETParal= histRecoilParal[binWpt-1]->GetRandom();
  // Shoot perpendicular component
  double dataSampledMETPerp = hRecoilPerp->GetRandom();
  
  // Rotate back from system in which the boson is in the x axis
  double cosW = W_p.px()/W_p.pt();
  double sinW = W_p.py()/W_p.pt();
  
  double dataSampledMEx = 
    cosW*dataSampledMETParal - sinW*dataSampledMETPerp;
  double dataSampledMEy = 
    sinW*dataSampledMETParal + cosW*dataSampledMETPerp;
  
  hadronicMETcached.Set(dataSampledMEx, dataSampledMEy);
  setHadronicMETCalculated(true);
  return hadronicMETcached;    

}

// true if current file under processing contains "data" in name/description
void WPrimeUtil::setRunningOnData()
{
  runningOnData_ = (getSampleName().find("data") != string::npos);
}

//Check if Run/Evt is in Debug list
bool WPrimeUtil::DebugEvent(edm::EventBase const& event) const
{
  const edm::EventID & evtToCheck = event.id();
  for(uint i=0; i<vEventsToDebug_.size(); ++i){
    if(evtToCheck == vEventsToDebug_[i])
      return true;
  }
  return false;
}

/////////////////////
void WPrimeUtil::convertElectrons(const vector<pat::Electron>& patElectrons,
                                  ElectronV & electrons){
  electrons.clear();
  for (size_t i = 0; i < patElectrons.size(); i++) {
    electrons.push_back(heep::Ele(patElectrons[i]));   
  }
}

void WPrimeUtil::convertMuons(const vector<pat::Muon>& patMuons, const uint& muAlgo,
                              MuonV & muons){
  muons.clear();
  for (size_t i = 0; i < patMuons.size(); i++) {
    muons.push_back(TeVMuon(patMuons[i], muAlgo));   
  }
}

void WPrimeUtil::getElectrons(edm::EventBase const & event, const edm::InputTag& label, ElectronV & electrons){
  const vector<pat::Electron> patElectrons = getProduct<vector<pat::Electron> >(event, label);
  convertElectrons(patElectrons, electrons);
}

void WPrimeUtil::getMuons(edm::EventBase const & event, const edm::InputTag& label, const uint& muAlgo, MuonV & muons){
  const vector<pat::Muon> patMuons = getProduct<vector<pat::Muon> >(event, label);
  convertMuons(patMuons, muAlgo, muons);
}

inline void WPrimeUtil::getPFCands(edm::EventBase const & event, const edm::InputTag& label, vector<pat::PFParticle> & pfCands){
  edm::Handle<vector<pat::PFParticle> > handle;
  event.getByLabel(label, handle);
  if(handle.isValid()) pfCands = *handle.product();
  else std::cerr<<" Didn't find pfCands in event, skipping\n";
}

inline void WPrimeUtil::getMET(edm::EventBase const & event, const edm::InputTag& label, pat::MET & met){
  met = getProduct<METV>(event, label)[0];
}

void
WPrimeUtil::AdjustMET(edm::EventBase const & event, 
                      const ElectronV & electrons, const MuonV & muons,
                      const edm::InputTag& pfLabel,  pat::MET & met){
  vector<pat::PFParticle> pfCands;
  getPFCands(event, pfLabel, pfCands);
  AdjustMET(electrons, pfCands, met);
  AdjustMET(muons, pfCands, met);
}

/*
void WPrimeUtil::getElectronsMET(edm::EventBase const & event,
                                 const edm::InputTag& eLabel, ElectronV & electrons,
                                 const edm::InputTag& metLabel, const bool & adjMET, pat::MET & met,
                                 const edm::InputTag& pfLabel){
  getElectrons(event, eLabel, electrons);
  getMET(event, metLabel, met);
  if(adjMET) AdjustMET(event, pfLabel, met);
}

void WPrimeUtil::getMuonsMET(edm::EventBase const & event,
                             const edm::InputTag& muLabel, const uint& muAlgo, MuonV & muons,
                             const edm::InputTag& metLabel, const bool & adjMET, pat::MET & met,
                             const edm::InputTag& pfLabel){
  getMuons(event, muLabel, muAlgo, muons);
  getMET(event, metLabel, met);
  if(adjMET) AdjustMET(event, pfLabel, met);
}

void WPrimeUtil::getLeptonsMET(edm::EventBase const & event, 
                               const edm::InputTag& eLabel, ElectronV & electrons,
                               const edm::InputTag& muLabel, const uint& muAlgo, MuonV & muons,
                               const edm::InputTag& metLabel, const bool & adjMET, pat::MET & met,
                               const edm::InputTag& pfLabel){
  getElectrons(event, eLabel, electrons);
  getMuons(event, muLabel, muAlgo, muons);
  getMET(event, metLabel, met);
  if(adjMET) AdjustMET(event, pfLabel, met);
}
*/

void WPrimeUtil::getElectronsMET(edm::EventBase const & event,
                                 const vector<pat::Electron>& patElectrons, ElectronV & electrons,
                                 const edm::InputTag& metLabel, const bool & adjMET, pat::MET & met,
                                 const edm::InputTag& pfLabel){
  convertElectrons(patElectrons, electrons);
  getMET(event, metLabel, met);
  if(adjMET) AdjustMET(event, electrons, pfLabel, met);
}

void WPrimeUtil::getMuonsMET(edm::EventBase const & event,
                             const vector<pat::Muon    >& patMuons, const uint& muAlgo, MuonV & muons,
                             const edm::InputTag& metLabel, const bool & adjMET, pat::MET & met,
                             const edm::InputTag& pfLabel){
  convertMuons(patMuons, muAlgo, muons);
  getMET(event, metLabel, met);
  if(adjMET) AdjustMET(event, muons, pfLabel, met);
}

void WPrimeUtil::getLeptonsMET(edm::EventBase const & event, 
                               const vector<pat::Electron>& patElectrons, ElectronV & electrons,
                               const vector<pat::Muon    >& patMuons, const uint& muAlgo, MuonV & muons,
                               const edm::InputTag& metLabel, const bool & adjMET, pat::MET & met,
                               const edm::InputTag& pfLabel){
  convertElectrons(patElectrons, electrons);
  convertMuons(patMuons, muAlgo, muons);
  getMET(event, metLabel, met);
  if(adjMET) AdjustMET(event, electrons, muons, pfLabel, met);
}

void WPrimeUtil::getElectronsMET(const PatElectronVH & patElectronsH, ElectronV & electrons,
                                 const METVH & metH, const bool & adjMET, pat::MET & met,
                                 const PFCandidateVH & pfCandidatesH){
  convertElectrons(*patElectronsH.product(), electrons);
  met = (*metH.product())[0];
  if(adjMET) AdjustMET(electrons, *pfCandidatesH.product(), met);
}

void WPrimeUtil::getMuonsMET(const PatMuonVH & patMuonsH, const uint& muAlgo, MuonV & muons,
                             const METVH & metH, const bool & adjMET, pat::MET & met,
                             const PFCandidateVH & pfCandidatesH){
  convertMuons(*patMuonsH.product(), muAlgo, muons);
  met = (*metH.product())[0];
  if(adjMET) AdjustMET(muons, *pfCandidatesH.product(), met);
}

void WPrimeUtil::getLeptonsMET(const PatElectronVH & patElectronsH, ElectronV & electrons,
                               const PatMuonVH & patMuonsH, const uint& muAlgo, MuonV & muons,
                               const METVH & metH, const bool & adjMET, pat::MET & met,
                               const PFCandidateVH & pfCandidatesH){
  convertElectrons(*patElectronsH.product(), electrons);
  convertMuons(*patMuonsH.product(), muAlgo, muons);
  met = (*metH.product())[0];
  if(adjMET){
    AdjustMET(electrons, *pfCandidatesH.product(), met);
    AdjustMET(muons, *pfCandidatesH.product(), met);
  }
}

////////////////////////
///Tabulating Functions
////////////////////////

//Writing results to a txt file
//--------------------------------------------------------------------------
void WPrimeUtil::tabulateSummary(wprime::EffV results){
  for(uint i = 0; i < results.size(); ++i){
    //calculate efficiencies
    int idx = std::max((int)i-1,0);
    float num =   results[i].Nsurv_evt_cut_w;
    float denom = results[idx].Nsurv_evt_cut_w;
    WPrimeUtil::getEff(results[i].eff, results[i].deff, num, denom);
    WPrimeUtil::getEff(results[i].eff_abs, results[i].deff_abs, num, results[0].Nsurv_evt_cut_w);
  } // loop over different cuts
}//tabulateSummary

void WPrimeUtil::printSummary(const string& dir, const string& description, const vstring & Cuts, const wprime::EffV results, ofstream & out){ 
  out << setiosflags(std::ios::fixed) << std::setprecision(2);
  out<<"$$$$$$$$$$$$$$$$$$$$$$$ Type of sample: "<<dir<<" ( "<<description<<" )"<<endl;
//  out << " xsec*lumi expected events = " << Nthe_evt << endl;
//  out << " # of evt passing preselection = " << Nexp_evt << " per "<<Form("%.0f", wprimeUtil_->getLumi_ipb())<<" inv pb"<<endl;
  
  for(uint i = 0; i < Cuts.size(); ++i){
    out<<std::right<<"Cut " << std::setw(2) << i << "("
       <<std::left<< std::setw(15) << Cuts[i]
       <<std::right << "): " <<"Evts = " << std::setw(10) << results[i].Nsurv_evt_cut_w
       <<" (" << std::right << std::setw(7) << results[i].Nsurv_evt_cut << ")";
    
    out << std::setw(9) <<"\tRel eff: "<<std::setw(6)<<results[i].eff*100
        << " +/- " << std::setw(6)<<results[i].deff*100 << "%"
        << std::setw(9) <<"\tAbs eff: "<<std::setw(6)<<results[i].eff_abs *100
        << " +/- " << std::setw(6)<<results[i].deff_abs*100 << "%"
        << endl;
  } // loop over different cuts
}//printSummary

