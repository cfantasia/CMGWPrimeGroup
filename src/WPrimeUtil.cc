#include "UserCode/CMGWPrimeGroup/interface/WPrimeUtil.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

using std::cout; using std::cerr; using std::endl;
using std::string; using std::vector;


WPrimeUtil::WPrimeUtil(edm::InputTag genLabel, 
		       edm::InputTag pfLabel, string cross_sections)
{

  hRecoilPerp = 0;
  hRecoilParalvsVBPt = 0;
  histRecoilParal = NULL;
  lumi_ipb = -1;
  genLabel_ = genLabel;
  pfLabel_ = pfLabel;
  sample_cross_sections = cross_sections;
  setupZMETcorrection();

}
WPrimeUtil::~WPrimeUtil()
{
  int N = hRecoilParalvsVBPt->GetXaxis()->GetNbins();
  for(int bin_no = 0; bin_no < N; ++bin_no)
    delete histRecoilParal[bin_no];
  delete [] histRecoilParal;
  delete hRecoilPerp;
  delete hRecoilParalvsVBPt;
}

void WPrimeUtil::setupZMETcorrection()
{
  string filename = "ZMET_data.root";
  // open the Z data file with info about recoil
  TFile* File = TFile::Open(filename.c_str(), "READONLY");
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
      // get projection in the W pt bin
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
    //if starts with #, ignore the line
    if(new_line.find("#") == 0) continue;

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
             << " sample = " << new_file->weight << endl;
        if(new_file->isSignal())
          cout << " This is a signal sample with mass = " 
               << new_file->signalMass << " TeV " << endl;

        //This section will split the list of input files into a number of
        //smaller samples for faster processing in parallel.  
        //Caveat: Splitting is not recommended if you do not run in parallel
        //since the directory (sample) name is not being changed.
        //Default is not to split (splitInto = 1)
        vstring pathnames = new_file->pathnames;
        size_t splitInto = new_file->splitInto;
        size_t nPerFile = ceil((float) pathnames.size() / splitInto);
        if(splitInto > 1) 
          cout<<"Trying to split file with "<<pathnames.size()<<" files into "
              <<splitInto<<" parts, with "<<nPerFile<<" files per part"<<endl;
        for(size_t i=0; i<splitInto; ++i){
          size_t first = i*nPerFile;
          size_t last  = std::min(first+nPerFile, pathnames.size());
          if(first >= pathnames.size()) break;
          new_file->pathnames.assign(pathnames.begin()+first,
                                     pathnames.begin()+last);

          // all info should now be logged in; check!
          new_file->checkFile();
          // if we made it here, everything looks good: 
          // add to vector of input files
          inputFiles.push_back(*new_file);
        }
        //////////
        cout << endl;
        // release memory
        delete new_file;
      }
    }
  }
}

void WPrimeUtil::setLumiWeights(const string & MCFile, const string & DataFile,
                                const string & MCHist, const string & DataHist){
  LumiWeights3D_ = edm::Lumi3DReWeighting(MCFile, DataFile, MCHist, DataHist);
  LumiWeights3D_.weight3D_init( 1.0);
  LumiWeights_ = edm::LumiReWeighting(MCFile, DataFile, MCHist, DataHist);
  //LumiWeights_.weight3D_init("UserCode/CMGWPrimeGroup/root_macros/Weight3D.root");
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
float WPrimeUtil::getPUWeight3BX(const edm::EventBase & event, const std::string& label){
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

////////////
float WPrimeUtil::getPUWeight3D(const std::vector< PileupSummaryInfo > & PupInfo){
  std::vector<PileupSummaryInfo>::const_iterator PVI;
  int nm1 = -1; int n0 = -1; int np1 = -1;
  for(PVI = PupInfo.begin(); PVI != PupInfo.end(); ++PVI) {
    int BX = PVI->getBunchCrossing();
    
    if     (BX == -1) nm1 = PVI->getPU_NumInteractions();
    else if(BX ==  0)  n0 = PVI->getPU_NumInteractions();
    else if(BX ==  1) np1 = PVI->getPU_NumInteractions();
  }
  return LumiWeights3D_.weight3D( nm1,n0,np1);
}
////////////

void WPrimeUtil::CheckStream(const ofstream& stream, const std::string & s){
  if(!stream) { 
    std::cout << "Cannot open file " << s << std::endl; 
    abort();
  } 
}

// calculate efficiencies
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

  i = new_line.find("signalMass = ");
  if(i != string::npos)
    {
      in_file->signalMass = atof(new_line.substr(13, new_line.length() - 13).c_str());
      return;
    }


  i = new_line.find("Nprod_evt = ");
  if(i != string::npos)
    {
      in_file->Nprod_evt = atoi(new_line.substr(12, new_line.length() - 12).c_str());
      return;
    }
  
  i = new_line.find("splitInto = ");
  if(i != string::npos)
    {
      in_file->splitInto = atoi(new_line.substr(12, new_line.length() - 12).c_str());
      return;
    }
}

//////////////////
//print Functions/
//////////////////
void WPrimeUtil::printEvent(const edm::EventBase & event){
  cout<<"run #: "<<event.id().run()
      <<" lumi: "<<event.id().luminosityBlock()
      <<" eventID: "<<event.id().event()<<endl;
}

//////////////////
//Trigger Fns/////
//////////////////
bool WPrimeUtil::passTriggersCut(const edm::EventBase & event, std::string label, const std::vector<std::string>& triggerNames){
  pat::TriggerEvent triggerEvent = getProduct<pat::TriggerEvent>(event, label);
  return passTriggersCut(triggerEvent,triggerNames);
}

bool WPrimeUtil::passTriggersCut(const pat::TriggerEvent & triggerEvent,const std::vector<std::string>& triggerNames){
  const pat::TriggerPathCollection* paths = triggerEvent.paths();
  for (size_t i = 0; i < paths->size(); i++){
    if(FoundAndpassed(triggerEvent, paths->at(i), triggerNames)) return true;
  }//Paths loop
  return false;
}

void
WPrimeUtil::printPassingTriggers(const pat::TriggerEvent & triggerEvent,const std::vector<std::string>& triggerNames){
  const pat::TriggerPathCollection* paths = triggerEvent.paths();
  for (size_t i = 0; i < paths->size(); i++){
    if(FoundAndpassed(triggerEvent, paths->at(i), triggerNames))
      cout<<"passed path: "<<paths->at(i).name()<<endl;
  }
}

inline bool
WPrimeUtil::FoundAndpassed(const pat::TriggerEvent & triggerEvent,const pat::TriggerPath& path,const std::vector<std::string>& triggerNames){
  return passed(triggerEvent,path) && FindTrigger(path, triggerNames);
}

const bool ignoreL1prescale = true;
inline bool
WPrimeUtil::passed(const pat::TriggerEvent & triggerEvent, const pat::TriggerPath& path){
  return (path.wasAccept() && path.prescale() == 1 && (ignoreL1prescale || MaxL1Prescale(triggerEvent,path)==1) );
}

inline unsigned
WPrimeUtil::L1Prescale(const pat::TriggerEvent & triggerEvent, const pat::TriggerPath& path){
  assert(path.l1Seeds().size() == 1);
  return triggerEvent.algorithm(path.l1Seeds()[0].second)->prescale();
}

unsigned
WPrimeUtil::MaxL1Prescale(const pat::TriggerEvent & triggerEvent,const pat::TriggerPath& path){
  unsigned ps = 1;
  pat::L1SeedCollection l1algos = path.l1Seeds();
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
WPrimeUtil::FindTrigger(const pat::TriggerPath& path, const std::vector<std::string>& triggerNames){
  for (size_t j = 0; j < triggerNames.size(); j++){
    if(SameTrigger(path.name(), triggerNames[j])) return true;
  }
  return false;
}


// return pointer to gen-particle with pdgId and mother pdgId_mother
const reco::Candidate * WPrimeUtil::getGenParticle(edm::EventBase const & event,int pdgId, int pdgId_mother)
{
  event.getByLabel(genLabel_, genParticles);
  for(size_t i = 0; i != genParticles->size(); ++i) {
    const reco::GenParticle & p = (*genParticles)[i];
    if(TMath::Abs(p.pdgId()) != pdgId)continue;
    if(p.status() != 3)continue;
    const reco::Candidate * mother  = findMother(&p);
    if(!mother)continue;
    if(TMath::Abs(mother->pdgId()) != pdgId_mother)continue;
    return &p;
  } // loop over genParticles
  return (const reco::Candidate *) 0;
}



// get hadronic MET component (that needs to be corrected 
// if applyMETCorrection=true)from Z data; this will be done according to hadronic 
// activity from Z->mumu reconstructed events
TVector2 WPrimeUtil::getHadronicMET(const edm::EventBase & event)
{
  if(hadronicMETcalculated_)
    return hadronicMETcached;

  assert(event.isRealData() == false);

  int W_index = -1;
  event.getByLabel(genLabel_, genParticles);
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
  setHadronicMETcalculated(true);
  return hadronicMETcached;    

}

//Check if Run/Evt is in Debug list
bool WPrimeUtil::DebugEvent(const edm::EventBase & event) const
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

void WPrimeUtil::getElectrons(const edm::EventBase & event, const edm::InputTag& label, ElectronV & electrons){
  const vector<pat::Electron>& patElectrons = getProduct<vector<pat::Electron> >(event, label);
  convertElectrons(patElectrons, electrons);
}

void WPrimeUtil::getMuons(const edm::EventBase & event, const edm::InputTag& label, const uint& muAlgo, MuonV & muons){
  const vector<pat::Muon>& patMuons = getProduct<vector<pat::Muon> >(event, label);
  convertMuons(patMuons, muAlgo, muons);
}

bool WPrimeUtil::warningShown_ = false;

inline void WPrimeUtil::getPFCands(const edm::EventBase & event, const edm::InputTag& label, vector<pat::PFParticle> & pfCands){
  edm::Handle<vector<pat::PFParticle> > handle;
  event.getByLabel(label, handle);
  if(handle.isValid()) pfCands = *handle.product();
  else 
    {
      if(!warningShown_)
	{
	  std::cerr << " WARNING!!! Didn't find pfCands in event, skipping corrections to pfMET\n";
	  warningShown_ = true;
	}
    }
}

void
WPrimeUtil::AdjustMET(const edm::EventBase & event, 
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

void WPrimeUtil::getElectronsMET(const edm::EventBase & event,
                                 const vector<pat::Electron>& patElectrons, ElectronV & electrons,
                                 const edm::InputTag& metLabel, const bool & adjMET, pat::MET & met,
                                 const edm::InputTag& pfLabel){
  convertElectrons(patElectrons, electrons);
  getMET(event, metLabel, met);
  if(adjMET) AdjustMET(event, electrons, pfLabel, met);
}

void WPrimeUtil::getMuonsMET(const edm::EventBase & event,
                             const vector<pat::Muon    >& patMuons, const uint& muAlgo, MuonV & muons,
                             const edm::InputTag& metLabel, const bool & adjMET, pat::MET & met,
                             const edm::InputTag& pfLabel){
  convertMuons(patMuons, muAlgo, muons);
  getMET(event, metLabel, met);
  if(adjMET) AdjustMET(event, muons, pfLabel, met);
}

void WPrimeUtil::getLeptonsMET(const edm::EventBase & event, 
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

bool WPrimeUtil::Match(const heep::Ele & p1, const heep::Ele & p2){
  return Match(p1.patEle(), p2.patEle());
}

// get GEN-level transverse mass for lepton + neutrino;
// using delta-R matching requirement between RECO-lepton and GEN-lepton
float WPrimeUtil::getGenWprimeMt(edm::EventBase const& event, int pdgId_lepton,
				 int pdgId_neutrino, const reco::Candidate * lepton)
{
  float ret = -9999;
  const reco::Candidate * genLepton = 0;
  const reco::Candidate * genNeutrino = 0;
  genLepton = getGenParticle(event, pdgId_lepton, PDG_ID_WPRIME);
  if(!genLepton)
    {
      cerr << " *** Oops! No lepton with pdgId = " << pdgId_lepton << endl;
      abort();
    }
  genNeutrino = getGenParticle(event, pdgId_neutrino, PDG_ID_WPRIME);
  if(!genNeutrino)
    {
      cerr << " *** Oops! No neutrino with pdgId = " << pdgId_neutrino << endl;
      abort();
    }
  float dR = reco::deltaR(*genLepton, *lepton);
  // do we really need a delta-R match?
  if(dR < 0.1)
    {    
      WCandidate gen(*genLepton, *genNeutrino);
      ret = gen.mt();
    }
  return ret;
}

