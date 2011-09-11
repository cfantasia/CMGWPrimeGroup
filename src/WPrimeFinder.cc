#include "UserCode/CMGWPrimeGroup/interface/WPrimeFinder.h"

#include <iostream>

using std::cout; using std::cerr; using std::endl; using std::vector;
using std::string;

#include <TH1F.h>

// constructor: needs configuration file to set things up
WPrimeFinder::WPrimeFinder(char * config_file, int fileToRun)
{
  getConfiguration(config_file,fileToRun);
  wprimeUtil->getInputFiles(inputFiles);
  if(fileToRun != -1){
    if(fileToRun < (int)inputFiles.size()){
      inputFiles.assign(1,inputFiles[fileToRun]);
    }else{
      cerr<<"You asked for sample "<<fileToRun
          <<" but only "<<inputFiles.size()
          <<" are listed!\n";
      inputFiles.clear();
      abort();
    }
  }
}

WPrimeFinder::~WPrimeFinder()
{
  if(wprimeAnalyzer) delete wprimeAnalyzer;
  if(wprimeUtil) delete wprimeUtil;
  outLogFile_.close(); 
}

// parse configuration, extract parameters
void WPrimeFinder::getConfiguration(char * cfg_file, int fileToRun)
{
  PythonProcessDesc builder(cfg_file);
  const edm::ParameterSet& cfg = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("WprimeAnalyzer");
  
  // now get each parameter
  /////  inputFiles_  = cfg.getParameter<vector<string> >("fileNames") ;
  outputFile_  = cfg.getParameter<string  >("outputFile" );
  logFile_  = cfg.getParameter<string  >("logFile" );
  reportAfter_ = cfg.getParameter<unsigned int>("reportAfter");
  maxEvents_   = cfg.getParameter<int>("maxEvents");
  useJSON_   = cfg.getParameter<bool>("useJSON") ;
  countGenEvts_ = cfg.getParameter<bool>("countGenEvts");
  genLabel_ = cfg.getParameter<edm::InputTag>("genParticles" );
  pfLabel_ = cfg.getParameter<edm::InputTag>("particleFlow" );
  pileupLabel_ = cfg.getParameter<edm::InputTag>("pileupTag" );
  doRecoilCorrectionForW_ = cfg.getParameter<bool>("doRecoilCorrectionForW");
  runMuMETAnalysis_ = cfg.getParameter<bool>("runMuMETAnalysis" );
  runElMETAnalysis_ = cfg.getParameter<bool>("runElMETAnalysis" );
  runWZAnalysis_ = cfg.getParameter<bool>("runWZAnalysis" );
  runHadVZAnalysis_ = cfg.getParameter<bool>("runHadVZAnalysis" );
  runHadVWAnalysis_ = cfg.getParameter<bool>("runHadVWAnalysis" );
  runTBAnalysis_ = cfg.getParameter<bool>("runTBAnalysis" );
  runWgammaAnalysis_ =cfg.getParameter<bool>("runWgammaAnalysis"); 
  string sample_cross_sections = cfg.getParameter<string>("sample_cross_sections");
  edm::ParameterSet const& inputs = cfg.getParameter<edm::ParameterSet>("inputs");
  if ( inputs.exists("lumisToProcess") ) 
    {
      vector<edm::LuminosityBlockRange> const & lumisTemp =
	inputs.getUntrackedParameter<vector<edm::LuminosityBlockRange> > ("lumisToProcess");
      jsonVector.resize( lumisTemp.size() );
      copy( lumisTemp.begin(), lumisTemp.end(), jsonVector.begin() );
    }

  if(fileToRun != -1){
    logFile_ = Form("Sample%i_%s",fileToRun,logFile_.c_str()); 
    outputFile_ = Form("Sample%i_%s",fileToRun,outputFile_.c_str()); 
  }


  outLogFile_.open(logFile_.c_str());
  WPrimeUtil::CheckStream(outLogFile_, logFile_);

  ctrNames_ = (cfg.getParameter<vstring>("eventCounters"));

  MCPUDistFile_   = cfg.getParameter<string>("MCPUDistFile" );
  MCPUDistHist_   = cfg.getParameter<string>("MCPUDistHist" );
  DataPUDistFile_ = cfg.getParameter<string>("DataPUDistFile" );
  DataPUDistHist_ = cfg.getParameter<string>("DataPUDistHist" );
  
  std::vector<edm::EventID> vEventsToDebug = cfg.getParameter<std::vector<edm::EventID> >("vEventsToDebug");
  
  wprimeUtil = new WPrimeUtil(outputFile_.c_str(), genLabel_, pfLabel_, sample_cross_sections);
  wprimeUtil->setLumiWeights(MCPUDistFile_, DataPUDistFile_, MCPUDistHist_, DataPUDistHist_);
  wprimeUtil->setEventsToDebug(vEventsToDebug);

  if     (runMuMETAnalysis_) wprimeAnalyzer = new MuMETAnalyzer(cfg, wprimeUtil);
  else if(runElMETAnalysis_) wprimeAnalyzer = new EleMETAnalyzer(cfg, wprimeUtil);
  else if(runWgammaAnalysis_)wprimeAnalyzer = new WgammaAnalyzer(cfg, wprimeUtil);
  else if(runWZAnalysis_)    wprimeAnalyzer = new WZAnalyzer(cfg, wprimeUtil);
  else if(runHadVZAnalysis_) wprimeAnalyzer = new HadronicVZAnalyzer(cfg, wprimeUtil);
  else if(runHadVWAnalysis_) wprimeAnalyzer = new HadronicVWAnalyzer(cfg, wprimeUtil);
  else if(runTBAnalysis_)    wprimeAnalyzer = new TBAnalyzer(cfg, wprimeUtil);
}

// operations to be done when changing input file (e.g. create new histograms)
void WPrimeFinder::beginFile(vector<wprime::InputFile>::iterator it)
{
  bool shouldCorrectMt = 
    ((it->samplename=="W" || it->samplename=="Wlowpt") 
     && doRecoilCorrectionForW_);
  wprimeUtil->setApplyHadronicRecoilCorrection(shouldCorrectMt);

  wprimeUtil->setSampleName(it->samplename);
  wprimeUtil->setSampleWeight(it->weight);
  wprimeUtil->setRunningOnData();
  wprimeUtil->resetWarnings();

  if(wprimeUtil->runningOnData())
    // Nprod_evt presumably contains the # of events before any filtering
    // that results in Nact_evt (< Nprod_evt) events contained in the file.
    // For data, we tend not to know how many events we started with,
    // so just assume pre-selection efficiency = 100%;
    // this affects only the efficiency calculations printed
    // at the end of the job - nothing else!
    it->Nprod_evt = it->Nact_evt;

  // call beginFile for each finder here
  wprimeAnalyzer->beginFile(it);
}

void WPrimeFinder::eventLoop(edm::EventBase const & event)
{
  wprimeUtil->setHadronicMETcalculated(false);

  if(wprimeUtil->runningOnData()){
    wprimeUtil->setWeight(wprimeUtil->getSampleWeight());
  }else{
    event.getByLabel(pileupLabel_, PupH_);
    float PU_Weight = wprimeUtil->getPUWeight3BX(*PupH_);
    wprimeUtil->setWeight(wprimeUtil->getSampleWeight() * PU_Weight);
  }

  wprimeAnalyzer->eventLoop(event);
}



void WPrimeFinder::run()
{
  int ievt_all=0;  int ievt_skipped = 0;
  unsigned i_sample = 1;
  vector<wprime::InputFile>::iterator it;

  for(it = inputFiles.begin(); it != inputFiles.end(); ++it, ++i_sample){
    int ievt=0;  
    cout << "\n Opening sample " << it->samplename 
         << "( " << it->description << " ) ... ";
    fwlite::ChainEvent ev(it->pathnames);
    it->Nact_evt = ev.size();
    cout<<" Done. \n";
  
    cout << " Opened sample " << it->samplename << " with " << it->Nact_evt
         << " events (Input file #" << i_sample << " out of " << inputFiles.size()
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
      
      if(useJSON_ && wprimeUtil->runningOnData() &&
         !jsonContainsEvent (jsonVector, event))
      {
        ++ievt_skipped;
        continue;
      }
      else
        ++ievt_all;
      
      eventLoop(event);
    } // loop over events
    if(countGenEvts_ && (int)nEvents_[0] != it->Nprod_evt) 
      cout<<"Weight Wrong: Found "<<nEvents_[0]<<" generated events and sample file lists "<<it->Nprod_evt<<endl;
    wprimeAnalyzer->endFile(it, outLogFile_);
    
  } // loop over input files
  cout<<"Done with Input Samples\n";

  wprimeUtil->getFileService()->cd(); 
  TH1F * h = new TH1F("lumi_ipb", "Integrated luminosity in pb^{-1}", 1, 0, 1);
  h->SetBinContent(1, wprimeUtil->getLumi_ipb());
  //  h->Write();
  TH1F * hFilecounter = new TH1F("hFilecounter", "counter indicates number of files merged", 1, 0, 1);
  hFilecounter->SetBinContent(1, 1);
  
  wprimeAnalyzer->endAnalysis(outLogFile_);

}

bool WPrimeFinder::jsonContainsEvent (const vector<edm::LuminosityBlockRange>&jsonVec, const edm::EventBase &event)
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
