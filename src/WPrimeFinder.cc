#include "UserCode/CMGWPrimeGroup/interface/WPrimeFinder.h"

#include <iostream>

using std::cout; using std::cerr; using std::endl; using std::vector;
using std::string;

#include <TH1F.h>

// constructor: needs configuration file to set things up
WPrimeFinder::WPrimeFinder(char * config_file)
{
  getConfiguration(config_file);
  wprimeUtil->getInputFiles(inputFiles);
}

WPrimeFinder::~WPrimeFinder()
{
  if(muMETAnalyzer) delete muMETAnalyzer;
  if(WmunugammaAnalyzer) delete WmunugammaAnalyzer;
  if(wzAnalyzer) delete wzAnalyzer;
  delete wprimeUtil;
  outLogFile_.close(); 
}

// parse configuration, extract parameters
void WPrimeFinder::getConfiguration(char * cfg_file)
{
  PythonProcessDesc builder(cfg_file);
  const edm::ParameterSet& cfg = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("WprimeAnalyzer");
  
  // now get each parameter
  /////  inputFiles_  = cfg.getParameter<vector<string> >("fileNames") ;
  outputFile_  = cfg.getParameter<string  >("outputFile" );
  logFile_  = cfg.getParameter<string  >("logFile" );
  reportAfter_ = cfg.getParameter<unsigned int>("reportAfter");
  maxEvents_   = cfg.getParameter<int>("maxEvents") ;
  useJSON_   = cfg.getParameter<bool>("useJSON") ;
  countGenEvts_ = cfg.getParameter<bool>("countGenEvts");
  genParticles_ = cfg.getParameter<edm::InputTag>("genParticles" );
  doRecoilCorrectionForW_ = cfg.getParameter<bool>("doRecoilCorrectionForW");
  runMuMETAnalysis_ = cfg.getParameter<bool>("runMuMETAnalysis" );
  runElMETAnalysis_ = cfg.getParameter<bool>("runElMETAnalysis" );
  runWZAnalysis_ = cfg.getParameter<bool>("runWZAnalysis" );
  runHadVZAnalysis_ = cfg.getParameter<bool>("runHadVZAnalysis" );
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

  outLogFile_.open(logFile_.c_str());
  WPrimeUtil::CheckStream(outLogFile_, logFile_);

  ctrNames_ = (cfg.getParameter<vstring>("eventCounters"));

  MCPUDistFile_   = cfg.getParameter<string>("MCPUDistFile" );
  MCPUDistHist_   = cfg.getParameter<string>("MCPUDistHist" );
  DataPUDistFile_ = cfg.getParameter<string>("DataPUDistFile" );
  DataPUDistHist_ = cfg.getParameter<string>("DataPUDistHist" );
  
  std::vector<edm::EventID> vEventsToDebug = cfg.getParameter<std::vector<edm::EventID> >("vEventsToDebug");
  
  wprimeUtil = new WPrimeUtil(outputFile_.c_str(), genParticles_, sample_cross_sections);
  wprimeUtil->SetLumiWeights(MCPUDistFile_, DataPUDistFile_, MCPUDistHist_, DataPUDistHist_);
  wprimeUtil->SetEventsToDebug(vEventsToDebug);

  if(runMuMETAnalysis_)
    muMETAnalyzer = new MuMETAnalyzer(cfg, wprimeUtil);
  else
    muMETAnalyzer = 0;

  if(runElMETAnalysis_)
    eleMETAnalyzer = new EleMETAnalyzer(cfg, wprimeUtil);
  else
    eleMETAnalyzer = 0;

  if(runWgammaAnalysis_)
    WmunugammaAnalyzer = new WgammaAnalyzer(cfg, wprimeUtil);
  else
    WmunugammaAnalyzer = 0;

  if(runWZAnalysis_)
    wzAnalyzer = new WZAnalyzer(cfg, wprimeUtil);
  else
    wzAnalyzer = 0;

  if(runHadVZAnalysis_)
    hadvzAnalyzer = new HadronicVZAnalyzer(cfg, wprimeUtil);
  else
    hadvzAnalyzer = 0;
}

// operations to be done when changing input file (e.g. create new histograms)
void WPrimeFinder::beginFile(vector<wprime::InputFile>::const_iterator it)
{
  bool shouldCorrectMt = 
    ((it->samplename=="W" || it->samplename=="Wlowpt") 
     && doRecoilCorrectionForW_);
  wprimeUtil->setApplyMETCorrection(shouldCorrectMt);

  wprimeUtil->setSampleName(it->samplename);
  wprimeUtil->setWeight(it->weight);
  wprimeUtil->setRunningOnData();

  // call beginFile for each finder here
  if(runMuMETAnalysis_)
      muMETAnalyzer->beginFile(it);
  if(runElMETAnalysis_)
      eleMETAnalyzer->beginFile(it);
  if(runWgammaAnalysis_) 
      WmunugammaAnalyzer->beginFile(it);
  if(runWZAnalysis_) 
      wzAnalyzer->beginFile(it);
  if(runHadVZAnalysis_)
    hadvzAnalyzer->beginFile(it);

}

void WPrimeFinder::eventLoop(edm::EventBase const & event)
{
  wprimeUtil->setHadronicMETCalculated(false);

  if(runMuMETAnalysis_)
    muMETAnalyzer->eventLoop(event);

  if(runElMETAnalysis_)
    eleMETAnalyzer->eventLoop(event);

  if(runWgammaAnalysis_)
      WmunugammaAnalyzer->eventLoop(event);

  if(runWZAnalysis_)
      wzAnalyzer->eventLoop(event);

  if(runHadVZAnalysis_)
    hadvzAnalyzer->eventLoop(event);
}



void WPrimeFinder::run()
{
  int ievt_all=0;  int ievt_skipped = 0;
  unsigned i_sample = 1;
  vector<wprime::InputFile>::iterator it;

  for(it = inputFiles.begin(); it != inputFiles.end(); ++it, ++i_sample){
    int ievt=0;  
    cout << "\n Opening sample " << it->samplename << " ... ";
    fwlite::ChainEvent ev(it->pathnames);
    it->Nact_evt = ev.size();
    cout<<" Done. \n";
  
    if(wprimeUtil->runningOnData())
      // Nprod_evt presumably contains the # of events before any filtering
      // that results in Nact_evt (< Nprod_evt) events contained in the file.
      // For data, we tend not to know how many events we started with,
      // so just assume pre-selection efficiency = 100%;
      // this affects only the efficiency calculations printed
      // at the end of the job - nothing else!
      it->Nprod_evt = it->Nact_evt;
    
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
        updateEventCounts(ev, nEvents_, 
                          runNumber, lumiID, 
                          ctrNames_, false);
      
      // skip event if maximal number of events per input file is reached 
      if(maxEvents_>0 &&  ievt > maxEvents_) continue;
      
      // simple event counter
      if(reportAfter_!=0 ? (ievt>0 && ievt%reportAfter_==0) : false) 
        cout << " Processing event: " << ievt << " or " 
             << 100.*ievt/it->Nact_evt << "% of input file"
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
    endFile(it);
    
  } // loop over input files
  cout<<"Done with Input Samples\n";

  wprimeUtil->getFileService()->cd(); 
  TH1F * h = new TH1F("lumi_ipb", "Integrated luminosity in pb^{-1}", 1, 0, 1);
  h->SetBinContent(1, wprimeUtil->getLumi_ipb());
  //  h->Write();
  
  endAnalysis();

}

// operations to be done when closing input file 
// (e.g. save histograms, print summary)
void WPrimeFinder::endFile(vector<wprime::InputFile>::const_iterator it)
{
   // call endFile for each finder here
  if(runMuMETAnalysis_)
    muMETAnalyzer->endFile(it, outLogFile_);
  if(runElMETAnalysis_)
    eleMETAnalyzer->endFile(it, outLogFile_);
  if(runWgammaAnalysis_) 
      WmunugammaAnalyzer->endFile(it, outLogFile_);  
  if(runWZAnalysis_) 
      wzAnalyzer->endFile(it, outLogFile_);  
  if(runHadVZAnalysis_)
    hadvzAnalyzer->endFile(it, outLogFile_);  
}

// e.g. print summmary of expected events for all samples
void WPrimeFinder::endAnalysis()
{
  if(runMuMETAnalysis_)
    muMETAnalyzer->endAnalysis(outLogFile_);
  if(runElMETAnalysis_)
    eleMETAnalyzer->endAnalysis(outLogFile_);
  if(runWgammaAnalysis_)
      WmunugammaAnalyzer->endAnalysis(outLogFile_);
  if(runWZAnalysis_)
      wzAnalyzer->endAnalysis(outLogFile_);
  if(runHadVZAnalysis_)
    hadvzAnalyzer->endAnalysis(outLogFile_);
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
