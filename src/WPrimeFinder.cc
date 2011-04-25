#include "UserCode/CMGWPrimeGroup/interface/WPrimeFinder.h"

#include <iostream>

using std::cout; using std::cerr; using std::endl;
using std::string;

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

// constructor: needs configuration file to set things up
WPrimeFinder::WPrimeFinder(char * config_file)
  : out_("event_counts.txt")
{
  if(!out_) { 
    cout << " Cannot open output text file " << endl;
    abort();
  } 
  
  getConfiguration(config_file);
  wprimeUtil->getInputFiles(inputFiles);
}

WPrimeFinder::~WPrimeFinder()
{
  if(muMETAnalyzer) delete muMETAnalyzer;
  if(WmunugammaAnalyzer) delete WmunugammaAnalyzer;
  delete wprimeUtil;
  out_.close();
}

// parse configuration, extract parameters
void WPrimeFinder::getConfiguration(char * cfg_file)
{
  PythonProcessDesc builder(cfg_file);
  const edm::ParameterSet& cfg = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("WprimeAnalyzer");
  
  // now get each parameter
  /////  inputFiles_  = cfg.getParameter<std::vector<string> >("fileNames") ;
  outputFile_  = cfg.getParameter<string  >("outputFile" );
  reportAfter_ = cfg.getParameter<unsigned int>("reportAfter");
  maxEvents_   = cfg.getParameter<int>("maxEvents") ;
  genParticles_ = cfg.getParameter<edm::InputTag>("genParticles" );
  doRecoilCorrectionForW_ = cfg.getParameter<bool>("doRecoilCorrectionForW");
  runMuMETAnalysis_ = cfg.getParameter<bool>("runMuMETAnalysis" );
  runElMETAnalysis_ = cfg.getParameter<bool>("runElMETAnalysis" );
  runWZAnalysis_ = cfg.getParameter<bool>("runWZAnalysis" );
  runTBAnalysis_ = cfg.getParameter<bool>("runTBAnalysis" );
  runWgammaAnalysis_ =cfg.getParameter<bool>("runWgammaAnalysis"); 
  string sample_cross_sections = cfg.getParameter<string>("sample_cross_sections");

  wprimeUtil = new WPrimeUtil(outputFile_.c_str(), genParticles_, sample_cross_sections);

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

}

// operations to be done when changing input file (e.g. create new histograms)
void WPrimeFinder::beginFile(std::vector<wprime::InputFile>::const_iterator it)
{

  bool shouldCorrectMt = 
    ((it->samplename=="W" || it->samplename=="Wlowpt") 
     && doRecoilCorrectionForW_);
  wprimeUtil->setApplyMETCorrection(shouldCorrectMt);

  wprimeUtil->setSampleName(it->samplename);
  wprimeUtil->setWeight(it->weight);

  // call beginFile for each finder here
  if(runMuMETAnalysis_)
      muMETAnalyzer->beginFile(it);
  if(runElMETAnalysis_)
      eleMETAnalyzer->beginFile(it);
  if(runWgammaAnalysis_) 
      WmunugammaAnalyzer->beginFile(it);
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
}


void WPrimeFinder::run()
{
  int ievt=0;  
  std::vector<wprime::InputFile>::iterator it;
  for(it = inputFiles.begin(); it != inputFiles.end(); ++it){
    // loop over input files
    TFile * inFile = TFile::Open(it->pathname.c_str());
    if(!inFile || inFile->IsZombie())
      {
	cout << " *** Can't find file " << it->pathname
	     << " ! Aborting... " << endl;
	abort();
      }

    it->file = inFile;
    TTree * events = (TTree *) inFile->Get("Events");
    assert(events);

    it->Nact_evt = events->GetEntries();
    if(it->samplename == "data")
        // Nprod_evt presumably contains the # of events before any filtering
        // that results in Nact_evt (< Nprod_evt) events contained in the file.
        // For data, we tend not to know how many events we started with,
        // so just assume pre-selection efficiency = 100%;
        // this affects only the efficiency calculations printed
        // at the end of the job - nothing else!
        it->Nprod_evt = it->Nact_evt;
 	  	 
    cout << " Opened file " << it->pathname << " with " << it->Nact_evt
         << " events" << endl;

    beginFile(it);
    fwlite::Event ev(inFile);
    for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){// loop over events
      edm::EventBase const & event = ev;
      // break loop if maximal number of events is reached 
      if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
      // simple event counter
      if(reportAfter_!=0 ? (ievt>0 && ievt%reportAfter_==0) : false) 
	cout << "  Processing event: " << ievt << endl;
      eventLoop(event);
    } // loop over events
    endFile(it);
    
     // break loop if maximal number of events is reached:
    // this has to be done twice to stop the file loop as well
    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
    
  } // loop over input files

  wprimeUtil->getFileService()->cd(); 
  TH1F * h = new TH1F("lumi_ipb", "Integrated luminosity in pb^{-1}", 1, 0, 1);
  h->SetBinContent(1, wprimeUtil->getLumi_ipb());
  //  h->Write();
  
  endAnalysis();

}

// operations to be done when closing input file 
// (e.g. save histograms, print summary)
void WPrimeFinder::endFile(std::vector<wprime::InputFile>::const_iterator it)
{
  // close input file
  it->file->Close();
  
   // call endFile for each finder here
  if(runMuMETAnalysis_)
    muMETAnalyzer->endFile(it, out_);
  if(runElMETAnalysis_)
    eleMETAnalyzer->endFile(it, out_);
  if(runWgammaAnalysis_) 
      WmunugammaAnalyzer->endFile(it, out_);  
}

// e.g. print summmary of expected events for all samples
void WPrimeFinder::endAnalysis()
{
  if(runMuMETAnalysis_)
    muMETAnalyzer->endAnalysis(out_);
  if(runElMETAnalysis_)
    eleMETAnalyzer->endAnalysis(out_);
  if(runWgammaAnalysis_)
      WmunugammaAnalyzer->endAnalysis(out_);
}
