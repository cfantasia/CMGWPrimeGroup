#ifndef _wprime_finder_h_
#define _wprime_finder_h_

#include "FWCore/ParameterSet/interface/ProcessDesc.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "UserCode/CMGWPrimeGroup/interface/WPrimeUtil.h"
#include "UserCode/CMGWPrimeGroup/interface/MuMETAnalyzer.h"
#include "UserCode/CMGWPrimeGroup/interface/EleMETAnalyzer.h"
#include "UserCode/CMGWPrimeGroup/interface/WgammaAnalyzer.h"
#include "UserCode/CMGWPrimeGroup/interface/WZAnalyzer.h"
#include "UserCode/CMGWPrimeGroup/interface/HadronicVZAnalyzer.h"

class WPrimeFinder
{
 public:
  // constructor: needs configuration file to set things up
  explicit WPrimeFinder(char * config_file, int fileToRun);
  ~WPrimeFinder();
  
  void run();

 private: 
  // parse configuration, extract parameters
  void getConfiguration(char * config_file, int fileToRun);

  // operations to be done when changing input file (e.g. create new histograms)
  void beginFile(std::vector<wprime::InputFile>::const_iterator it);

  // operations to be done when closing input file 
  // (e.g. save histograms, print summary)
  void endFile(std::vector<wprime::InputFile>::const_iterator it);

  // e.g. print summary of expected events for all samples
  void endAnalysis();

  ofstream outLogFile_;

  // print out event # 
  unsigned int reportAfter_;
  // maximum # of events to process (set to <0 for processing all events)
  int maxEvents_;
  // Should we use the json file
  bool useJSON_;
  // Should we count the number of gen evts in patTuple?
  bool countGenEvts_;

  void eventLoop(edm::EventBase const & event);

  std::vector<wprime::InputFile> inputFiles; 

  std::string outputFile_;
  std::string logFile_;

  //Variables for counting # of gen events
  std::vector<uint> nEvents_;
  std::vector<std::string> ctrNames_;

  // enable/disable analysis in specific channels
  bool runMuMETAnalysis_;
  bool runElMETAnalysis_;
  bool runWZAnalysis_;
  bool runHadVZAnalysis_;
  bool runTBAnalysis_;
  bool runWgammaAnalysis_;
  bool doRecoilCorrectionForW_;

  edm::InputTag genParticles_;

  MuMETAnalyzer * muMETAnalyzer;
  EleMETAnalyzer * eleMETAnalyzer;
  WgammaAnalyzer * WmunugammaAnalyzer;
  WZAnalyzer * wzAnalyzer;
  HadronicVZAnalyzer * hadvzAnalyzer;
  WPrimeUtil * wprimeUtil;

  std::vector< edm::LuminosityBlockRange > jsonVector;

  bool jsonContainsEvent (const std::vector<edm::LuminosityBlockRange>&jsonVec,
			  const edm::EventBase &event);

  std::string MCPUDistFile_;
  std::string MCPUDistHist_;
  std::string DataPUDistFile_;
  std::string DataPUDistHist_;
};

#endif // #define _wprime_finder_h_


