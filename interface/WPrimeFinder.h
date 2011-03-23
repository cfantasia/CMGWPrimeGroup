#ifndef _wprime_finder_h_
#define _wprime_finder_h_

#include "FWCore/ParameterSet/interface/ProcessDesc.h"

#include "UserCode/CMGWPrimeGroup/interface/WPrimeUtil.h"
#include "UserCode/CMGWPrimeGroup/interface/MuMETAnalyzer.h"


class WPrimeFinder
{
 public:
  // constructor: needs configuration file to set things up
  explicit WPrimeFinder(char * config_file);
  ~WPrimeFinder();
  
  void run();

 private: 
  // parse configuration, extract parameters
  void getConfiguration(char * config_file);

  // operations to be done when changing input file (e.g. create new histograms)
  void beginFile(std::vector<wprime::InputFile>::const_iterator it);

  // operations to be done when closing input file 
  // (e.g. save histograms, print summary)
  void endFile(std::vector<wprime::InputFile>::const_iterator it);

  // e.g. print summary of expected events for all samples
  void endAnalysis();

  ofstream out_;

  // print out event # 
  unsigned int reportAfter_;
  // maximum # of events to process (set to <0 for processing all events)
  int maxEvents_;

  void eventLoop(edm::EventBase const & event);

  std::vector<wprime::InputFile> inputFiles; 

  std::string outputFile_;

  // enable/disable analysis in specific channels
  bool runMuMETAnalysis_;
  bool runElMETAnalysis_;
  bool runWZAnalysis_;
  bool runTBAnalysis_;

  bool doRecoilCorrectionForW_;

  edm::InputTag genParticles_;

  MuMETAnalyzer * muMETAnalyzer;
  WPrimeUtil * wprimeUtil;

};

#endif // #define _wprime_finder_h_
