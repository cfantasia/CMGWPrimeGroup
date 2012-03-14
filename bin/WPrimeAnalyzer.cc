#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <TSystem.h>

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "FWCore/ParameterSet/interface/ProcessDesc.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "UserCode/CMGWPrimeGroup/interface/WPrimeUtil.h"
#include "UserCode/CMGWPrimeGroup/interface/MuMETAnalyzer.h"
#include "UserCode/CMGWPrimeGroup/interface/EleMETAnalyzer.h"
#include "UserCode/CMGWPrimeGroup/interface/WgammaAnalyzer.h"
#include "UserCode/CMGWPrimeGroup/interface/WZAnalyzer.h"
#include "UserCode/CMGWPrimeGroup/interface/HadronicVZAnalyzer.h"
#include "UserCode/CMGWPrimeGroup/interface/HadronicVWAnalyzer.h"
#include "UserCode/CMGWPrimeGroup/interface/TBAnalyzer.h"
#include "UserCode/CMGWPrimeGroup/interface/TTbarAnalyzer.h"

using std::cout; using std::cerr; using std::endl;

AnalyzerBase* getAnalyzer(char * cfg_file, int fileToRun);

int main(int argc, char* argv[]) 
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  TH1::SetDefaultSumw2();
  
  // only allow one argument for this simple example which should be the
  // the python cfg file
  if ( argc < 2 ) {
    std::cout << "Usage : " << argv[0] << " [parameters.py] <fileToRun>" << std::endl;
    return -1;
  }

  char * config_file = argv[1];
  int fileToRun = argc > 2 ? atoi(argv[2]) : -1;

  
  AnalyzerBase* wprimeAnalyzer = getAnalyzer(config_file, fileToRun);
  wprimeAnalyzer->run();

  if(wprimeAnalyzer) delete wprimeAnalyzer;

  return 0;
}

// parse configuration, extract parameters
AnalyzerBase* getAnalyzer(char * cfg_file, int fileToRun)
{
  PythonProcessDesc builder(cfg_file);
  const edm::ParameterSet& cfg = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("WprimeAnalyzer");
  
  // now get each parameter
  bool runMuMETAnalysis  = cfg.getParameter<bool>("runMuMETAnalysis" );
  bool runElMETAnalysis  = cfg.getParameter<bool>("runElMETAnalysis" );
  bool runTBAnalysis     = cfg.getParameter<bool>("runTBAnalysis"    );
  bool runWgammaAnalysis = cfg.getParameter<bool>("runWgammaAnalysis"); 
  bool runWZAnalysis     = cfg.getParameter<bool>("runWZAnalysis"    );
  bool runHadVZAnalysis  = cfg.getParameter<bool>("runHadVZAnalysis" );
  bool runHadVWAnalysis  = cfg.getParameter<bool>("runHadVWAnalysis" );
  bool runTTbarAnalysis  = cfg.getParameter<bool>("runTTbarAnalysis" );

  if     (runMuMETAnalysis) return new MuMETAnalyzer (cfg, fileToRun);
  else if(runElMETAnalysis) return new EleMETAnalyzer(cfg, fileToRun);
  else if(runTBAnalysis)    return new TBAnalyzer    (cfg, fileToRun);
  else if(runWgammaAnalysis)return new WgammaAnalyzer(cfg, fileToRun);
  else if(runWZAnalysis)    return new WZAnalyzer    (cfg, fileToRun);
  else if(runHadVZAnalysis) return new HadronicVZAnalyzer(cfg, fileToRun);
  else if(runHadVWAnalysis) return new HadronicVWAnalyzer(cfg, fileToRun);
  else if(runTTbarAnalysis) return new TTbarAnalyzer(cfg, fileToRun);
  
  cerr<<" You haven't enabled any analysis modes!\n";
  abort();
  return NULL;
}
