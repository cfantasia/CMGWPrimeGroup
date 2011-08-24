#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <TSystem.h>

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "UserCode/CMGWPrimeGroup/interface/WPrimeFinder.h"

using std::cout; using std::endl;

int main(int argc, char* argv[]) 
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  TH1::SetDefaultSumw2();
  
  // only allow one argument for this simple example which should be the
  // the python cfg file
  if ( argc < 2 ) {
    std::cout << "Usage : " << argv[0] << " [parameters.py] <fileToRun" << std::endl;
    return -1;
  }

  char * config_file = argv[1];
  int fileToRun = argc > 2 ? atoi(argv[2]) : -1;
  WPrimeFinder a(config_file, fileToRun);

  a.run();


  return 0;
}
