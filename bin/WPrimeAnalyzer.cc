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
  
  // only allow one argument for this simple example which should be the
  // the python cfg file
  if ( argc < 2 ) {
    std::cout << "Usage : " << argv[0] << " [parameters.py]" << std::endl;
    return -1;
  }

  char * config_file = argv[1];
  WPrimeFinder a(config_file);

  a.run();


  return 0;
}
