#include "UserCode/CMGWPrimeGroup/root_macros/input_file.h"

#include <vector>
using std::vector;

int loadCrossSections(vector<wprime::InputFile> & qcd_files, 
		      vector<wprime::InputFile> & z_files,
		      vector<wprime::InputFile> & w_files,
		      vector<wprime::InputFile> & top_files,
		      vector<wprime::InputFile> & wprime10_files,
		      vector<wprime::InputFile> & wprime15_files,
		      vector<wprime::InputFile> & wprime20_files
		      )
{
  // first argument: PYTHIA cross-section (pb), second argument: # of produced evts

  qcd_files.push_back(wprime::InputFile(590000  , 450000)); // 100_150
  qcd_files.push_back(wprime::InputFile( 83000  , 425000)); // 150_200
  qcd_files.push_back(wprime::InputFile( 24000  , 440000)); // 200_300
  qcd_files.push_back(wprime::InputFile(  3000  , 395000)); // 300_400
  qcd_files.push_back(wprime::InputFile(   730  , 135000)); // 400_600
  qcd_files.push_back(wprime::InputFile(    66  , 310000)); // 600_800
  qcd_files.push_back(wprime::InputFile(    12  , 130000)); // 800_1200
  qcd_files.push_back(wprime::InputFile(    0.63,  30000)); // 1200_1600
  qcd_files.push_back(wprime::InputFile(    0.54,  25000)); // 1600_2000

  z_files.push_back(wprime::InputFile(1400     , 19100)); //   30-110
  z_files.push_back(wprime::InputFile(  19     ,  4100)); //  110-200
  z_files.push_back(wprime::InputFile(   1.2   ,  1000)); //  200-300
  z_files.push_back(wprime::InputFile(   0.24  ,   900)); //  300-400
  z_files.push_back(wprime::InputFile(   0.075 ,   900)); //  400-500
  z_files.push_back(wprime::InputFile(   0.029 ,   900)); //  500-600
  z_files.push_back(wprime::InputFile(   0.013 ,   900)); //  600-700
  z_files.push_back(wprime::InputFile(   0.0063,   700)); //  700-800
  z_files.push_back(wprime::InputFile(   0.0033,   700)); //  800-900
  z_files.push_back(wprime::InputFile(   0.0019,   900)); //  900-1000
  z_files.push_back(wprime::InputFile(   0.0029,   600)); // 1000-infinity

  w_files.push_back(wprime::InputFile(11680, 1000000)); //     0- 200
  w_files.push_back(wprime::InputFile(8.679E-2, 10000)); //  200- 250
  w_files.push_back(wprime::InputFile(3.253E-2, 10000)); //  250- 300
  w_files.push_back(wprime::InputFile(1.431E-2, 10000)); //  300- 350
  w_files.push_back(wprime::InputFile(6.946E-3, 10000)); //  350- 400
  w_files.push_back(wprime::InputFile(3.653E-3, 10000)); //  400- 450
  w_files.push_back(wprime::InputFile(2.026E-3, 10000)); //  450- 500
  w_files.push_back(wprime::InputFile(1.865E-3, 5000)); //  500- 600
  w_files.push_back(wprime::InputFile(7.134E-4, 5000)); //  600- 700
  w_files.push_back(wprime::InputFile(2.949E-4, 5000)); //  700- 800
  w_files.push_back(wprime::InputFile(1.312E-4, 5000)); //  800- 900
  w_files.push_back(wprime::InputFile(6.080E-5, 5000)); //  900-1000
  w_files.push_back(wprime::InputFile(2.905E-5, 3000)); // 1000-1100
  w_files.push_back(wprime::InputFile(1.418E-5, 3000)); // 1100-1200
  w_files.push_back(wprime::InputFile(7.205E-6, 3000)); // 1200-1300
  w_files.push_back(wprime::InputFile(3.616E-6, 3000)); // 1300-1400
  w_files.push_back(wprime::InputFile(1.843E-6, 3000)); // 1400-1500
  w_files.push_back(wprime::InputFile(9.588E-7, 3000)); // 1500-1600
  w_files.push_back(wprime::InputFile(4.932E-7, 3000)); // 1600-1700
  w_files.push_back(wprime::InputFile(2.552E-7, 3000)); // 1700-1800
  w_files.push_back(wprime::InputFile(1.304E-7, 3000)); // 1800-1900
  w_files.push_back(wprime::InputFile(6.635E-8, 3000)); // 1900-2000

  top_files.push_back(wprime::InputFile(450, 259000)); // top

  wprime10_files.push_back(wprime::InputFile(1.5, 6000)); // 1 TeV
  wprime15_files.push_back(wprime::InputFile(0.24, 6000)); // 1.5 TeV
  wprime20_files.push_back(wprime::InputFile(0.0535, 6000)); // 2.0 TeV

  return 0;
}

