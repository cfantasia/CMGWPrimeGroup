#ifndef _load_cuts_h__
#define _load_cuts_h__

#include "UserCode/CMGWPrimeGroup/interface/wprimeEvent.h"

// returns # of (global) muons with tracker-pt above <tracker_muon_pt>
unsigned NmuAboveThresh(float tracker_muon_pt, const wprime::Event * ev);

// returns # of jets with Et above threshold
unsigned NjetAboveThresh(float threshold, const wprime::Event * ev);

// returns # of jets with Et above threshold with angle greater than delta_phi from muon
unsigned NjetAboveThresh(float threshold, float delta_phi, 
			 const wprime::Muon * mu, const wprime::Event * ev);


#endif // #define _load_cuts_h__
