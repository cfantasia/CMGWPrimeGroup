#ifndef _common_fit_h_
#define _common_fit_h_

// signal-free, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0 TeV
const unsigned mass_points = 9; 
const unsigned num_ref_plots = mass_points + 1; // equals to wprime points + bgd 

// used for filename storing fit results
string data_desc[mass_points] = {"NoSignal", "0.8TeV", "1.0TeV", "1.1TeV",
				 "1.2TeV", "1.3TeV", "1.4TeV", "1.5TeV",
				 "2.0TeV"};

// used for description of the tree of the fit results
const string histo_desc[num_ref_plots] =
  {" Signal-free", " W ' (0.8 TeV)", " W ' (1.0 TeV)", " W ' (1.1 TeV)", 
   " W ' (1.2 TeV)", " W ' (1.3 TeV)", " W ' (1.4 TeV)", " W ' (1.5 TeV)", 
   " W ' (2.0 TeV)", " SM background"}; 

// this is the mass point to be used for fits of (early) data 
// or background-only distributions. The value determines the intial value of 
// the mass, but more importantly, the resolution function. Should we set to 
// expected limit of Wprime reach? Something to think about...
const unsigned mass_point_for_data = 3; // 1.1 TeV

// resolution for signal-free (or early data) sample to be specified later
// in the code via mass_point_for_data parameter
string gname_pre[mass_points] = {"NotSetHere", "g08", "g10", "g11", "g12",
				 "g13", "g14", "g15", "g20"};

// signal-free (mass_option=0) corresponds to be specified later
// in the code via mass_point_for_data parameter
float all_masses[mass_points] = {-9999, 800, 1000, 1100, 1200, 1300, 1400,
				 1500, 2000};

// ==================================
// Fit settings for users begin here
// ==================================

// Set integrated luminosity for toy-MC ensemble studies 
// Will use this value to scale histograms in Wprime_analysis.root 
// by comparing with reference contained in "lumi_ipb" histogram
// (NB: W MC-statistics are enough for up to ~100 pb^-1)
//const float integ_lumi = 100.0; // in pb^-1
const float integ_lumi = 21.78; // in pb^-1

// select algorithm option
const unsigned algo_option = 3; // 0: "glb", 1: "trk", 2: "tpfms", 3: "ckt", "pic", "tmr"

// choose whether to perform fit
const bool doFits = true; // if false: produce only plots, no fitting

// choose whether to use ROOT's built-in log-likelihood (false)
//  or custom chi^2 (true)
const bool custom_chi2 = true;

// if using custom chi^2: 
// choose between Almeida-Barbi-do Vale (true) or Baker-Cousis (false) methods
const bool ABV_chi2 = false; 

// choose whether to fix total # of bgd+sig events to number of histogram entries 
const bool fixNtot = true;

// choose whether to apply limits on wprime mass
const bool use_wprime_mass_limits = false;
const float upper_wprime_mass_limit = 2500; // in GeV

// if true, background is modeled by Landau, otherwise by RBW
const bool isBgdLandau = false;

// choose whether to fix fudge factor (see CMS AN-2009/157 for details)
const bool fixFudge = true; // relative increase of W' width in fit
// initial (or fixed) value for fudge factor
const float Fudge = 1.0; // relative increase of W' width

#endif // #define _common_fit_h_

