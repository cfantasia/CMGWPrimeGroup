#ifndef _wprimeFitter_signalDescriptions_h_
#define _wprimeFitter_signalDescriptions_h_

const unsigned Nsignal = 9; // 1 bgd-only + 8 true-signal mass points

// the first point corresponds to bgd-only distributions (ie. not a real W')
// the first mass value is used for the S+B hypothesis fit;
// may have to revisit this argument
const string dirname[Nsignal] = {"wprime2.5", "wprime1.5", "wprime1.6", "wprime1.7", "wprime1.8", "wprime1.9", "wprime2.0", "wprime2.1", "wprime2.2"};


const string desc[Nsignal] = {"SM", "W' (1.5 TeV)", "W' (1.6 TeV)", 
			      "W' (1.7 TeV)", "W' (1.8 TeV)", "W' (1.9 TeV)",
			      "W' (2.0 TeV)", "W' (2.1 TeV)", "W' (2.2 TeV)"};

const float WprimeMass[Nsignal] = {2500, 1500, 1600, 1700, 1800, 1900, 2000, 
				   2100, 2200};

// SSM cross-sections for leptonic channels on pb;
// these will be divided by scale factor
// first point corresponds to bgd-only case: cross-section = 0
const float xsec_lep[Nsignal] = {
  0, 0.0949, 0.0633, 0.0424, 0.0285, 0.0194, 0.0135, 0.009373, 0.006605};


// color for plotting Z-values of different signals
Color_t color[Nsignal] = {kBlack, kRed, kBlue, kViolet, kGreen, kOrange, kYellow, kViolet, kOrange};


// this structure is used to keep track of expected bgd, sig 
// and total(=sig+bgd) of events for each W'-mass & W' cross-section scenario;
// useful for debugging
struct Nexp{
  float Ntot;
  float Nbgd;
  float Nsig;
};

//Nexp Nevt[Nsignal];

const float M_W = 80.399;
const float G_W = 2.085;

// strings to construct histogram names for signal distributions

// electron channel
const string resHist_name_el = "Res_el";
const string mtHist_name_el = "hTM_met";
const string genMtHist_name_el = "GenMt_el";
const string bgd_name_el = "hbgd20";

const string file_SIG_el = "Wprime_ElMET_25GeV_MC_signal.root";
const string file_BGD_el = "Tbg.root";
const string file_data_el = "Wprime_data.root";

// muon channel
const string resHist_name_mu = "Res_mu";
const string mtHist_name_mu = "hTMcmkt_met";
const string genMtHist_name_mu = "GenMt_mu";
const string bgd_name_mu = "MT_bgd";

const string file_SIG_mu = "Wprime_analysis_MuMET.root";
const string file_BGD_mu = "Wprime_bgd.root";
const string file_data_mu = "Wprime_data.root";

enum channel {wprime_MuMET = 0, wprime_ElMET};

#endif // #define _wprimeFitter_signalDescriptions_h_
