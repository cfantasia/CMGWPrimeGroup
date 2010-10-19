#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>
#include <map>
#include <cmath>
#include <vector>

#include "fit_wprime.h"

using std::string; using std::cout; using std::endl;

// ****** SETTINGS FOR ADVANCED USERS - BEGIN *********

// define RBW bin-size for convolving with RBW
// 5.0 GeV works fine for most fits
// In certain cases the plot appears with some weird resonance that goes away
// for a finer bin-size (e.g. 1.9 GeV). Haven't been able to debug this yet
const float rbw_bin_size = 5.0;
//const Double_t rbw_bin_size = 1.9;

// lower, upper limits for (non-smeared or MC-truth) muon Pt distributions
const float PTMIN = 100; const float PTMAX = 1500;

// ****** SETTINGS FOR ADVANCED USERS - END *********


const Double_t epsilon_dp = 0.5; // 0.01;
const Double_t epsilon_dm = 0.01; // 0.01;
const Double_t epsilon_dg = 0.01; // 0.01;


// cache here value returned by mySig (modulo some parameterization)
// given pt, mass and width, in order to speed up MINUIT minimization
struct rbw_point{ 
  unsigned p; // muon pt
  unsigned m; // Wprime mass
  unsigned g; // Wprime width fudge factor
  rbw_point(unsigned p_, unsigned m_, unsigned g_)
  {p = p_; m = m_; g = g_;}
};

struct lt_rbw_point
{
  bool operator()(const rbw_point & s1, const rbw_point &s2) const
  {
    return ( (s1.p < s2.p) || ((s1.p == s2.p) && (s1.m < s2.m)) 
	     || ((s1.p == s2.p) && (s1.m == s2.m) && (s1.g < s2.g)));
  }
};

static std::map<rbw_point, Double_t, lt_rbw_point> rbw_weight;
typedef std::map<rbw_point, Double_t, lt_rbw_point>::iterator It;

// must be defined in global scope because they are needed by static function
Double_t minRBW = 0;
Double_t maxRBW = 3000;

// to be set by caller, should be equal to histogram bin-size we fit
float BIN_SIZE = 0; 
void setBinSize(float bin_size){BIN_SIZE = bin_size;}
// this is the resolution function (histogram), also to be set by user
TH1F * gsmear; 
// set resolution function
void setResolution(TH1F * g){gsmear = g;}

bool landauFlag = false;
// set landau (true) or RBW (false) background
void setLandauBgd(bool flag){landauFlag = flag;}

// static functions; to be initialized in initFunc
// auxiliary RBW function for background
TF1 * rbw_aux_bgd = 0;
// auxiliary RBW function for signal
TF1 * rbw_aux = 0;
// this is the "auxiliary" laundau function; it will be used to calculate
// the integral between fXMIN, fXMAX, so we can normalize myLandau
TF1 * landau_aux = 0;

// fast inverse square root (from wikipedia!)
float InvSqrt(float x)
{
  union {
    float f;
    int i;
  } tmp;
  tmp.f = x;
  tmp.i = 0x5f3759df - (tmp.i >> 1);
  float y = tmp.f;
  return y * (1.5f - 0.5f * x * y * y);
}


const Double_t epsilon = 0.00001;
// this is the "auxiliary", non-normalized relativistic Breight Wigner function; 
// (ie. it is normalized for [0, infinity) only)
Double_t RBW_aux(Double_t * x, Double_t * par)
{
  // protect numerator (and integral of RBW)
  if(par[1] < epsilon)par[1] = epsilon;
  if(par[2] < epsilon)par[2] = epsilon;

  // par[0]: # of events (=1), par[1]: Mass, par[2]: Gamma
  Double_t ret = 0;
  Double_t part1 = x[0]*x[0] - par[1]*par[1];
  Double_t part2 = par[1]*par[2];
  Double_t denom = part1*part1 + part2*part2;
  if (denom != 0)
    ret = 2*par[0]*par[1]*part2*TMath::InvPi()/denom;
  return ret;
}

// if dN/dcos(theta) = (1 + costheta^2) or (1 +- costheta)^2
// use this to normalize part-(a) of function mySig below
// x[0] -> muon-pt , par[1] -> W' energy at CM
Double_t mySig_aux(Double_t * x, Double_t * par)
{
  Double_t ret = 0;
  // par[0]: scale, par[1]: energy
  if(par[1] > 0)
    {
      Double_t invpar1 = 1.0/par[1];
      Double_t arg = 2*x[0]*invpar1; // factor of 2 is because M_W = 2*E_mu
      if(TMath::Abs(arg) <= 1)
	{
	  Double_t sinsq = arg*arg; Double_t arg2 = 1 - sinsq;
	  if(arg2 > 0)
	    // (1.5/par[1]) is normalization factor so that integral is 1
	    ret = 1.5*par[0]*arg*(2 - sinsq)*invpar1*InvSqrt(arg2);
	}
    }
  return ret;
}

// relativistic Breit Wigner (normalized) function;
// to be used for description of background (between fXMIN, fXMAX)
Double_t myRBW(Double_t * x, Double_t * par)
{
  // par[0]: # of events, par[1]: M, par[2]: Gamma
  for(unsigned i = 0; i != 3; ++i)
    if(par[i] < 0)return 0;

  static TF1 rbw_aux_bgd_tmp("rbw_aux_bgd",RBW_aux,fXMIN,fXMAX, 3);
  rbw_aux_bgd = &rbw_aux_bgd_tmp;
  rbw_aux_bgd->SetParameters(1, par[1], par[2]);
  Double_t Integ = rbw_aux_bgd->Integral(fXMIN, fXMAX);
  Double_t par_aux[3] = {1, par[1], par[2]};
  // par[0] gives the # of events between fXMIN, fXMAX
  return par[0]*BIN_SIZE*rbw_aux_bgd->EvalPar(x, par_aux)/Integ;
}

static TF1 * mysig_aux = new TF1("mysig_aux", mySig_aux, PTMIN, PTMAX, 2);

// convolution of (a) dN/dcos(theta) = (1 + costheta^2) [or (1 +- costheta)^2]
// and (b) relativistic Breit Wigner
// this is considered unnormalized, except if it is integrated between 0, infinity
Double_t mySig(Double_t * x, Double_t * par)
{
  // par[0]: scale, par[1]: Mass, par[2]: "fudge" factor to increase width
  for(unsigned i = 0; i != 3; ++i)
    if(par[i] < 0)return 0;

  // Gamma is determined from Mass via width(W') = (4/3)* width(W) * M_w'/M_w
  // par[2] is the "fudge" factor
  Double_t width = par[2]*(4./3.) * width_W * (par[1]/mass_W);
  //  width = par[2]; if(par[2] < 0)return 0;

#if 0
  float delta = std::sqrt(float(30*par[1]*width));
  minRBW = TMath::Max(0.0, par[1] - delta);
  maxRBW = par[1]+delta;
  //  cout << " minRBW = " << minRBW << " maxRBW =  " << maxRBW << endl;
#endif

  unsigned Nbins_rbw = int((maxRBW - minRBW)/rbw_bin_size);
  static TF1 rbw_aux_tmp("rbw_aux", RBW_aux, minRBW, maxRBW, 3);
  rbw_aux = &rbw_aux_tmp;
  rbw_aux->SetParameters(1, par[1], width);
  Double_t par_aux[3] = {1, par[1], width};

  Double_t ret = 0;

  int x1 = int(x[0]/epsilon_dp); int y1 = int(par[1]/epsilon_dm); 
  int z1 = int(par[2]/epsilon_dg);
  rbw_point a(x1, y1, z1);
  It it = rbw_weight.find(a);
  if(it == rbw_weight.end())
    {
      Double_t Energy = minRBW - 0.5*rbw_bin_size;
      // loop over Energy bins
      for(unsigned i = 1; i <= Nbins_rbw; ++i)
	{
	  // find Energy at middle of bin
	  Energy += rbw_bin_size;
	  
	  // weight contribution according to RBW function
	  Double_t weight =rbw_bin_size*rbw_aux->EvalPar(&Energy, par_aux);
	  // weight contribution according to Jacobian-edge-like function
	  Double_t par_aux2[2] = {1, Energy};
	  ret += mysig_aux->EvalPar(x, par_aux2)*weight;
	}
      rbw_weight[a] = ret;

    }
  else
    ret = it->second;
  
  return par[0]*ret*BIN_SIZE;
}

// normalized Landau function
Double_t myLandau(Double_t * x, Double_t * par)
{
  // par[0]: # of events, par[1]: most probably value, par[2]: sigma

  // auxiliary landau function used to calculate integral and normalize myLandau
  static TF1 landau_aux_tmp("landau_aux", "landau", fXMIN, fXMAX);
  landau_aux = &landau_aux_tmp;
  landau_aux->SetParameters(1, par[1], par[2]);
  Double_t Integ = landau_aux->Integral(fXMIN, fXMAX);
  
  // par[0] gives the # of events between fXMIN, fXMAX
  return par[0]*BIN_SIZE*TMath::Landau(x[0], par[1], par[2])/Integ;

  // par[0] gives the # of events in full range (0, inf)
  //  return par[0]*BIN_SIZE*TMath::Landau(x[0], par[1], par[2])/par[2];
}

Double_t myExp(Double_t * x, Double_t * par)
{
  // par[0]: # of events, par[1]: lifetime = 1/Gamma
  if(par[1] < 0) par[1] = -par[1];

  // calculate integral between fXMIN, fXMAX
  Double_t Integ = par[1] * (TMath::Exp(-fXMIN/par[1]) - 
			     TMath::Exp(-fXMAX/par[1]));
  // par[0] gives the # of events between fXMIN, fXMAX
  return par[0]*BIN_SIZE*TMath::Exp(-x[0]/par[1])/Integ;

  // par[0] gives the # of events in full range (0, inf)
  //  return par[0]*BIN_SIZE*TMath::Exp(-x[0]/par[1])/par[1];
}

Double_t my_gauss(Double_t * x, Double_t * par)
{
  // par[0]: # of events, par[1]: mean, par[2]: sigma
  Double_t arg = 0;
  if (par[2]<0) par[2]=-par[2];
  if (par[2] != 0) arg = (x[0] - par[1])/par[2]; 
  return par[0]*BIN_SIZE*TMath::Exp(-0.5*arg*arg)/
    (TMath::Sqrt(2*TMath::Pi())*par[2]);
}

Double_t my_double_gauss(Double_t * x, Double_t * par)
{
  Double_t par_1[3] = {par[0], par[1], par[2]};
  Double_t par_2[3] = {par[3], par[4], par[5]};
  return my_gauss(x, par_1) + my_gauss(x, par_2);
}

// function to smear (global variable, be careful when setting!)
TF1 * func2smear = 0;

/* smear function "func2smear" (global pointer); 
   need to set pointer before calling smear_func 
   (should be done by spefic smearing function, e.g. smeared_exp) */
Double_t smear_func(Double_t * x, Double_t * par)
{
  Double_t func_value;
  Int_t bin_no;
  Axis_t xval;
  Int_t total_weight = int(gsmear->Integral());
  Double_t sum = 0;
  // this is just the implementation of Integral[ g(x-x')*f(x')*dx']
  for(bin_no = 1; bin_no <= gsmear->GetNbinsX(); ++bin_no) 
    {
      xval = x[0] - gsmear->GetBinCenter(bin_no);
      if(!func2smear->IsInside(&xval))
	continue;
      func_value = func2smear->EvalPar(&xval,par)*
	gsmear->GetBinContent(bin_no)/total_weight;
      sum += func_value;
    }
  return sum;
}


Double_t smeared_sig(Double_t * x, Double_t * par)
{
  // par[0]: scale, par[1]: Mass
  static TF1 mysig_tmp("mysig_tmp", mySig, PTMIN, PTMAX, 3);
  mysig_tmp.SetParameters(par[0], par[1], par[2]);
  func2smear = &mysig_tmp;
  return (smear_func(x, par));
}


Double_t myBgd(Double_t * x, Double_t * par)
{
  if(landauFlag)
    return myLandau(x, par);
  else
    return myRBW(x, par);
}

Double_t mySigBgd(Double_t * x, Double_t * par)
{
  for(int i = 0; i !=5 ; ++i)
    if(par[i] < 0)return 0;

  if(par[0] < par[3])
    return 0;
  
  // upper limit on mass
  if(par[4]  > (2*fXMAX - 10))
    return 0;

  Double_t par_bgd[3] = {par[0] - par[3], par[1], par[2]};
  Double_t par_sig[3] = {par[3], par[4], par[5]};
  return myBgd(x, par_bgd) + smeared_sig(x, par_sig);
}

