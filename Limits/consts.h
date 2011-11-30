#include "TROOT.h"
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string>
#include <map>
#include <fstream>
#include <stdlib.h>
#include <algorithm>
#include "TFile.h"
#include "TTree.h" 
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TMath.h"

#include "../root_macros/common.h"
const float sLumiFrac = 0.045;

vector<string>
SampleName(int code){
  vector<string> name;
  switch(code){
  case 3:
    name.push_back( "WprimeToWZTo3LNu_M-300"); break;
  case 4:
    name.push_back( "WprimeToWZTo3LNu_M-400"); break;
  case 5:
    name.push_back( "WprimeToWZTo3LNu_M-500"); break;
  case 6:
    name.push_back( "WprimeToWZTo3LNu_M-600"); break;
  case 7:
    name.push_back( "WprimeToWZTo3LNu_M-700"); break;
  case 8:
    name.push_back( "WprimeToWZTo3LNu_M-800"); break;
  case 9:
    name.push_back( "WprimeToWZTo3LNu_M-900"); break;
  case 10:
    name.push_back( "WprimeToWZTo3LNu_M-1000"); break;
  case 11:
    name.push_back( "WprimeToWZTo3LNu_M-1100"); break;
  case 12:
    name.push_back( "WprimeToWZTo3LNu_M-1200"); break;
  case 13:
    name.push_back( "WprimeToWZTo3LNu_M-1300"); break;
  case 14:
    name.push_back( "WprimeToWZTo3LNu_M-1400"); break;
  case 15:
    name.push_back( "WprimeToWZTo3LNu_M-1500"); break;
  case 103:
    name.push_back( "TC_WZ_300"); break;
  case 104:
    name.push_back( "TC_WZ_400"); break;
  case 105:
    name.push_back( "TC_WZ_500"); break;
  case 207:
    name.push_back( "Summer11_RSZZmmjj_750"); 
    name.push_back( "Summer11_RSZZeejj_750"); 
    break;
  case 210:
    name.push_back( "Summer11_RSZZmmjj_1000"); 
    name.push_back( "Summer11_RSZZeejj_1000"); 
    break;
  case 212:
    name.push_back( "Summer11_RSZZmmjj_1250"); 
    name.push_back( "Summer11_RSZZeejj_1250"); 
    break;
  case 215:
    name.push_back( "Summer11_RSZZmmjj_1500"); 
    name.push_back( "Summer11_RSZZeejj_1500"); 
    break;
  case 217:
    name.push_back( "Summer11_RSZZmmjj_1750"); 
    name.push_back( "Summer11_RSZZeejj_1750"); 
    break;
  case 220:
    name.push_back( "Summer11_RSZZmmjj_2000"); 
    name.push_back( "Summer11_RSZZeejj_2000"); 
    break;
  default:
    cout<<"Failed looking for code "<<code<<endl;
    abort();
    break;
  }
  return name;

}

float
SysErr(string sample){
  if(!sample.find("GV")) return 0.13;
  if(!sample.find("ZZ")) return 0.075;
  if(!sample.find("WZ")) return 0.17;
  if(!sample.find("WprimeToWZTo3LNu_M-300")) return 0.0447;
  if(!sample.find("WprimeToWZTo3LNu_M-400")) return 0.0469;
  if(!sample.find("WprimeToWZTo3LNu_M-500")) return 0.0495;
  if(!sample.find("WprimeToWZTo3LNu_M-600")) return 0.0530;
  if(!sample.find("WprimeToWZTo3LNu_M-700")) return 0.0566;
  if(!sample.find("WprimeToWZTo3LNu_M-800")) return 0.0598;
  if(!sample.find("WprimeToWZTo3LNu_M-900")) return 0.0635;
  if(!sample.find("TC_WZ_300"))  return 0.0447;
  if(!sample.find("TC_WZ_400"))  return 0.0469;
  if(!sample.find("TC_WZ_500"))  return 0.0495;
  if(!sample.find("WprimeToWZTo3LNu_M-")) return 0.0635;//Assume >900 is same as 900
  return 0;
}

float
XSecWprimeWZ(float m){
  TTree* tXsec = new TTree("tXsec", "W' Cross Sections");
  tXsec->ReadFile("xSec_WZ.dat");
  tXsec->Draw("Xsec", Form("Mass==%f", m), "goff");
  
  assert(tXsec->GetSelectedRows() == 1);
  float retval = tXsec->GetV1()[0];
  delete tXsec;
  return retval;
}

float
XSecTCWZ(float rho, float pi){
  TTree* tXsec = new TTree("tXsec", "TC Cross Sections");
  tXsec->ReadFile("xSec_TCWZ.dat");
  tXsec->Draw("Xsec", Form("Rho==%f && Pi==%f", 
                           rho, pi), "goff");
  
  //std::cout<<Form("Rho==%f && Pi==%f", rho, pi)<<std::endl;
  assert(tXsec->GetSelectedRows() == 1);
  float retval = tXsec->GetV1()[0];
  delete tXsec;
  return retval;
}

float
MassPi(float rho){
  return 3*rho/4 - 25.;
}

float
XSecRSZZ(float m){
  TTree* tXsec = new TTree("tXsec", "RS Cross Sections");
  tXsec->ReadFile("xSec_RSZZ.dat");
  tXsec->Draw("Xsec", Form("Mass==%f", m), "goff");
  
  assert(tXsec->GetSelectedRows() == 1);
  float retval = tXsec->GetV1()[0];
  delete tXsec;
  return retval;
}

float
XSec(std::string sample, float mass){
  if(!sample.find("WprimeToWZTo3LNu_M-")) return XSecWprimeWZ(mass);
  if(!sample.find("TC_WZ_")) return XSecTCWZ(mass, MassPi(mass));
  if(sample.find("RSZZ") != string::npos)  return XSecRSZZ(mass);

  std::cout<<"Didn't find sample named "<<sample<<std::endl;
  abort();
  return 0;
}

float
AddInQuad(float a, float b){
  return sqrt(a*a + b*b);
}

const float FitWind_low = 70;
const float FitWind_high = 110;

const float ZWind_low = 80;
const float ZWind_high = 100;

Bool_t reject;
Double_t fline(Double_t *x, Double_t *par)
{
    if (reject && x[0] > ZWind_low && x[0] < ZWind_high) {
      TF1::RejectPoint();
      return 0;
   }
   return par[0] + par[1]*x[0];
}

float Slope(const float x1, const float y1,
            const float x2, const float y2){
  return (y2-y1)/(x2-x1);
}
float Intercept(const float x1, const float y1,
                const float x2, const float y2){
  return y1 - x1*(y2-y1)/(x2-x1);
}

float
IntersectionX(const float s1, const float i1,
              const float s2, const float i2){
  return (i2-i1)/(s1-s2);
}
