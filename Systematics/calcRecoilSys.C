//Usage: root -b -l -q calcRecoilSys.C+
#include "TROOT.h"
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string>
#include <stdlib.h>
#include "TFile.h"
#include "TTree.h" 
#include "TH1F.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include <fstream>

#include "../root_macros/common.h"

float NEvtsCorMET(TTree* tEvts, int seed);


void
calcRecoilSys(){

  //We use recoil method which is described in our AN-11-259 (chapter 6).
  //The code extracts gen level WZ vector from MC truth (the transverse
  //component), by traversing monte carlo truth particles and getting
  //components of W and Z
  // (this is done on tree level, so I am not sure how useful it is for you).

//Then it is filled to TVector2 containing the transversal components.
//Same is done for reconstructed W and Z
  ofstream fSigRes("SysSigMETRes.dat");
  ofstream fSigScale("SysSigMETScale.dat");
  ofstream fBkgRes("SysBkgMETRes.dat");
  ofstream fBkgScale("SysBkgMETScale.dat");
  fSigRes  << setiosflags(ios::fixed) << setprecision(4) << setiosflags(ios::left);
  fSigRes<<"Mass/F:SysErr/F"<<endl;
  fSigScale  << setiosflags(ios::fixed) << setprecision(4) << setiosflags(ios::left);
  fSigScale<<"Mass/F:SysErr/F"<<endl;
  fBkgRes  << setiosflags(ios::fixed) << setprecision(4) << setiosflags(ios::left);
  fBkgRes<<"Mass/F:SysErr/F"<<endl;
  fBkgScale  << setiosflags(ios::fixed) << setprecision(4) << setiosflags(ios::left);
  fBkgScale<<"Mass/F:SysErr/F"<<endl;

  vector<string> names;
  for(int mass=200; mass<=2000; mass+=100) names.push_back(Form("WprimeToWZTo3LNu_M-%i", mass));
  names.push_back("WZJetsTo3LNu");

  TFile *fIn = TFile::Open("../../../WprimeWZ.root", "read"); assert(fIn);

  for(int iname=0; iname<(int)names.size(); ++iname){
    const string & name = names[iname];

    ofstream & fRes   = (name.find("Wprime") == string::npos) ? fBkgRes : fSigRes;
    ofstream & fScale = (name.find("Wprime") == string::npos) ? fBkgScale : fSigScale;

    float mass = (name.find("M-") == string::npos) ? 0 : atof(name.substr(name.find("M-")+2).c_str());
    //cout<<"read mass of "<<mass<<" from "<<inputFile<<endl;

    //get number of base events 
    TTree* tMET = getTree(fIn, name, "tEvts_MET"); assert(tMET);
    float nEvtsUncorMET = 0;
    tMET->Draw("MET", "weight", "para goff");
    double n = tMET->GetSelectedRows(); 
    for(int ievt=0; ievt<n; ++ievt){
      float weight = tMET->GetW()[ievt];
      nEvtsUncorMET += weight;
    }
    cout<<"For "<<name<<" found "<<nEvtsUncorMET<<" events passing uncorrected met "<<endl;

    //Now count number of mod met events
    TTree* tValidW = getTree(fIn, name, "tEvts_ValidW"); assert(tValidW);

    //get number of scaled met events
    float nEvtsScaledMET = NEvtsCorMET(tValidW, 0);
    cout<<" Error taken on "<<name<<" due to met scale is "
        <<fabs(nEvtsScaledMET - nEvtsUncorMET) / nEvtsUncorMET
        <<endl;

    //get number of smeared met events
    TH1F hNEvts("hNEvts", "Yield due to various smearing seeds", 100, 0, 100.);
    hNEvts.StatOverflows(kTRUE);
    for(int seed=1; seed<=1000; ++seed){
      float nEvtsSmearedMET = NEvtsCorMET(tValidW, seed);
      hNEvts.Fill(nEvtsSmearedMET);
    }

    printf(" Error taken on %s due to met smearing is %.4f with mean = %.4f and gaus=%.4f\n", name.c_str(), hNEvts.GetRMS()/hNEvts.GetMean(), hNEvts.GetMean(), hNEvts.GetRMS());

    //Print Line for Latex Table
    printf("  %s & %.2f%% & %.2f%% \\\\\n", name.c_str(), fabs(nEvtsScaledMET - nEvtsUncorMET) / nEvtsUncorMET*100, hNEvts.GetRMS()/hNEvts.GetMean()*100);//Latex Line
    
    //print sys file for limits
    fScale<<mass<<"\t"<<fabs(nEvtsScaledMET - nEvtsUncorMET) / nEvtsUncorMET<<endl;
    fRes  <<mass<<"\t"<<hNEvts.GetRMS()/hNEvts.GetMean()<<endl;
    
  }
}

float
NEvtsCorMET(TTree* tEvts, int seed){
  //cout<<" Seed is "<<seed<<endl;
  tEvts->Draw("ZLep1PtGen:ZLep1PhiGen:ZLep2PtGen:ZLep2PhiGen:WLepPtGen:WLepPhiGen:WNeuPtGen:WNeuPhiGen:MET:METPhi", "weight", "para goff");
  double n = tEvts->GetSelectedRows(); 
  //if(1) cout<<"Found "<<n<<" input gen events "<<endl;
  float NEvts = 0;
  for(int ievt=0; ievt<n; ++ievt){
    int treeIdx = 0;
    const double ZLep1PtGen  = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep1PhiGen = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2PtGen  = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2PhiGen = tEvts->GetVal(treeIdx++)[ievt];
    const double WLepPtGen   = tEvts->GetVal(treeIdx++)[ievt];
    const double WLepPhiGen  = tEvts->GetVal(treeIdx++)[ievt];
    const double WNeuPtGen   = tEvts->GetVal(treeIdx++)[ievt];
    const double WNeuPhiGen  = tEvts->GetVal(treeIdx++)[ievt];
    const float weight       = tEvts->GetW()[ievt];

    if(ZLep1PtGen == 0 || 
       ZLep2PtGen == 0 ||
       WLepPtGen  == 0){
      cout<<"Skipping event b/c not 3 leptons"<<endl;
      continue;
    }
    
    if(WNeuPtGen  == 0){
      cout<<"Skipping event b/c no gen neutrino"<<endl;
      continue;
    }

    TVector2 zLep1, zLep2, wLep, met;
    zLep1.SetMagPhi(ZLep1PtGen, ZLep1PhiGen);
    zLep2.SetMagPhi(ZLep2PtGen, ZLep2PhiGen);
    wLep.SetMagPhi(WLepPtGen, WLepPhiGen);
    met.SetMagPhi(WNeuPtGen, WNeuPhiGen);

    TVector2 zrec( zLep1.Px() + zLep2.Px(), zLep1.Py() + zLep2.Py());
    TVector2 wrec(  wLep.Px() +   met.Px(),  wLep.Py() +   met.Py());
    TVector2 wz_gen = zrec + wrec;

    //TVector2 zrec( z.Px(), z.Py());
    //TVector2 wrec( w.Px(), w.Py());

    //Then ux and uy components of recoil are calculated:

    double ux= -wrec.Mod()*cos(wrec.Phi()-wz_gen.Phi())-zrec.Mod()*cos(zrec.Phi()-wz_gen.Phi());
    double uy= -wrec.Mod()*sin(wrec.Phi()-wz_gen.Phi())-zrec.Mod()*sin(zrec.Phi()-wz_gen.Phi());

    //cout<<"Ux="<<ux<<" and uy="<<uy<<endl;
    //double recoilScale = wz_gen->Mod(); //this is no longer necessary

    //Before we used values based on the recoil scale, but for Type1 pf MET,
    //we only use 0.5 GeV for resolution smearing (both of ux and uy
    //components),
    //and 2% for scale uncertainty (applied only to ux)

    double ux_corr = ux;
    double uy_corr = uy;

    if(seed>0){//smearing:
      TRandom random(seed);
      ux_corr += random.Gaus(0, 0.5);
      uy_corr += random.Gaus(0, 0.5);
    }else{//scale:
      ux_corr *= 0.98;
      uy_corr *= 0.98;
    }

    //After it is scaled or smeared, it is applied back on gen-level 2-vector:
  
    TVector2 ux_v2corr( ux_corr*cos(wz_gen.Phi()), ux_corr*sin(wz_gen.Phi()));
    TVector2 uy_v2corr(-uy_corr*sin(wz_gen.Phi()), uy_corr*cos(wz_gen.Phi()));

    //And we calculate smeared/shifted MET:

    double Met_Px_corr = -ux_v2corr.Px() - uy_v2corr.Px() - zrec.Px() - wLep.Px();
    double Met_Py_corr = -ux_v2corr.Py() - uy_v2corr.Py() - zrec.Py() - wLep.Py();

    double Met_modified = sqrt( pow(Met_Px_corr,2) + pow(Met_Py_corr,2) );

    if(Met_modified > 30.) NEvts += weight;

//  cout<<" Orig met: "<<MET
//      <<" Corr met: "<<Met_modified
//      <<endl;
  }//loop over nevts
  //cout<<" Number of events passing modified met is "<<NEvts<<endl;
  return NEvts;
/*
//For scaling, it is enough to calculate yield difference after cutting
//on modified and non-modified MET to get uncertainty effect (on either
//yield or acceptance).
  
//   For smearing, I calculate it 500 times for each event with different
//   random smearing, and keep 500 different sums,
//   then I look at RMS of yield distribution as the uncertainty on yield
//   (or acceptance).

//   Caveat is that at present we do this only on signal because it uses WZ
//                                                     to compute the recoil..
                                                    
//                                                     For PU systematic, I now use truth method. I think you can straight
//                                                     use the existing LumiReweighting package and 1D part, if the truth
//   number of vertices input is provided. I just forwarded you answer from
//   Javier and Jordi where they explained me how to compute data PU
//   distributions that are input to this (there is also twiki link in that
//                                         mail). I can send you those distributions for 2011 period as well.
  
//   Let me know if you need more information about this.. :)
*/  
}
