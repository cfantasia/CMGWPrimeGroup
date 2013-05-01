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

const int nch = 4;

void CalcMETSys(TFile* fIn, const vector<string> & names, const int mass, ofstream & fRes, ofstream & fScale);
void CalcMETSys(TFile* fIn, const string & name, const int mass, ofstream & fRes, ofstream & fScale){
  vector<string> names(1,name);
  CalcMETSys(fIn, names, mass, fRes, fScale);
}
Value NEvtsCorMET(TTree* tEvts, int seed, const string & moreCuts);


void
calcRecoilSys(){

  //We use recoil method which is described in our AN-11-259 (chapter 6).
  //The code extracts gen level WZ vector from MC truth (the transverse
  //component), by traversing monte carlo truth particles and getting
  //components of W and Z
  // (this is done on tree level, so I am not sure how useful it is for you).

//Then it is filled to TVector2 containing the transversal components.
//Same is done for reconstructed W and Z
  TFile *fIn = TFile::Open("../../../WprimeWZ.root", "read"); assert(fIn);

  ofstream fSigRes("SysSigMETRes.dat");
  ofstream fSigScale("SysSigMETScale.dat");
  fSigRes  << setiosflags(ios::fixed) << setprecision(4) << setiosflags(ios::left);
  fSigRes<<"Mass/F:SysErr/F"<<endl;
  fSigScale  << setiosflags(ios::fixed) << setprecision(4) << setiosflags(ios::left);
  fSigScale<<"Mass/F:SysErr/F"<<endl;
  for(int mass=200; mass<=2000; mass+=100){
    CalcMETSys(fIn, Form("WprimeToWZTo3LNu_M-%i", mass), mass, fSigRes, fSigScale);
  }

  ofstream fBkgRes("SysBkgMETRes.dat");
  ofstream fBkgScale("SysBkgMETScale.dat");
  fBkgRes  << setiosflags(ios::fixed) << setprecision(4) << setiosflags(ios::left);
  fBkgRes<<"Mass/F:SysErr/F"<<endl;
  fBkgScale  << setiosflags(ios::fixed) << setprecision(4) << setiosflags(ios::left);
  fBkgScale<<"Mass/F:SysErr/F"<<endl;
  for(int mass=200; mass<=2000; mass+=100){
    CalcMETSys(fIn, BkgSamples(), mass, fBkgRes, fBkgScale);
  }

}

void
CalcMETSys(TFile* fIn, const vector<string> & names, const int mass, ofstream & fRes, ofstream & fScale){
  string analyiscuts = AnalysisCuts(mass);
  const string name = (names[0] == "WZJetsTo3LNu") ? "Bkg" : names[0];

  cout<<name;
  fScale<<mass;
  fRes  <<mass;
  for(int ch=0; ch<nch; ++ch){
    string cuts = analyiscuts + Form(" && EvtType == %i", ch);
    //get number of base events
    TTree* tMET = getTree(fIn, names, "tEvts_MET"); assert(tMET);
    Value vEvtsUncorMET = GetNEvtsAndError(tMET, Form("weight*(%s)", cuts.c_str()));
    cout<<"For "<<name<<" found "<<vEvtsUncorMET<<" events passing uncorrected met "<<endl;

    //Now count number of mod met events
    TTree* tValidW = getTree(fIn, names, "tEvts_ValidW"); assert(tValidW);
    
    //get number of scaled met events
    Value vEvtsScaledMET = NEvtsCorMET(tValidW, 0, cuts);
    Value vDiffScaledMET = ShiftErr(vEvtsUncorMET, vEvtsScaledMET);
    cout<<" Error taken on "<<name<<" due to met scale is "
        <<vDiffScaledMET<<" with nevts="<<vEvtsScaledMET<<endl;
    
    //get number of smeared met events
    TH1F hNEvts("hNEvts", "Yield due to various smearing seeds", 100, 0, 100.);
    hNEvts.StatOverflows(kTRUE);
    for(int seed=1; seed<=1000; ++seed){
      Value vEvtsSmearedMET = NEvtsCorMET(tValidW, seed, cuts);
      hNEvts.Fill(vEvtsSmearedMET.val);
    }
    Value vEvtsSmearedMET(hNEvts.GetMean()+hNEvts.GetRMS(), hNEvts.GetRMSError());
    Value vDiffSmearedMET = ShiftErr(vEvtsUncorMET, vEvtsSmearedMET);
    printf(" Error taken on %s due to met smearing is %.4f with mean = %.4f and gaus=%.4f\n", name.c_str(), vDiffSmearedMET.val, hNEvts.GetMean(), hNEvts.GetRMS());

    //Print Line for Latex Table
    printf(" & %.2f%% & %.2f%% ", vDiffScaledMET.val*100, vDiffSmearedMET.val*100);//Latex Line
    
    //print sys file for limits
    fScale<<"\t"<<vDiffScaledMET.val<<"\t"<<vDiffScaledMET.err;
    fRes  <<"\t"<<vDiffSmearedMET.val<<"\t"<<vDiffSmearedMET.err;
  }//ch loop
  printf("\\\\\n");
  fScale<<endl;
  fRes  <<endl;
}

Value
NEvtsCorMET(TTree* tEvts, int seed, const string & moreCuts){
  TH1F hist("hist", "Dummy Hist", 1, 0, 10);
  hist.Sumw2();
  //cout<<" Seed is "<<seed<<endl;
  tEvts->Draw("ZLep1PtGen:ZLep1PhiGen:ZLep2PtGen:ZLep2PhiGen:WLepPtGen:WLepPhiGen:WNeuPtGen:WNeuPhiGen:ZLep1Pt:ZLep1Phi:ZLep2Pt:ZLep2Phi:WLepPt:WLepPhi:MET:METPhi", Form("weight*(%s)", moreCuts.c_str()), "para goff");
  double n = tEvts->GetSelectedRows(); 
  //if(1) cout<<"Found "<<n<<" input gen events "<<endl;
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
    
    const double ZLep1Pt  = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep1Phi = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2Pt  = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2Phi = tEvts->GetVal(treeIdx++)[ievt];
    const double WLepPt   = tEvts->GetVal(treeIdx++)[ievt];
    const double WLepPhi  = tEvts->GetVal(treeIdx++)[ievt];
    const double MET   = tEvts->GetVal(treeIdx++)[ievt];
    const double METPhi  = tEvts->GetVal(treeIdx++)[ievt];
    const float weight       = tEvts->GetW()[ievt];

    if(ZLep1PtGen == 0 || 
       ZLep2PtGen == 0 ||
       WLepPtGen  == 0){
      //cout<<"Skipping event b/c not 3 leptons"<<endl;
      continue;
    }
    
    if(WNeuPtGen  == 0){
      cout<<"Skipping event b/c no gen neutrino"<<endl;
      continue;
    }

    TVector2 zLep1, zLep2, wLep, met;
    zLep1.SetMagPhi(ZLep1Pt, ZLep1Phi);
    zLep2.SetMagPhi(ZLep2Pt, ZLep2Phi);
    wLep.SetMagPhi(WLepPt, WLepPhi);
    met.SetMagPhi(MET, METPhi);
    TVector2 zrec = zLep1 + zLep2;
    TVector2 wrec = wLep + met;


    TVector2 zLep1Gen, zLep2Gen, wLepGen, wNeuGen;
    zLep1Gen.SetMagPhi(ZLep1PtGen, ZLep1PhiGen);
    zLep2Gen.SetMagPhi(ZLep2PtGen, ZLep2PhiGen);
    wLepGen.SetMagPhi(WLepPtGen, WLepPhiGen);
    wNeuGen.SetMagPhi(WNeuPtGen, WNeuPhiGen);
    TVector2 wz_gen = zLep1Gen + zLep2Gen + wLepGen + wNeuGen;

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

    if(Met_modified > 30.) hist.Fill(1, weight);
    //if(MET > 30.) hist.Fill(1, weight);//Idiot check

    //printf(" Orig met: %.2f Corr met: %.2f\n", MET, Met_modified);

  }//loop over nevts
  //cout<<" Number of events passing modified met is "<<hist.GetBinContent(1)<<endl;
  return Value(hist.GetBinContent(1), hist.GetBinError(1));
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


