//Usage: root -b -l -q 'FixResolution.C+(2000)'
#include "UserCode/CMGWPrimeGroup/root_macros/common.h"
#include "TPaveStats.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include <fstream>

void
Analyze(float k1, float k2, 
        TLorentzVector & l1_gen, TLorentzVector & l2_gen, TLorentzVector &sum_gen,
        TLorentzVector & l1_old, TLorentzVector & l2_old, TLorentzVector &sum_old,
        TLorentzVector & l1_new, TLorentzVector & l2_new, TLorentzVector &sum_new,
        bool print=false);
float pt(float dr);

void FixTry(TLorentzVector & lA, TLorentzVector & lB, double lAPtErr, double lBPtErr, double bestProd, double MassMother);
void FixBil(TLorentzVector & lA, TLorentzVector & lB, double lAPtErr, double lBPtErr, double bestProd, double MassMother);
void FixFil(TLorentzVector & lA, TLorentzVector & lB, double lAPtErr, double lBPtErr, double bestProd, double MassMother);
void FixHal(TLorentzVector & lA, TLorentzVector & lB, double lAPtErr, double lBPtErr, double bestProd, double MassMother);
void FixCat(TLorentzVector & lA, TLorentzVector & lB, double lAPtErr, double lBPtErr, double bestProd, double MassMother);

void ApplyScale(TLorentzVector & v, const double & scale);
void Summerize(const TLorentzVector & l1_gen, const TLorentzVector & l2_gen,
               const TLorentzVector & l1_new, const TLorentzVector & l2_new,
               TH1F* htryRes1, TH1F* htryRes2, TH1F* htryRes, TH2F* htryRes2D,
               const string & name);
void Print(const TLorentzVector & lA, const TLorentzVector & lB, const double & kA, const double & kB);
TLorentzVector
UpdateMET(const TLorentzVector & lA_old, const TLorentzVector & lB_old, 
          const TLorentzVector & lA_new, const TLorentzVector & lB_new, 
          const TLorentzVector & lOldMET);


void
FixResolution(int mass=2000){
  TH1::StatOverflows(kTRUE);

  string sample = Form("WprimeToWZTo3LNu_M-%i", mass);
  string outName = Form("WprimeWZ%i", mass);

  TFile *fIn = TFile::Open("../../../2013-04-01-FullDataset/WprimeWZ.root");

  TFile *fout = TFile::Open("WprimeWZ-Resolution-test.root", "recreate"); assert(fout);

  string outfile("fixRes.txt");
  ofstream ftxt(outfile.c_str());
  if(!ftxt) { 
    cout << "Cannot open file " << outfile << endl; 
    abort();
  } 
  ftxt  << setiosflags(ios::fixed) << setprecision(4) << setiosflags(ios::left);
    
  //Write to text file
  ftxt<<"EvtType/F:"
      <<"ZLep1PtGen/F:ZLep1EtaGen/F:ZLep1PhiGen/F:"
      <<"ZLep2PtGen/F:ZLep2EtaGen/F:ZLep2PhiGen/F:"
      <<"ZLep1Pt/F:ZLep1Eta/F:ZLep1Phi/F:ZLep1PtErr/F:"
      <<"ZLep2Pt/F:ZLep2Eta/F:ZLep2Phi/F:ZLep2PtErr/F:"
      <<"tryZLep1Pt/F:tryZLep2Pt/F:"
      <<"bilZLep1Pt/F:bilZLep2Pt/F:"
      <<"filZLep1Pt/F:filZLep2Pt/F:"
      <<"halZLep1Pt/F:halZLep2Pt/F:"
      <<"catZLep1Pt/F:catZLep2Pt/F:"
      <<"WLepPtGen/F:WLepEtaGen/F:WLepPhiGen/F:"
      <<"WNeuPtGen/F:WNeuEtaGen/F:WNeuPhiGen/F:"
      <<"WLepPt/F:WLepEta/F:WLepPhi/F:WLepPtErr/F:"
      <<"WNeuPt/F:WNeuEta/F:WNeuPhi/F:WNeuPtErr/F:"
      <<"tryWLepPt/F:tryWNeuPt/F:"
      <<"bilWLepPt/F:bilWNeuPt/F:"
      <<"filWLepPt/F:filWNeuPt/F:"
      <<"halWLepPt/F:halWNeuPt/F:"
      <<"catWLepPt/F:catWNeuPt/F"
      <<endl;  
  
  TH1F* holdRes1  = new TH1F("holdRes1", "", 40, -200, 200);
  TH1F* holdRes2  = new TH1F("holdRes2", "", 40, -200, 200);
  TH1F* holdResZ  = new TH1F("holdResZ" , "", 40, -200, 200);
  TH2F* holdResZ2D= new TH2F("holdResZ2D","", 40, -200, 200, 20, 0, 2000);

  TH1F* htryRes1  = new TH1F("htryRes1", "", 40, -200, 200);  htryRes1->SetLineColor(kGreen);
  TH1F* htryRes2  = new TH1F("htryRes2", "", 40, -200, 200);  htryRes2->SetLineColor(kGreen);
  TH1F* htryResZ  = new TH1F("htryResZ" , "", 40, -200, 200); htryResZ->SetLineColor(kGreen);
  TH2F* htryResZ2D= new TH2F("htryResZ2D","", 40, -200, 200, 20, 0, 2000); htryResZ2D ->SetLineColor(kGreen);

  TH1F* hbilRes1  = new TH1F("hbilRes1", "", 40, -200, 200);  hbilRes1->SetLineColor(kOrange);
  TH1F* hbilRes2  = new TH1F("hbilRes2", "", 40, -200, 200);  hbilRes2->SetLineColor(kOrange);
  TH1F* hbilResZ  = new TH1F("hbilResZ" , "", 40, -200, 200); hbilResZ->SetLineColor(kOrange);
  TH2F* hbilResZ2D= new TH2F("hbilResZ2D","", 40, -200, 200, 20, 0, 2000); hbilResZ2D ->SetLineColor(kOrange);

  TH1F* hfilRes1  = new TH1F("hfilRes1", "", 40, -200, 200);  hfilRes1->SetLineColor(kRed);
  TH1F* hfilRes2  = new TH1F("hfilRes2", "", 40, -200, 200);  hfilRes2->SetLineColor(kRed);
  TH1F* hfilResZ  = new TH1F("hfilResZ" , "", 40, -200, 200); hfilResZ->SetLineColor(kRed);
  TH2F* hfilResZ2D= new TH2F("hfilResZ2D","", 40, -200, 200, 20, 0, 2000); hfilResZ2D ->SetLineColor(kRed);

  TH1F* hhalRes1  = new TH1F("hhalRes1", "", 40, -200, 200);  hhalRes1->SetLineColor(kBlue);
  TH1F* hhalRes2  = new TH1F("hhalRes2", "", 40, -200, 200);  hhalRes2->SetLineColor(kBlue);
  TH1F* hhalResZ  = new TH1F("hhalResZ" , "", 40, -200, 200); hhalResZ->SetLineColor(kBlue);
  TH2F* hhalResZ2D= new TH2F("hhalResZ2D","", 40, -200, 200, 20, 0, 2000); hhalResZ2D ->SetLineColor(kBlue);

  TH1F* hcatRes1  = new TH1F("hcatRes1", "", 40, -200, 200);  hcatRes1->SetLineColor(kCyan);
  TH1F* hcatRes2  = new TH1F("hcatRes2", "", 40, -200, 200);  hcatRes2->SetLineColor(kCyan);
  TH1F* hcatResZ  = new TH1F("hcatResZ" , "", 40, -200, 200); hcatResZ->SetLineColor(kCyan);
  TH2F* hcatResZ2D= new TH2F("hcatResZ2D","", 40, -200, 200, 20, 0, 2000); hcatResZ2D ->SetLineColor(kCyan);

  //W Histos
  TH1F* holdRes3  = new TH1F("holdRes3", "", 40, -200, 200);
  TH1F* holdRes4  = new TH1F("holdRes4", "", 40, -200, 200);
  TH1F* holdResW  = new TH1F("holdResW" , "", 40, -200, 200);
  TH2F* holdResW2D= new TH2F("holdResW2D","", 40, -200, 200, 20, 0, 2000);

  TH1F* htryRes3  = new TH1F("htryRes3", "", 40, -200, 200); htryRes3->SetLineColor(kGreen);
  TH1F* htryRes4  = new TH1F("htryRes4", "", 40, -200, 200); htryRes4->SetLineColor(kGreen);
  TH1F* htryResW  = new TH1F("htryResW" , "", 40, -200, 200); htryResW ->SetLineColor(kGreen);
  TH2F* htryResW2D= new TH2F("htryResW2D","", 40, -200, 200, 20, 0, 2000); htryResW2D ->SetLineColor(kGreen);

  TH1F* hbilRes3 = new TH1F("hbilRes3", "", 40, -200, 200); hbilRes3->SetLineColor(kOrange);
  TH1F* hbilRes4 = new TH1F("hbilRes4", "", 40, -200, 200); hbilRes4->SetLineColor(kOrange);
  TH1F* hbilResW  = new TH1F("hbilResW" , "", 40, -200, 200); hbilResW ->SetLineColor(kOrange);
  TH2F* hbilResW2D= new TH2F("hbilResW2D","", 40, -200, 200, 20, 0, 2000); hbilResW2D ->SetLineColor(kOrange);

  TH1F* hfilRes3 = new TH1F("hfilRes3", "", 40, -200, 200); hfilRes3->SetLineColor(kRed);
  TH1F* hfilRes4 = new TH1F("hfilRes4", "", 40, -200, 200); hfilRes4->SetLineColor(kRed);
  TH1F* hfilResW  = new TH1F("hfilResW" , "", 40, -200, 200); hfilResW ->SetLineColor(kRed);
  TH2F* hfilResW2D= new TH2F("hfilResW2D","", 40, -200, 200, 20, 0, 2000); hfilResW2D ->SetLineColor(kRed);

  TH1F* hhalRes3 = new TH1F("hhalRes3", "", 40, -200, 200); hhalRes3->SetLineColor(kBlue);
  TH1F* hhalRes4 = new TH1F("hhalRes4", "", 40, -200, 200); hhalRes4->SetLineColor(kBlue);
  TH1F* hhalResW  = new TH1F("hhalResW" , "", 40, -200, 200); hhalResW ->SetLineColor(kBlue);
  TH2F* hhalResW2D= new TH2F("hhalResW2D","", 40, -200, 200, 20, 0, 2000); hhalResW2D ->SetLineColor(kBlue);

  TH1F* hcatRes3 = new TH1F("hcatRes3", "", 40, -200, 200); hcatRes3->SetLineColor(kCyan);
  TH1F* hcatRes4 = new TH1F("hcatRes4", "", 40, -200, 200); hcatRes4->SetLineColor(kCyan);
  TH1F* hcatResW  = new TH1F("hcatResW" , "", 40, -200, 200); hcatResW ->SetLineColor(kCyan);
  TH2F* hcatResW2D= new TH2F("hcatResW2D","", 40, -200, 200, 20, 0, 2000); hcatResW2D ->SetLineColor(kCyan);


  TTree* tEvts = getTree(fIn, sample, "tEvts_MET"); assert(tEvts);
  //string cuts = Form("weight*(max(ZLep1PtGen,ZLep2PtGen)>500)*(EvtType>=2)");
  string cuts = Form("weight*(max(ZLep1PtGen,ZLep2PtGen)>500)*(EvtType==2)*(WFlavorGen>0)*(ZFlavorGen>0)");
  int n = tEvts->Draw("EvtType:ZLep1PtGen:ZLep1EtaGen:ZLep1PhiGen:ZLep2PtGen:ZLep2EtaGen:ZLep2PhiGen:ZLep1Pt:ZLep1PtErr:ZLep1Eta:ZLep1Phi:ZLep2Pt:ZLep2PtErr:ZLep2Eta:ZLep2Phi:WLepPtGen:WLepEtaGen:WLepPhiGen:WNeuPtGen:WNeuEtaGen:WNeuPhiGen:WLepPt:WLepPtErr:WLepEta:WLepPhi:MET:MET/METSig:asinh(METPz/MET):METPhi", 
                        cuts.c_str(), "para goff");
  printf("There are %i events total\n", n);
  for(int ievt=0; ievt<min(n,500); ++ievt){  //loop over events
    printf("--------------------\n");
    int treeIdx = 0;
    const int EvtType  = round(tEvts->GetVal(treeIdx++)[ievt]);

    const double ZLep1PtGen  = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep1EtaGen = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep1PhiGen = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2PtGen  = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2EtaGen = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2PhiGen = tEvts->GetVal(treeIdx++)[ievt];
    
    const double ZLep1Pt  = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep1PtErr= tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep1Eta = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep1Phi = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2Pt  = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2PtErr= tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2Eta = tEvts->GetVal(treeIdx++)[ievt];
    const double ZLep2Phi = tEvts->GetVal(treeIdx++)[ievt];

    const double WLepPtGen  = tEvts->GetVal(treeIdx++)[ievt];
    const double WLepEtaGen = tEvts->GetVal(treeIdx++)[ievt];
    const double WLepPhiGen = tEvts->GetVal(treeIdx++)[ievt];
    const double WNeuPtGen  = tEvts->GetVal(treeIdx++)[ievt];
    const double WNeuEtaGen = tEvts->GetVal(treeIdx++)[ievt];
    const double WNeuPhiGen = tEvts->GetVal(treeIdx++)[ievt];
    
    const double WLepPt  = tEvts->GetVal(treeIdx++)[ievt];
    const double WLepPtErr= tEvts->GetVal(treeIdx++)[ievt];
    const double WLepEta = tEvts->GetVal(treeIdx++)[ievt];
    const double WLepPhi = tEvts->GetVal(treeIdx++)[ievt];
    const double WNeuPt  = tEvts->GetVal(treeIdx++)[ievt];
    const double WNeuPtErr= tEvts->GetVal(treeIdx++)[ievt];
    const double WNeuEta = tEvts->GetVal(treeIdx++)[ievt];
    const double WNeuPhi = tEvts->GetVal(treeIdx++)[ievt];

    const float weight       = tEvts->GetW()[ievt];
    
    //Gen
    TLorentzVector l1_gen, l2_gen, zsum_gen;
    if(ZLep1PtGen > ZLep2PtGen){
      l1_gen.SetPtEtaPhiM(   ZLep1PtGen,ZLep1EtaGen,ZLep1PhiGen,0.); 
      l2_gen.SetPtEtaPhiM(   ZLep2PtGen,ZLep2EtaGen,ZLep2PhiGen,0.); 
    }else{
      l2_gen.SetPtEtaPhiM(   ZLep1PtGen,ZLep1EtaGen,ZLep1PhiGen,0.); 
      l1_gen.SetPtEtaPhiM(   ZLep2PtGen,ZLep2EtaGen,ZLep2PhiGen,0.); 
    }
    zsum_gen = l1_gen + l2_gen;
    float mZ_gen = zsum_gen.M();

    //Old
    TLorentzVector l1_old, l2_old, zsum_old;
    l1_old.SetPtEtaPhiM(   ZLep1Pt,ZLep1Eta,ZLep1Phi,0.); 
    l2_old.SetPtEtaPhiM(   ZLep2Pt,ZLep2Eta,ZLep2Phi,0.); 
    zsum_old = l1_old + l2_old;
    printf(" lAErr %.2f lBErr %.2f\n", ZLep1PtErr, ZLep2PtErr);
    Summerize(l1_gen, l2_gen, l1_old, l2_old, holdRes1, holdRes2, holdResZ, holdResZ2D, "Old Z");
    float mZ_old = zsum_old.M();

    //Require good match
    if(l1_gen.DeltaR(l1_old) > 0.1 || l2_gen.DeltaR(l2_old) > 0.1) continue;

    //Ugly Print Outs
    printf("  GEN: mz %.2f zpt %.2f zpz %.2f zp %.2f ze %.2f zeta %.2f zphi %.2f dr %.2f l1pt %.2f l2pt %.2f l1eta %.2f l2eta %.2f l1phi %.2f l2phi %.2f cosdphi %.4f\n", mZ_gen, zsum_gen.Pt(),zsum_gen.Pz(),zsum_gen.P(),zsum_gen.E(),zsum_gen.Eta(), zsum_gen.Phi(), l1_gen.DeltaR(l2_gen), l1_gen.Pt(),l2_gen.Pt(),l1_gen.Eta(),l2_gen.Eta(),l1_gen.Phi(),l2_gen.Phi(), (zsum_gen.Perp2() - l1_gen.Perp2() - l2_gen.Perp2())/(2.*l1_gen.Perp()*l2_gen.Perp()));
    l1_gen.Print();
    l2_gen.Print();
    zsum_gen.Print();
    printf("  OLD: mz %.2f zpt %.2f zpz %.2f zp %.2f ze %.2f zeta %.2f zphi %.2f dr %.2f l1pt %.2f l2pt %.2f l1eta %.2f l2eta %.2f l1phi %.2f l2phi %.2f l1pterr %.2f l2pterr %.2f cosdphi %.4f\n", mZ_old, zsum_old.Pt(),zsum_old.Pz(),zsum_old.P(),zsum_old.E(),zsum_old.Eta(),zsum_old.Phi(),l1_old.DeltaR(l2_old), l1_old.Pt(),l2_old.Pt(),l1_old.Eta(),l2_old.Eta(),l1_old.Phi(),l2_old.Phi(), ZLep1PtErr, ZLep2PtErr, (zsum_old.Perp2() - l1_old.Perp2() - l2_old.Perp2())/(2.*l1_old.Perp()*l2_old.Perp()));
    l1_old.Print();
    l2_old.Print();
    zsum_old.Print();

    //Basic Quantity
    float bestProd = pow(91.188 / zsum_old.M(), 2);

    //Want
    TLorentzVector l1_wnt, l2_wnt, zsum_wnt;
    float wantk1 = l1_gen.Pt() / l1_old.Pt();
    float wantk2 = l2_gen.Pt() / l2_old.Pt();
    cout<<"Z Want was "<<wantk1<<" and "<<wantk2<<" with prod "<<wantk1*wantk2<<endl;

    //Try Me!   
    TLorentzVector l1_try = l1_old;
    TLorentzVector l2_try = l2_old;
    FixTry(l1_try, l2_try, ZLep1PtErr, ZLep2PtErr, bestProd, 91.188);
    double tryk1 = l1_try.Pt() / l1_old.Pt();
    double tryk2 = l2_try.Pt() / l2_old.Pt();
    Print(l1_try, l2_try, tryk1, tryk2);
    Summerize(l1_gen, l2_gen, l1_try, l2_try, htryRes1, htryRes2, htryResZ, htryResZ2D, "Try Z");

    //Bil Me!
    TLorentzVector l1_bil = l1_old;
    TLorentzVector l2_bil = l2_old;
    FixBil(l1_bil, l2_bil, ZLep1PtErr, ZLep2PtErr, bestProd, 91.188);
    double bilk1 = l1_bil.Pt() / l1_old.Pt();
    double bilk2 = l2_bil.Pt() / l2_old.Pt();
    Print(l1_bil, l2_bil, bilk1, bilk2);
    Summerize(l1_gen, l2_gen, l1_bil, l2_bil, hbilRes1, hbilRes2, hbilResZ, hbilResZ2D, "Bil Z");

    //Fil Me!
    TLorentzVector l1_fil = l1_old;
    TLorentzVector l2_fil = l2_old;
    FixFil(l1_fil, l2_fil, ZLep1PtErr, ZLep2PtErr, bestProd, 91.188);
    double filk1 = l1_fil.Pt() / l1_old.Pt();
    double filk2 = l2_fil.Pt() / l2_old.Pt();
    Print(l1_fil, l2_fil, filk1, filk2);
    Summerize(l1_gen, l2_gen, l1_fil, l2_fil, hfilRes1, hfilRes2, hfilResZ, hfilResZ2D, "Fil Z");

    //Hal Me!, weight by gaus prob of pt err
    TLorentzVector l1_hal = l1_old;
    TLorentzVector l2_hal = l2_old;
    FixHal(l1_hal, l2_hal, ZLep1PtErr, ZLep2PtErr, bestProd, 91.188);
    double halk1 = l1_hal.Pt() / l1_old.Pt();
    double halk2 = l2_hal.Pt() / l2_old.Pt();
    Print(l1_hal, l2_hal, halk1, halk2);
    Summerize(l1_gen, l2_gen, l1_hal, l2_hal, hhalRes1, hhalRes2, hhalResZ, hhalResZ2D, "Hal Z");


    //Cat Me!
    TLorentzVector l1_cat = l1_old;
    TLorentzVector l2_cat = l2_old;
    FixCat(l1_cat, l2_cat, ZLep1PtErr, ZLep2PtErr, bestProd, 91.188);
    double catk1 = l1_cat.Pt() / l1_old.Pt();
    double catk2 = l2_cat.Pt() / l2_old.Pt();
    Print(l1_cat, l2_cat, catk1, catk2);
    Summerize(l1_gen, l2_gen, l1_cat, l2_cat, hcatRes1, hcatRes2, hcatResZ, hcatResZ2D, "Cat Z");

    ////////////////////
    // Now for the W ///
    ////////////////////
    printf("======== Working on the W now ========\n");

    TLorentzVector l3_gen, l4_gen, wsum_gen;
    l3_gen.SetPtEtaPhiM(   WLepPtGen,WLepEtaGen,WLepPhiGen,0.); 
    l4_gen.SetPtEtaPhiM(   WNeuPtGen,WNeuEtaGen,WNeuPhiGen,0.); 
    Print(l3_gen, l4_gen, 1., 1.);
    wsum_gen = l3_gen + l4_gen;
    float mW_gen = wsum_gen.M();
    //Summarize
    printf("  GEN: mB %.2f Bpt %.2f Bpz %.2f Bp %.2f Be %.2f Beta %.2f Bphi %.2f dr %.2f l3pt %.2f l4pt %.2f l3eta %.2f l4eta %.2f l3phi %.2f l4phi %.2f cosdphi %.4f\n", mW_gen, wsum_gen.Pt(),wsum_gen.Pz(),wsum_gen.P(),wsum_gen.E(),wsum_gen.Eta(), wsum_gen.Phi(), l3_gen.DeltaR(l4_gen), l3_gen.Pt(),l4_gen.Pt(),l3_gen.Eta(),l4_gen.Eta(),l3_gen.Phi(),l4_gen.Phi(), (wsum_gen.Perp2() - l3_gen.Perp2() - l4_gen.Perp2())/(2.*l3_gen.Perp()*l4_gen.Perp()));
    l3_gen.Print();
    l4_gen.Print();
    wsum_gen.Print();

    //Old
    TLorentzVector l3_old, l4_old, wsum_old;
    l3_old.SetPtEtaPhiM(   WLepPt,WLepEta,WLepPhi,0.); 
    l4_old.SetPtEtaPhiM(   WNeuPt,WNeuEta,WNeuPhi,0.); 
    Summerize(l3_gen, l4_gen, l3_old, l4_old, 
              holdRes3, holdRes4, holdResW, holdResW2D, "Old W");
    Print(l3_old, l4_old, 1., 1.);
    printf(" lAErr %.2f lBErr %.2f\n", WLepPtErr, WNeuPtErr);
    wsum_old = l3_old + l4_old;
    float mW_old = wsum_old.M();
    //Summerize
    printf("  OLD: mB %.2f Bpt %.2f Bpz %.2f Bp %.2f Be %.2f zBta %.2f Bphi %.2f dr %.2f l3pt %.2f l4pt %.2f l3eta %.2f l4eta %.2f l3phi %.2f l4phi %.2f l3pterr %.2f l4pterr %.2f cosdphi %.4f\n", mW_old, wsum_old.Pt(),wsum_old.Pz(),wsum_old.P(),wsum_old.E(),wsum_old.Eta(),wsum_old.Phi(),l3_old.DeltaR(l4_old), l3_old.Pt(),l4_old.Pt(),l3_old.Eta(),l4_old.Eta(),l3_old.Phi(),l4_old.Phi(), ZLep1PtErr, ZLep2PtErr, (wsum_old.Perp2() - l3_old.Perp2() - l4_old.Perp2())/(2.*l3_old.Perp()*l4_old.Perp()));
    l3_old.Print();
    l4_old.Print();
    wsum_old.Print();
    
    //Require good match
    if(l3_gen.DeltaR(l3_old) > 0.1) continue;


    //Basic Quantity
    bestProd = pow(80.398 / wsum_old.M(), 2);

    //Want
    TLorentzVector l3_wnt, l4_wnt, sumW_wnt;
    float wantk3 = l3_gen.Pt() / l3_old.Pt();
    float wantk4 = l4_gen.Pt() / l4_old.Pt();
    cout<<"W Want was "<<wantk3<<" and "<<wantk4<<" with prod "<<wantk3*wantk4<<endl;

    //Try
    TLorentzVector l3_try = l3_old;
    TLorentzVector l4_try = UpdateMET(l1_old, l2_old, l1_try, l2_try, l4_old);
    FixTry(l3_try, l4_try, WLepPtErr, WNeuPtErr, bestProd, 80.398);
    double tryk3 = l3_try.Pt() / l3_old.Pt();
    double tryk4 = l4_try.Pt() / l4_old.Pt();
    Summerize(l3_gen, l4_gen, l3_try, l4_try, 
              htryRes3, htryRes4, htryResW, htryResW2D, "Try W");
    Print(l3_try, l4_try, tryk3, tryk4);

    //Bil
    TLorentzVector l3_bil = l3_old;
    TLorentzVector l4_bil = UpdateMET(l1_old, l2_old, l1_bil, l2_bil, l4_old);
    FixBil(l3_bil, l4_bil, WLepPtErr, WNeuPtErr, bestProd, 80.398);
    double bilk3 = l3_bil.Pt() / l3_old.Pt();
    double bilk4 = l4_bil.Pt() / l4_old.Pt();
    Summerize(l3_gen, l4_gen, l3_bil, l4_bil, 
              hbilRes3, hbilRes4, hbilResW, hbilResW2D, "Bil W");
    Print(l3_bil, l4_bil, bilk3, bilk4);
    
    //Fil
    TLorentzVector l3_fil = l3_old;
    TLorentzVector l4_fil = UpdateMET(l1_old, l2_old, l1_fil, l2_fil, l4_old);
    FixFil(l3_fil, l4_fil, WLepPtErr, WNeuPtErr, bestProd, 80.398);
    double filk3 = l3_fil.Pt() / l3_old.Pt();
    double filk4 = l4_fil.Pt() / l4_old.Pt();
    Summerize(l3_gen, l4_gen, l3_fil, l4_fil, 
              hfilRes3, hfilRes4, hfilResW, hfilResW2D, "Fil W");
    Print(l3_fil, l4_fil, filk3, filk4);

    //Hal
    TLorentzVector l3_hal = l3_old;
    TLorentzVector l4_hal = UpdateMET(l1_old, l2_old, l1_hal, l2_hal, l4_old);
    FixHal(l3_hal, l4_hal, WLepPtErr, WNeuPtErr, bestProd, 80.398);
    double halk3 = l3_hal.Pt() / l3_old.Pt();
    double halk4 = l4_hal.Pt() / l4_old.Pt();
    Summerize(l3_gen, l4_gen, l3_hal, l4_hal, 
              hhalRes3, hhalRes4, hhalResW, hhalResW2D, "Hal W");
    Print(l3_hal, l4_hal, halk3, halk4);
    
    //Cat
    TLorentzVector l3_cat = l3_old;
    TLorentzVector l4_cat = UpdateMET(l1_old, l2_old, l1_cat, l2_cat, l4_old);
    FixCat(l3_cat, l4_cat, WLepPtErr, WNeuPtErr, bestProd, 80.398);
    double catk3 = l3_cat.Pt() / l3_old.Pt();
    double catk4 = l4_cat.Pt() / l4_old.Pt();
    Summerize(l3_gen, l4_gen, l3_cat, l4_cat, 
              hcatRes3, hcatRes4, hcatResW, hcatResW2D, "Cat W");
    Print(l3_cat, l4_cat, catk1, catk2);

    //Write to text file
    ftxt<<EvtType<<" "
        <<ZLep1PtGen<<" "<<ZLep1EtaGen<<" "<<ZLep1PhiGen<<" "
        <<ZLep2PtGen<<" "<<ZLep2EtaGen<<" "<<ZLep2PhiGen<<" "
        <<ZLep1Pt   <<" "<<ZLep1Eta   <<" "<<ZLep1Phi   <<" "<<ZLep1PtErr<<" "
        <<ZLep2Pt   <<" "<<ZLep2Eta   <<" "<<ZLep2Phi   <<" "<<ZLep2PtErr<<" "
        <<tryk1*ZLep1Pt<<" "<<tryk2*ZLep2Pt<<" "
        <<bilk1*ZLep1Pt<<" "<<bilk2*ZLep2Pt<<" "
        <<filk1*ZLep1Pt<<" "<<filk2*ZLep2Pt<<" "
        <<halk1*ZLep1Pt<<" "<<halk2*ZLep2Pt<<" "
        <<catk1*ZLep1Pt<<" "<<catk2*ZLep2Pt<<" "
        <<WLepPtGen<<" "<<WLepEtaGen<<" "<<WLepPhiGen<<" "
        <<WNeuPtGen<<" "<<WNeuEtaGen<<" "<<WNeuPhiGen<<" "
        <<WLepPt   <<" "<<WLepEta   <<" "<<WLepPhi   <<" "<<WLepPtErr<<" "
        <<WNeuPt   <<" "<<WNeuEta   <<" "<<WNeuPhi   <<" "<<WNeuPtErr<<" "
        <<tryk3*WLepPt<<" "<<tryk4*WNeuPt<<" "
        <<bilk3*WLepPt<<" "<<bilk4*WNeuPt<<" "
        <<filk3*WLepPt<<" "<<filk4*WNeuPt<<" "
        <<halk3*WLepPt<<" "<<halk4*WNeuPt<<" "
        <<catk3*WLepPt<<" "<<catk4*WNeuPt<<" "
        <<endl;


  }//Evt loop


  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextSize(0.05);
  latexLabel.SetTextFont(42);
  TCanvas *c1 = new TCanvas("c1","Graph2D example",0,0,600,600);
  //c1->SetLogy();
  TLegend *leg = new TLegend(0.7, 0.43,0.9, 0.89,"");
  prepLegend(leg);
  leg->AddEntry(holdRes1, "Old", "le");
//  leg->AddEntry(htryRes1, "Try", "le");
//  leg->AddEntry(htryRes1, "Bil", "le");
//  leg->AddEntry(hfilRes1, "Fil", "le");
//  leg->AddEntry(hcatRes1, "Cat", "le");
  leg->AddEntry(hhalRes1, "Hal", "le");

  THStack *mg1 = new THStack("mg1", ";p_{T}^{reco}-p_{T}^{gen} (GeV);a.u");
  mg1->Add(holdRes1);
//  mg1->Add(htryRes1);
//  mg1->Add(hbilRes1);
//  mg1->Add(hfilRes1);
//  mg1->Add(hcatRes1);
  mg1->Add(hhalRes1);
  mg1->Draw("e nostack");
  leg->Draw();
  latexLabel.DrawLatex(0.15, 0.85, Form("Old #mu=%.0f #sigma=%.0f",holdRes1->GetMean(), holdRes1->GetRMS()));
//  latexLabel.DrawLatex(0.15, 0.80, Form("Try #mu=%.0f #sigma=%.0f",htryRes1->GetMean(), htryRes1->GetRMS()));
//  latexLabel.DrawLatex(0.15, 0.75, Form("Try #mu=%.0f #sigma=%.0f",hbilRes1->GetMean(), hbilRes1->GetRMS()));
//  latexLabel.DrawLatex(0.15, 0.70, Form("Fil #mu=%.0f #sigma=%.0f",hfilRes1->GetMean(), hfilRes1->GetRMS()));
//  latexLabel.DrawLatex(0.15, 0.70, Form("Cat #mu=%.0f #sigma=%.0f",hcatRes1->GetMean(), hcatRes1->GetRMS()));
  latexLabel.DrawLatex(0.15, 0.65, Form("Hal #mu=%.0f #sigma=%.0f",hhalRes1->GetMean(), hhalRes1->GetRMS()));

  c1->SaveAs("NewRes_Pt1.png");

  THStack *mg2 = new THStack("mg2", ";p_{T}^{reco}-p_{T}^{gen} (GeV);a.u");
  mg2->Add(holdRes2);
//  mg2->Add(htryRes2);
//  mg2->Add(hbilRes2);
//  mg2->Add(hfilRes2);
//  mg2->Add(hcatRes2);
  mg2->Add(hhalRes2);
  mg2->Draw("e nostack");
  leg->Draw();
  latexLabel.DrawLatex(0.15, 0.85, Form("Old #mu=%.0f #sigma=%.0f",holdRes2->GetMean(), holdRes2->GetRMS()));
//  latexLabel.DrawLatex(0.15, 0.80, Form("Try #mu=%.0f #sigma=%.0f",htryRes2->GetMean(), htryRes2->GetRMS()));
//  latexLabel.DrawLatex(0.15, 0.75, Form("Bil #mu=%.0f #sigma=%.0f",hbilRes2->GetMean(), hbilRes2->GetRMS()));
//  latexLabel.DrawLatex(0.15, 0.70, Form("Fil #mu=%.0f #sigma=%.0f",hfilRes2->GetMean(), hfilRes2->GetRMS()));
//  latexLabel.DrawLatex(0.15, 0.70, Form("Cat #mu=%.0f #sigma=%.0f",hcatRes2->GetMean(), hcatRes2->GetRMS()));
  latexLabel.DrawLatex(0.15, 0.65, Form("Hal #mu=%.0f #sigma=%.0f",hhalRes2->GetMean(), hhalRes2->GetRMS()));
  c1->SaveAs("NewRes_Pt2.png");

  THStack *mg = new THStack("mg", ";p_{T}^{reco}-p_{T}^{gen} (GeV);a.u");
  mg->Add(holdResZ);
//  mg->Add(htryResZ);
//  mg->Add(hbilResZ);
//  mg->Add(hfilResZ);
//  mg->Add(hcatResZ);
  mg->Add(hhalResZ);
  mg->Draw("e nostack");
  leg->Draw();
  latexLabel.DrawLatex(0.15, 0.85, Form("Old #mu=%.0f #sigma=%.0f",holdResZ->GetMean(), holdResZ->GetRMS()));
//  latexLabel.DrawLatex(0.15, 0.80, Form("Try #mu=%.0f #sigma=%.0f",htryResZ->GetMean(), htryResZ->GetRMS()));
//  latexLabel.DrawLatex(0.15, 0.80, Form("Bil #mu=%.0f #sigma=%.0f",hbilResZ->GetMean(), hbilResZ->GetRMS()));
//  latexLabel.DrawLatex(0.15, 0.75, Form("Fil #mu=%.0f #sigma=%.0f",hfilResZ->GetMean(), hfilResZ->GetRMS()));
//  latexLabel.DrawLatex(0.15, 0.75, Form("Cat #mu=%.0f #sigma=%.0f",hcatResZ->GetMean(), hcatResZ->GetRMS()));
  latexLabel.DrawLatex(0.15, 0.70, Form("Hal #mu=%.0f #sigma=%.0f",hhalResZ->GetMean(), hhalResZ->GetRMS()));
  c1->SaveAs("NewRes_Pt.png");

  ftxt.close();
  
  fout->cd();

  //Z Plots
  htryRes1->Write();
  htryRes2->Write();
  htryResZ ->Write();
  htryResZ2D->Write();

  hbilRes1->Write();
  hbilRes2->Write();
  hbilResZ ->Write();
  hbilResZ2D->Write();

  holdRes1->Write();
  holdRes2->Write();
  holdResZ ->Write();
  holdResZ2D->Write();

  hhalRes1->Write();
  hhalRes2->Write();
  hhalResZ ->Write();
  hhalResZ2D->Write();

  hfilRes1->Write();
  hfilRes2->Write();
  hfilResZ ->Write();
  hfilResZ2D->Write();
 
  hcatRes1->Write();
  hcatRes2->Write();
  hcatResZ ->Write();
  hcatResZ2D->Write();

  //W Plots
  htryRes3->Write();
  htryRes4->Write();
  htryResW ->Write();
  htryResW2D->Write();

  hbilRes3->Write();
  hbilRes4->Write();
  hbilResW ->Write();
  hbilResW2D->Write();

  holdRes3->Write();
  holdRes4->Write();
  holdResW ->Write();
  holdResW2D->Write();

  hhalRes3->Write();
  hhalRes4->Write();
  hhalResW ->Write();
  hhalResW2D->Write();

  hfilRes3->Write();
  hfilRes4->Write();
  hfilResW ->Write();
  hfilResW2D->Write();
 
  hcatRes3->Write();
  hcatRes4->Write();
  hcatResW ->Write();
  hcatResW2D->Write();

  fout->Close();

}

void
Analyze(float k1, float k2, 
        TLorentzVector & l1_gen, TLorentzVector & l2_gen, TLorentzVector &sum_gen,
        TLorentzVector & l1_old, TLorentzVector & l2_old, TLorentzVector &sum_old,
        TLorentzVector & l1_new, TLorentzVector & l2_new, TLorentzVector &sum_new,
        bool print){
      float mZ_new = sum_new.M();

      float prod = k1*k2;
      
      float dMass = sum_new.M() - sum_gen.M();
      float dPt   = sum_new.Pt() - sum_gen.Pt();

      float fn_pt = pt(l1_new.DeltaR(l2_new));
      //float fn_pt = pt(l1_new.DeltaPhi(l2_new));

      //g->AddPoint(ipoint++, dMass, dPt, k1);
//        if(fabs(mZ_new - 91.188) < 0.1 && fabs(dMass) < 5 && fabs(dPt) < 10 ){
      //if(print || (fabs(sum_new.Pt() - pt(l1_new.DeltaR(l2_new))) < 10 && fabs(mZ_new - 91.188) < 1.)){
      if(print){
        printf("k1: %.2f k2: %.2f prod:%.1f atan: %.2f dMass:%.1f dPt:%.1f ratio %.1f fn_pt %.1f \n", k1, k2, prod, TMath::ATan(k1/k2), dMass,dPt,sum_new.M() /  sum_old.M() * sum_old.Pt(), fn_pt);
        printf("  NEW: mz %.2f zpt %.2f zpz %.2f zp %.2f ze %.2f dr %.2f l1pt %.2f l2pt %.2f l1eta %.2f l2eta %.2f l1phi %.2f l2phi %.2f cosdphi %.4f\n", mZ_new, sum_new.Pt(),sum_new.Pz(),sum_new.P(),sum_new.E(),l1_new.DeltaR(l2_new), l1_new.Pt(),l2_new.Pt(),l1_new.Eta(),l2_new.Eta(),l1_new.Phi(),l2_new.Phi(), (sum_new.Perp2() - l1_new.Perp2() - l2_new.Perp2())/(2.*l1_new.Perp()*l2_new.Perp()));
      }
      //if(fabs(mZ_new - 91.188) < 0.1) break; //this k2 was good enough for this k1
 
      return;
      //return pow(sum_new.Pt() - fn_pt,2) + pow(50*(mZ_new - 91.188),2) + 100*(k1*k1 + k2*k2);
}
 
void
Print(const TLorentzVector & lA, const TLorentzVector & lB, const double & kA, const double & kB){
  TLorentzVector sum = lA + lB;
  float mB = sum.M();
  printf(" kA: %.3f kB: %.3f\n", kA, kB);
  lA.Print();
  lB.Print();
  
  printf("  mB %.2f Bpt %.2f Bpz %.2f Bp %.2f Be %.2f dr %.2f lApt %.2f lBpt %.2f lAeta %.2f lBeta %.2f lAphi %.2f lBphi %.2f\n", 
         mB, sum.Pt(),sum.Pz(),sum.P(),sum.E(),
         lA.DeltaR(lB), 
         lA.Pt(),lB.Pt(),
         lA.Eta(),lB.Eta(),
         lA.Phi(),lB.Phi());

}


float pt(float dr){

  //return 607. / (0.450 + dr);
  //return 172. / (0.048 + fabs(dr)) - 34.;

  return ( 0.108 + asinh(0.722 / (dr - 0.304))) / 0.0038;

  if( dr < 0.17 ) return  925;
  if( dr < 0.22 ) return  925;
  if( dr < 0.28 ) return  725;
  if( dr < 0.32 ) return  625;
  if( dr < 0.38 ) return  525;
  if( dr < 0.43 ) return  475;
  if( dr < 0.47 ) return  425;
  if( dr < 0.52 ) return  375;
  if( dr < 0.57 ) return  375;
  if( dr < 0.62 ) return  325;
  if( dr < 0.68 ) return  325;
  if( dr < 0.73 ) return  275;
  if( dr < 0.77 ) return  275;
  if( dr < 0.82 ) return  275;
  if( dr < 0.88 ) return  225;
  if( dr < 0.93 ) return  225;
  if( dr < 0.98 ) return  225;
  if( dr < 1.02 ) return  225;
  if( dr < 1.08 ) return  225;
  if( dr < 1.12 ) return  225;
  if( dr < 1.17 ) return  225;
  if( dr < 1.23 ) return  175;
  if( dr < 1.27 ) return  175;
  if( dr < 1.33 ) return  175;
  if( dr < 1.38 ) return  175;
  if( dr < 1.42 ) return  175;
  if( dr < 1.48 ) return  175;
  if( dr < 1.67 ) return  175;
  if( dr < 1.73 ) return  175;
  if( dr < 1.77 ) return  125;
  if( dr < 1.83 ) return  125;
  if( dr < 1.88 ) return  125;
  if( dr < 1.98 ) return  125;
  else            return  100;;
}

void
FixTry(TLorentzVector & lA, TLorentzVector & lB, double lAPtErr, double lBPtErr, double bestProd, double MassMother){
  //Weight by Pt
  float trykA = pow(bestProd, lA.Pt() / (lA.Pt() + lB.Pt()));
  float trykB = pow(bestProd, lB.Pt() / (lA.Pt() + lB.Pt()));
  
  ApplyScale(lA, trykA);
  ApplyScale(lB, trykB);
}

void
FixBil(TLorentzVector & lA, TLorentzVector & lB, double lAPtErr, double lBPtErr, double bestProd, double MassMother){
  //Weight by Pt Err
  float bilkA = pow(bestProd, lAPtErr / (lAPtErr + lBPtErr));
  float bilkB = pow(bestProd, lBPtErr / (lAPtErr + lBPtErr));
  
  ApplyScale(lA, bilkA);
  ApplyScale(lB, bilkB);
}

void
FixFil(TLorentzVector & lA, TLorentzVector & lB, double lAPtErr, double lBPtErr, double bestProd, double MassMother){
  //Count #sigma diff from 1., 1.
  float filkA = 1.;
  float filkB = 1.;

  float bestFil = 9e9;
  for(float kA = 0.2; kA<5.; kA+=0.001){
    float kB = bestProd / kA;
    float C = fabs( 1. - kA ) * lA.Pt() / lAPtErr + fabs( 1. - kB ) * lB.Pt() / lBPtErr;
    if(C < bestFil){
      bestFil = C;
      //printf("fil: kA %.2f kB %.2f A %.2f B: %.2f C: %.2f diff %.4f\n", kA, kB, A, B, C, fabs( A - B ));
      filkA = kA;
      filkB = kB;
    }
  }
  
  ApplyScale(lA, filkA);
  ApplyScale(lB, filkB);
}

void
FixHal(TLorentzVector & lA, TLorentzVector & lB, double lAPtErr, double lBPtErr, double bestProd, double MassMother){
  //Maximize 3 gaus (mZ,k1,k2)
  float halkA = 1.;
  float halkB = 1.;
  
  TLorentzVector sum_old = lA + lB;
  float bestProb = -1;
  TF1 zprob ( "zprob","gaus",60, 120);  zprob.SetParameters(1,MassMother, 1.);
  TF1 kAprob("kAprob","gaus",0,10); kAprob.SetParameters(1,1., lAPtErr/lA.Pt()); 
  TF1 kBprob("kBprob","gaus",0,10); kBprob.SetParameters(1,1., lBPtErr/lB.Pt()); 
  for(float kA = 0.1; kA<5.; kA+=0.001){
    for(float kB = 0.1; kB<5.; kB+=0.001){
      float prob = zprob(sum_old.M()*sqrt(kA*kB)) * kAprob(kA) * kBprob(kB);
      if(prob > bestProb){
        //printf("prob is %.4f kA: %.4f kB: %.4f\n", bestProb, kA, kB);
        bestProb = prob;
        halkA = kA;
        halkB = kB;
      }
    }
  }
  
  ApplyScale(lA, halkA);
  ApplyScale(lB, halkB);
}

void
FixCat(TLorentzVector & lA, TLorentzVector & lB, double lAPtErr, double lBPtErr, double bestProd, double MassMother){
  //Dist in sigma
  float catkA = 1.;
  float catkB = 1.;
  
  float minDist = 999.;
  float kASigma = lAPtErr/lA.Pt();
  float kBSigma = lBPtErr/lB.Pt();
  for(float kA = 0.1; kA<5.; kA+=0.001){
    float kB = bestProd / kA;
    float dist = pow((kA - 1)/kASigma,2) + pow((kB-1)/kBSigma,2);
    if(dist < minDist){
      //printf("prob is %.4f kA: %.4f kB: %.4f\n", bestProb, kA, kB);
      minDist = dist;
      catkA = kA;
      catkB = kB;
    }
  }

  ApplyScale(lA, catkA);
  ApplyScale(lB, catkB);
}

TLorentzVector
UpdateMET(const TLorentzVector & lA_old, const TLorentzVector & lB_old, 
          const TLorentzVector & lA_new, const TLorentzVector & lB_new, 
          const TLorentzVector & lOldMET){
  TLorentzVector lModMET = lOldMET;
  
  //Modify for Z
  lModMET.SetPx( lModMET.Px() + 
                 (lA_old.Px() - lA_new.Px()) + 
                 (lB_old.Px() - lB_new.Px()) );
  lModMET.SetPy( lModMET.Py() + 
                 (lA_old.Py() - lA_new.Py()) + 
                 (lB_old.Py() - lB_new.Py()) );
  printf("Updating MET: \n");
  lA_old.Print();
  lA_new.Print();
  lB_old.Print();
  lB_new.Print();
  lOldMET.Print();
  lModMET.Print();
  

  return lModMET;
  
}

void
ApplyScale(TLorentzVector & v, const double & scale){
  v.SetPtEtaPhiM(scale*v.Pt(),v.Eta(), v.Phi(), v.M());
}


void
Summerize(const TLorentzVector & l1_gen, const TLorentzVector & l2_gen,
          const TLorentzVector & l1_new, const TLorentzVector & l2_new,
          TH1F* hRes1, TH1F* hRes2, TH1F* hRes, TH2F* hRes2D,
          const string & name){
  cout<<"For "<<name<<": "<<endl;
  hRes1->Fill(l1_new.Pt() - l1_gen.Pt());
  hRes2->Fill(l2_new.Pt() - l2_gen.Pt());
  hRes->Fill(l1_new.Pt() - l1_gen.Pt());
  hRes->Fill(l2_new.Pt() - l2_gen.Pt());
  hRes2D->Fill(l1_new.Pt() - l1_gen.Pt(), l1_gen.Pt());
  hRes2D->Fill(l2_new.Pt() - l2_gen.Pt(), l2_gen.Pt());
  cout<<" Res1 is now "<<hRes1->GetRMS()<<" +- "<<hRes1->GetRMSError()<<" after iteration "<<hRes1->GetEntries()<<endl;
  cout<<" Res2 is now "<<hRes2->GetRMS()<<" +- "<<hRes2->GetRMSError()<<" after iteration "<<hRes2->GetEntries()<<endl;
  cout<<" Res  is now "<<hRes ->GetRMS()<<" +- "<<hRes ->GetRMSError()<<" after iteration "<<hRes ->GetEntries()<<endl;

}
