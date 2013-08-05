//Usage: root -b -l -q 'FixWResolution.C+(2000)'
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

void
FixWResolution(int mass=2000.){
  TH1::StatOverflows(kTRUE);

  string sample = Form("WprimeToWZTo3LNu_M-%i", mass);
  string outName = Form("WprimeWZ%i", mass);

  TFile *fIn = TFile::Open("../../../2013-04-01-FullDataset/WprimeWZ.root");

  TFile *fout = TFile::Open("WprimeWZ-Resolution-WLep.root", "recreate"); assert(fout);

  string outfile("fixRes-WLep.txt");
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
      <<"WLepPtGen/F:WLepEtaGen/F:WLepPhiGen/F:"
      <<"WNeuPtGen/F:WNeuEtaGen/F:WNeuPhiGen/F:"
      <<"WLepPt/F:WLepEta/F::WLepPhi/F:WLepPtErr/F:"
      <<"WNeuPt/F:WNeuEta/F::WNeuPhi/F:WNeuPtErr/F:"
      <<"tryWLepPt/F:tryWNeuPt/F:"
      <<"bilWLepPt/F:bilWNeuPt/F:"
      <<"filWLepPt/F:filWNeuPt/F:"
      <<"halWLepPt/F:halWNeuPt/F:"
      <<"catWLepPt/F:catWNeuPt/F"
      <<endl;  
  
  TH1F* holdRes1 = new TH1F("holdRes1", "", 40, -200, 200);
  TH1F* holdRes2 = new TH1F("holdRes2", "", 40, -200, 200);
  TH1F* holdRes  = new TH1F("holdRes" , "", 40, -200, 200);
  TH2F* holdRes2D= new TH2F("holdRes2D","", 40, -200, 200, 20, 0, 2000);

  TH1F* htryRes1 = new TH1F("htryRes1", "", 40, -200, 200); htryRes1->SetLineColor(kGreen);
  TH1F* htryRes2 = new TH1F("htryRes2", "", 40, -200, 200); htryRes2->SetLineColor(kGreen);
  TH1F* htryRes  = new TH1F("htryRes" , "", 40, -200, 200); htryRes ->SetLineColor(kGreen);
  TH2F* htryRes2D= new TH2F("htryRes2D","", 40, -200, 200, 20, 0, 2000); htryRes2D ->SetLineColor(kGreen);

  TH1F* hbilRes1 = new TH1F("hbilRes1", "", 40, -200, 200); hbilRes1->SetLineColor(kOrange);
  TH1F* hbilRes2 = new TH1F("hbilRes2", "", 40, -200, 200); hbilRes2->SetLineColor(kOrange);
  TH1F* hbilRes  = new TH1F("hbilRes" , "", 40, -200, 200); hbilRes ->SetLineColor(kOrange);
  TH2F* hbilRes2D= new TH2F("hbilRes2D","", 40, -200, 200, 20, 0, 2000); hbilRes2D ->SetLineColor(kOrange);

  TH1F* hfilRes1 = new TH1F("hfilRes1", "", 40, -200, 200); hfilRes1->SetLineColor(kRed);
  TH1F* hfilRes2 = new TH1F("hfilRes2", "", 40, -200, 200); hfilRes2->SetLineColor(kRed);
  TH1F* hfilRes  = new TH1F("hfilRes" , "", 40, -200, 200); hfilRes ->SetLineColor(kRed);
  TH2F* hfilRes2D= new TH2F("hfilRes2D","", 40, -200, 200, 20, 0, 2000); hfilRes2D ->SetLineColor(kRed);

  TH1F* hhalRes1 = new TH1F("hhalRes1", "", 40, -200, 200); hhalRes1->SetLineColor(kBlue);
  TH1F* hhalRes2 = new TH1F("hhalRes2", "", 40, -200, 200); hhalRes2->SetLineColor(kBlue);
  TH1F* hhalRes  = new TH1F("hhalRes" , "", 40, -200, 200); hhalRes ->SetLineColor(kBlue);
  TH2F* hhalRes2D= new TH2F("hhalRes2D","", 40, -200, 200, 20, 0, 2000); hhalRes2D ->SetLineColor(kBlue);

  TH1F* hcatRes1 = new TH1F("hcatRes1", "", 40, -200, 200); hcatRes1->SetLineColor(kCyan);
  TH1F* hcatRes2 = new TH1F("hcatRes2", "", 40, -200, 200); hcatRes2->SetLineColor(kCyan);
  TH1F* hcatRes  = new TH1F("hcatRes" , "", 40, -200, 200); hcatRes ->SetLineColor(kCyan);
  TH2F* hcatRes2D= new TH2F("hcatRes2D","", 40, -200, 200, 20, 0, 2000); hcatRes2D ->SetLineColor(kCyan);

  TTree* tEvts = getTree(fIn, sample, "tEvts_MET"); assert(tEvts);
  string cuts = Form("weight*(max(WLepPtGen,WLepPtGen)>500)*(EvtType == 1 || EvtType == 1)");
  int n = tEvts->Draw("EvtType:ZLep1PtGen:ZLep1EtaGen:ZLep1PhiGen:ZLep2PtGen:ZLep2EtaGen:ZLep2PhiGen:ZLep1Pt:ZLep1PtErr:ZLep1Eta:ZLep1Phi:ZLep2Pt:ZLep2PtErr:ZLep2Eta:ZLep2Phi:WLepPtGen:WLepEtaGen:WLepPhiGen:WNeuPtGen:WNeuEtaGen:WNeuPhiGen:WLepPt:WLepPtErr:WLepEta:WLepPhi:MET:MET/METSig:asinh(METPz/MET):METPhi", 
                      cuts.c_str(), "para goff");
  printf("There are %i events total\n", n);
  for(int ievt=0; ievt<min(n,200); ++ievt){  //loop over events
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
    
    TLorentzVector l1_gen, l2_gen, sum_gen;
    l1_gen.SetPtEtaPhiM(   WLepPtGen,WLepEtaGen,WLepPhiGen,0.); 
    l2_gen.SetPtEtaPhiM(   WNeuPtGen,WNeuEtaGen,WNeuPhiGen,0.); 
    sum_gen = l1_gen + l2_gen;
    float mZ_gen = sum_gen.M();

    TLorentzVector l1_old, l2_old, sum_old;
    l1_old.SetPtEtaPhiM(   WLepPt,WLepEta,WLepPhi,0.); 
    l2_old.SetPtEtaPhiM(   WNeuPt,WNeuEta,WNeuPhi,0.); 
    sum_old = l1_old + l2_old;
    float mZ_old = sum_old.M();

    if(l1_gen.DeltaR(l1_old) > 0.1) continue;

    printf("  GEN: mz %.2f zpt %.2f zpz %.2f zp %.2f ze %.2f zeta %.2f zphi %.2f dr %.2f l1pt %.2f l2pt %.2f l1eta %.2f l2eta %.2f l1phi %.2f l2phi %.2f cosdphi %.4f\n", mZ_gen, sum_gen.Pt(),sum_gen.Pz(),sum_gen.P(),sum_gen.E(),sum_gen.Eta(), sum_gen.Phi(), l1_gen.DeltaR(l2_gen), l1_gen.Pt(),l2_gen.Pt(),l1_gen.Eta(),l2_gen.Eta(),l1_gen.Phi(),l2_gen.Phi(), (sum_gen.Perp2() - l1_gen.Perp2() - l2_gen.Perp2())/(2.*l1_gen.Perp()*l2_gen.Perp()));
    l1_gen.Print();
    l2_gen.Print();
    sum_gen.Print();
    printf("  OLD: mz %.2f zpt %.2f zpz %.2f zp %.2f ze %.2f zeta %.2f zphi %.2f dr %.2f l1pt %.2f l2pt %.2f l1eta %.2f l2eta %.2f l1phi %.2f l2phi %.2f l1pterr %.2f l2pterr %.2f cosdphi %.4f\n", mZ_old, sum_old.Pt(),sum_old.Pz(),sum_old.P(),sum_old.E(),sum_old.Eta(),sum_old.Phi(),l1_old.DeltaR(l2_old), l1_old.Pt(),l2_old.Pt(),l1_old.Eta(),l2_old.Eta(),l1_old.Phi(),l2_old.Phi(), WLepPtErr, WNeuPtErr, (sum_old.Perp2() - l1_old.Perp2() - l2_old.Perp2())/(2.*l1_old.Perp()*l2_old.Perp()));
    l1_old.Print();
    l2_old.Print();
    sum_old.Print();
    printf("Want: k1: %.2f k2: %.2f arctan(k1/k2) = %.2f \n", l1_gen.Pt() / l1_old.Pt(), l2_gen.Pt() / l2_old.Pt(), TMath::ATan( (l1_gen.Pt()/l1_old.Pt()) / (l2_gen.Pt()/l2_old.Pt()) ) );

    printf("prod: gen Z+l1 %.4f gen z+l2 %.4f old z+l1 %.4f old z+l2 %.4f \n", 
           (sum_gen - l1_gen).Mag2(), (sum_gen - l2_gen).Mag2(), (sum_old - l1_old).Mag2(), (sum_old - l2_old).Mag2());

    float bestProd = pow(80.398 / sum_old.M(), 2);
    /*
    float bestDiff = 999.;
    for(float k = 0.01; k<10.; k+=0.01){
      TLorentzVector temp1,temp2;
      temp1.SetPtEtaPhiM(k*WLepPt,WLepEta,WLepPhi,0.); 
      temp2.SetPtEtaPhiM(  WNeuPt,WNeuEta,WNeuPhi,0.); 
      TLorentzVector temp = temp1 + temp2;
      float diff = abs(temp.M() - 80.398);
      if(diff < bestDiff){
        bestProd = k;
        bestDiff = diff;
      }
    }
    printf("BestDiff = %.4f Best k1*k2 = %.2f and ideal = %.2f\n", bestDiff, bestProd, pow(80.398 / sum_old.M(), 2));
    */
    TLorentzVector l1_wnt, l2_wnt, sum_wnt;
    float wantk1 = l1_gen.Pt() / l1_old.Pt();
    float wantk2 = l2_gen.Pt() / l2_old.Pt();
    cout<<"Want was "<<wantk1<<" and "<<wantk2<<" with prod "<<wantk1*wantk2<<endl;
    l1_wnt.SetPtEtaPhiM(wantk1*WLepPt,WLepEta,WLepPhi,0.); 
    l2_wnt.SetPtEtaPhiM(wantk2*WNeuPt,WNeuEta,WNeuPhi,0.); 
    sum_wnt = l1_wnt + l2_wnt;
    Analyze(wantk1,wantk2,l1_gen, l2_gen, sum_gen, l1_old, l2_old, sum_old, l1_wnt, l2_wnt, sum_wnt, true);

    //Try Me!     //Assume softer muon is better measured
    float tryk1 = pow(bestProd, WLepPt / (WLepPt + WNeuPt)) ; //Weight by pt
    float tryk2 = pow(bestProd, WNeuPt / (WLepPt + WNeuPt)) ; //Weight by pt

    cout<<"Try was "<<tryk1<<" and "<<tryk2<<endl;
    TLorentzVector l1_try, l2_try, sum_try;
    l1_try.SetPtEtaPhiM(tryk1*WLepPt,WLepEta,WLepPhi,0.); 
    l2_try.SetPtEtaPhiM(tryk2*WNeuPt,WNeuEta,WNeuPhi,0.); 
    sum_try = l1_try + l2_try;
    Analyze(tryk1,tryk2,l1_gen, l2_gen, sum_gen, l1_old, l2_old, sum_old, l1_try, l2_try, sum_try, true);
    htryRes1->Fill(l1_try.Pt() - l1_gen.Pt());
    htryRes2->Fill(l2_try.Pt() - l2_gen.Pt());
    htryRes->Fill(l1_try.Pt() - l1_gen.Pt());
    htryRes->Fill(l2_try.Pt() - l2_gen.Pt());
    htryRes2D->Fill(l1_try.Pt() - l1_gen.Pt(), l1_gen.Pt());
    htryRes2D->Fill(l2_try.Pt() - l2_gen.Pt(), l2_gen.Pt());
    cout<<"Try Res1 is now "<<htryRes1->GetRMS()<<" +- "<<htryRes1->GetRMSError()<<endl;
    cout<<"Try Res2 is now "<<htryRes2->GetRMS()<<" +- "<<htryRes2->GetRMSError()<<endl;
    cout<<"Try Res  is now "<<htryRes ->GetRMS()<<" +- "<<htryRes ->GetRMSError()<<endl;

    //Bil Me! Weight by res
    float bilk1 = pow(bestProd, WLepPtErr / (WLepPtErr + WNeuPtErr)) ; //Weight by pt err
    float bilk2 = pow(bestProd, WNeuPtErr / (WLepPtErr + WNeuPtErr)) ; //Weight by pt err 

    cout<<"Bil was "<<bilk1<<" and "<<bilk2<<endl;
    TLorentzVector l1_bil, l2_bil, sum_bil;
    l1_bil.SetPtEtaPhiM(bilk1*WLepPt,WLepEta,WLepPhi,0.); 
    l2_bil.SetPtEtaPhiM(bilk2*WNeuPt,WNeuEta,WNeuPhi,0.); 
    sum_bil = l1_bil + l2_bil;
    Analyze(bilk1,bilk2,l1_gen, l2_gen, sum_gen, l1_old, l2_old, sum_old, l1_bil, l2_bil, sum_bil, true);
    hbilRes1->Fill(l1_bil.Pt() - l1_gen.Pt());
    hbilRes2->Fill(l2_bil.Pt() - l2_gen.Pt());
    hbilRes->Fill(l1_bil.Pt() - l1_gen.Pt());
    hbilRes->Fill(l2_bil.Pt() - l2_gen.Pt());
    hbilRes2D->Fill(l1_bil.Pt() - l1_gen.Pt(), l1_gen.Pt());
    hbilRes2D->Fill(l2_bil.Pt() - l2_gen.Pt(), l2_gen.Pt());
    cout<<"Bil Res1 is now "<<hbilRes1->GetRMS()<<" +- "<<hbilRes1->GetRMSError()<<endl;
    cout<<"Bil Res2 is now "<<hbilRes2->GetRMS()<<" +- "<<hbilRes2->GetRMSError()<<endl;
    cout<<"Bil Res  is now "<<hbilRes ->GetRMS()<<" +- "<<hbilRes ->GetRMSError()<<endl;



    //Fil Me!
    float filk1 = pow(bestProd, WLepPt / (WLepPt + WNeuPt)) ; //Weight by sigma
    float filk2 = pow(bestProd, WNeuPt / (WLepPt + WNeuPt)) ; //Weight by sigma
    float A = (WLepPt / WLepPtErr) / (WNeuPt/WNeuPtErr);
    float bestFil = 9e9;
    for(float k1 = 0.1; k1<10.; k1+=0.01){
      float k2 = bestProd / k1;
      float B = fabs((1-k2) / (1-k1));
      float C = fabs( 1. - k1 ) * WLepPt / WLepPtErr + fabs( 1. - k2 ) * WNeuPt / WNeuPtErr;
      //if( C > 10) continue;
      if(C < bestFil){
        //if(fabs( A - B ) < bestFil){
        bestFil = C;//fabs( A - B );
        //printf("fil: k1 %.2f k2 %.2f A %.2f B: %.2f C: %.2f diff %.4f\n", k1, k2, A, B, C, fabs( A - B ));
        filk1 = k1;
        filk2 = k2;
      }
    }

    TLorentzVector l1_fil, l2_fil, sum_fil;
    cout<<"Fil was "<<filk1<<" and "<<filk2<<endl;
    l1_fil.SetPtEtaPhiM(filk1*WLepPt,WLepEta,WLepPhi,0.); 
    l2_fil.SetPtEtaPhiM(filk2*WNeuPt,WNeuEta,WNeuPhi,0.); 
    sum_fil = l1_fil + l2_fil;
    Analyze(filk1,filk2,l1_gen, l2_gen, sum_gen, l1_old, l2_old, sum_old, l1_fil, l2_fil, sum_fil, true);
    hfilRes1->Fill(l1_fil.Pt() - l1_gen.Pt());
    hfilRes2->Fill(l2_fil.Pt() - l2_gen.Pt());
    hfilRes->Fill(l1_fil.Pt() - l1_gen.Pt());
    hfilRes->Fill(l2_fil.Pt() - l2_gen.Pt());
    hfilRes2D->Fill(l1_fil.Pt() - l1_gen.Pt(), l1_gen.Pt());
    hfilRes2D->Fill(l2_fil.Pt() - l2_gen.Pt(), l2_gen.Pt());
    cout<<"Fil Res1 is now "<<hfilRes1->GetRMS()<<" +- "<<hfilRes1->GetRMSError()<<endl;
    cout<<"Fil Res2 is now "<<hfilRes2->GetRMS()<<" +- "<<hfilRes2->GetRMSError()<<endl;
    cout<<"Fil Res  is now "<<hfilRes ->GetRMS()<<" +- "<<hfilRes ->GetRMSError()<<endl;


    //Hal Me!, weight by gaus prob of pt err
    float halk1 = 1.;
    float halk2 = 1.;
    float bestProb = -1;
    TF1 zprob ( "zprob","gaus",60, 120);  zprob.SetParameters(1,80.398, 1.);
    TF1 k1prob("k1prob","gaus",0,10); k1prob.SetParameters(1,1., WLepPtErr/WLepPt); 
    TF1 k2prob("k2prob","gaus",0,10); k2prob.SetParameters(1,1., WNeuPtErr/WNeuPt); 
    for(float k1 = 0.1; k1<5.; k1+=0.001){
      for(float k2 = 0.1; k2<5.; k2+=0.001){
        //TLorentzVector temp1,temp2;
        //temp1.SetPtEtaPhiM(k1*WLepPt,WLepEta,WLepPhi,0.); 
        //temp2.SetPtEtaPhiM(k2*WNeuPt,WNeuEta,WNeuPhi,0.); 
        //TLorentzVector temp = temp1 + temp2;
        float prob = zprob(sum_old.M()*sqrt(k1*k2)) * k1prob(k1) * k2prob(k2);
        if(prob > bestProb){
          //printf("prob is %.4f k1: %.4f k2: %.4f\n", bestProb, k1, k2);
          bestProb = prob;
          halk1 = k1;
          halk2 = k2;
        }
      }
    }

    TLorentzVector l1_hal, l2_hal, sum_hal;
    cout<<"Hal was "<<halk1<<" and "<<halk2<<endl;
    l1_hal.SetPtEtaPhiM(halk1*WLepPt,WLepEta,WLepPhi,0.); 
    l2_hal.SetPtEtaPhiM(halk2*WNeuPt,WNeuEta,WNeuPhi,0.); 
    sum_hal = l1_hal + l2_hal;
    Analyze(halk1,halk2,l1_gen, l2_gen, sum_gen, l1_old, l2_old, sum_old, l1_hal, l2_hal, sum_hal, true);
    hhalRes1->Fill(l1_hal.Pt() - l1_gen.Pt());
    hhalRes2->Fill(l2_hal.Pt() - l2_gen.Pt());
    hhalRes->Fill(l1_hal.Pt() - l1_gen.Pt());
    hhalRes->Fill(l2_hal.Pt() - l2_gen.Pt());
    hhalRes2D->Fill(l1_hal.Pt() - l1_gen.Pt(), l1_gen.Pt());
    hhalRes2D->Fill(l2_hal.Pt() - l2_gen.Pt(), l2_gen.Pt());
    cout<<"Hal Res1 is now "<<hhalRes1->GetRMS()<<" +- "<<hhalRes1->GetRMSError()<<endl;
    cout<<"Hal Res2 is now "<<hhalRes2->GetRMS()<<" +- "<<hhalRes2->GetRMSError()<<endl;
    cout<<"Hal Res  is now "<<hhalRes ->GetRMS()<<" +- "<<hhalRes ->GetRMSError()<<endl;


    //Cat Me!, minimize # sigma added in quad
    float catk1 = 1.;
    float catk2 = 1.;
    float minDist = 999.;
    float k1Sigma = WLepPtErr/WLepPt;
    float k2Sigma = WNeuPtErr/WNeuPt;
    for(float k1 = 0.1; k1<5.; k1+=0.001){
      {
        float k2 = bestProd / k1;
        float dist = pow((k1 - 1)/k1Sigma,2) + pow((k2-1)/k2Sigma,2);
        if(dist < minDist){
          //printf("prob is %.4f k1: %.4f k2: %.4f\n", bestProb, k1, k2);
          minDist = dist;
          catk1 = k1;
          catk2 = k2;
        }
      }
    }

    TLorentzVector l1_cat, l2_cat, sum_cat;
    cout<<"Cat was "<<catk1<<" and "<<catk2<<endl;
    l1_cat.SetPtEtaPhiM(catk1*WLepPt,WLepEta,WLepPhi,0.); 
    l2_cat.SetPtEtaPhiM(catk2*WNeuPt,WNeuEta,WNeuPhi,0.); 
    sum_cat = l1_cat + l2_cat;
    Analyze(catk1,catk2,l1_gen, l2_gen, sum_gen, l1_old, l2_old, sum_old, l1_cat, l2_cat, sum_cat, true);
    hcatRes1->Fill(l1_cat.Pt() - l1_gen.Pt());
    hcatRes2->Fill(l2_cat.Pt() - l2_gen.Pt());
    hcatRes->Fill(l1_cat.Pt() - l1_gen.Pt());
    hcatRes->Fill(l2_cat.Pt() - l2_gen.Pt());
    hcatRes2D->Fill(l1_cat.Pt() - l1_gen.Pt(), l1_gen.Pt());
    hcatRes2D->Fill(l2_cat.Pt() - l2_gen.Pt(), l2_gen.Pt());
    cout<<"Cat Res1 is now "<<hcatRes1->GetRMS()<<" +- "<<hcatRes1->GetRMSError()<<endl;
    cout<<"Cat Res2 is now "<<hcatRes2->GetRMS()<<" +- "<<hcatRes2->GetRMSError()<<endl;
    cout<<"Cat Res  is now "<<hcatRes ->GetRMS()<<" +- "<<hcatRes ->GetRMSError()<<endl;


    //which boost gives half eneregy?

    //Old
    holdRes1->Fill(l1_old.Pt() - l1_gen.Pt());
    holdRes2->Fill(l2_old.Pt() - l2_gen.Pt());
    holdRes->Fill(l1_old.Pt() - l1_gen.Pt());
    holdRes->Fill(l2_old.Pt() - l2_gen.Pt());
    holdRes2D->Fill(l1_old.Pt() - l1_gen.Pt(), l1_gen.Pt());
    holdRes2D->Fill(l2_old.Pt() - l2_gen.Pt(), l2_gen.Pt());
    cout<<"Old Res1 is now "<<holdRes1->GetRMS()<<" +- "<<holdRes1->GetRMSError()<<" after iteration "<<holdRes1->GetEntries()<<endl;
    cout<<"Old Res2 is now "<<holdRes2->GetRMS()<<" +- "<<holdRes2->GetRMSError()<<" after iteration "<<holdRes2->GetEntries()<<endl;
    cout<<"Old Res  is now "<<holdRes ->GetRMS()<<" +- "<<holdRes ->GetRMSError()<<" after iteration "<<holdRes ->GetEntries()<<endl;


    //Write to text file
    ftxt<<EvtType<<" "
        <<ZLep1PtGen<<" "<<ZLep1EtaGen<<" "<<ZLep1PhiGen<<" "
        <<ZLep2PtGen<<" "<<ZLep2EtaGen<<" "<<ZLep2PhiGen<<" "
        <<ZLep1Pt   <<" "<<ZLep1Eta   <<" "<<ZLep1Phi   <<" "<<ZLep1PtErr<<" "
        <<ZLep2Pt   <<" "<<ZLep2Eta   <<" "<<ZLep2Phi   <<" "<<ZLep2PtErr<<" "
        <<WLepPtGen<<" "<<WLepEtaGen<<" "<<WLepPhiGen<<" "
        <<WNeuPtGen<<" "<<WNeuEtaGen<<" "<<WNeuPhiGen<<" "
        <<WLepPt   <<" "<<WLepEta   <<" "<<WLepPhi   <<" "<<WLepPtErr<<" "
        <<WNeuPt   <<" "<<WNeuEta   <<" "<<WNeuPhi   <<" "<<WNeuPtErr<<" "
        <<tryk1*WLepPt<<" "<<tryk2*WNeuPt<<" "
        <<bilk1*WLepPt<<" "<<bilk2*WNeuPt<<" "
        <<filk1*WLepPt<<" "<<filk2*WNeuPt<<" "
        <<halk1*WLepPt<<" "<<halk2*WNeuPt<<" "
        <<catk1*WLepPt<<" "<<catk2*WNeuPt<<" "
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
  mg->Add(holdRes);
//  mg->Add(htryRes);
//  mg->Add(hbilRes);
//  mg->Add(hfilRes);
//  mg->Add(hcatRes);
  mg->Add(hhalRes);
  mg->Draw("e nostack");
  leg->Draw();
  latexLabel.DrawLatex(0.15, 0.85, Form("Old #mu=%.0f #sigma=%.0f",holdRes->GetMean(), holdRes->GetRMS()));
//  latexLabel.DrawLatex(0.15, 0.80, Form("Try #mu=%.0f #sigma=%.0f",htryRes->GetMean(), htryRes->GetRMS()));
//  latexLabel.DrawLatex(0.15, 0.80, Form("Bil #mu=%.0f #sigma=%.0f",hbilRes->GetMean(), hbilRes->GetRMS()));
//  latexLabel.DrawLatex(0.15, 0.75, Form("Fil #mu=%.0f #sigma=%.0f",hfilRes->GetMean(), hfilRes->GetRMS()));
//  latexLabel.DrawLatex(0.15, 0.75, Form("Cat #mu=%.0f #sigma=%.0f",hcatRes->GetMean(), hcatRes->GetRMS()));
  latexLabel.DrawLatex(0.15, 0.70, Form("Hal #mu=%.0f #sigma=%.0f",hhalRes->GetMean(), hhalRes->GetRMS()));
  c1->SaveAs("NewRes_Pt.png");

  ftxt.close();
  
  fout->cd();

  htryRes1->Write();
  htryRes2->Write();
  htryRes ->Write();
  htryRes2D->Write();

  hbilRes1->Write();
  hbilRes2->Write();
  hbilRes ->Write();
  hbilRes2D->Write();

  holdRes1->Write();
  holdRes2->Write();
  holdRes ->Write();
  holdRes2D->Write();

  hhalRes1->Write();
  hhalRes2->Write();
  hhalRes ->Write();
  hhalRes2D->Write();

  hfilRes1->Write();
  hfilRes2->Write();
  hfilRes ->Write();
  hfilRes2D->Write();
 
  hcatRes1->Write();
  hcatRes2->Write();
  hcatRes ->Write();
  hcatRes2D->Write();

  fout->Close();

}

void
Analyze(float k1, float k2, 
        TLorentzVector & l1_gen, TLorentzVector & l2_gen, TLorentzVector &sum_gen,
        TLorentzVector & l1_old, TLorentzVector & l2_old, TLorentzVector &sum_old,
        TLorentzVector & l1_new, TLorentzVector & l2_new, TLorentzVector &sum_new,
        bool print){
      float mZ_new = sum_new.M();

      //if(fabs(sum_new.P() - pt(l1_new.DeltaR(l2_new))) > 10 ||
      //fabs(mZ_new - 80.398                        ) > 1.) return 9e9;

      
      float prod = k1*k2;
      
      float dMass = sum_new.M() - sum_gen.M();
      float dPt   = sum_new.Pt() - sum_gen.Pt();

      float fn_pt = pt(l1_new.DeltaR(l2_new));
      //float fn_pt = pt(l1_new.DeltaPhi(l2_new));

      //g->AddPoint(ipoint++, dMass, dPt, k1);
//        if(fabs(mZ_new - 80.398) < 0.1 && fabs(dMass) < 5 && fabs(dPt) < 10 ){
      //if(print || (fabs(sum_new.Pt() - pt(l1_new.DeltaR(l2_new))) < 10 && fabs(mZ_new - 80.398) < 1.)){
      if(print){
        printf("k1: %.2f k2: %.2f prod:%.1f atan: %.2f dMass:%.1f dPt:%.1f ratio %.1f fn_pt %.1f \n", k1, k2, prod, TMath::ATan(k1/k2), dMass,dPt,sum_new.M() /  sum_old.M() * sum_old.Pt(), fn_pt);
        printf("  NEW: mz %.2f zpt %.2f zpz %.2f zp %.2f ze %.2f dr %.2f l1pt %.2f l2pt %.2f l1eta %.2f l2eta %.2f l1phi %.2f l2phi %.2f cosdphi %.4f\n", mZ_new, sum_new.Pt(),sum_new.Pz(),sum_new.P(),sum_new.E(),l1_new.DeltaR(l2_new), l1_new.Pt(),l2_new.Pt(),l1_new.Eta(),l2_new.Eta(),l1_new.Phi(),l2_new.Phi(), (sum_new.Perp2() - l1_new.Perp2() - l2_new.Perp2())/(2.*l1_new.Perp()*l2_new.Perp()));
      }
      //if(fabs(mZ_new - 80.398) < 0.1) break; //this k2 was good enough for this k1

      return;
      //return pow(sum_new.Pt() - fn_pt,2) + pow(50*(mZ_new - 80.398),2) + 100*(k1*k1 + k2*k2);
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




