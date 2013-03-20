//Usage: root -b -l -q 'FixResolution.C+(2000)'
#include "UserCode/CMGWPrimeGroup/root_macros/common.h"
#include "TPaveStats.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLorentzVector.h"

float
Analyze(float k1, float k2, 
        TLorentzVector & l1_gen, TLorentzVector & l2_gen, TLorentzVector &sum_gen,
        TLorentzVector & l1_old, TLorentzVector & l2_old, TLorentzVector &sum_old,
        TLorentzVector & l1_new, TLorentzVector & l2_new, TLorentzVector &sum_new,
        bool print=false);
float pt(float dr);

void
FixResolution(int mass){
  TH1::StatOverflows(kTRUE);

  string sample = Form("WprimeToWZTo3LNu_M-%i", mass);
  string outName = Form("WprimeWZ%i", mass);

  TFile *fIn = TFile::Open("../../../WprimeWZ.root");
  
  TH1F* hnewRes1 = new TH1F("hnewRes1", "", 40, -200, 200);
  TH1F* hnewRes2 = new TH1F("hnewRes2", "", 40, -200, 200);

  TH1F* holdRes1 = new TH1F("holdRes1", "", 40, -200, 200);
  TH1F* holdRes2 = new TH1F("holdRes2", "", 40, -200, 200);

  TH1F* htryRes1 = new TH1F("htryRes1", "", 40, -200, 200);
  TH1F* htryRes2 = new TH1F("htryRes2", "", 40, -200, 200);

  TH1F* hfilRes1 = new TH1F("hfilRes1", "", 40, -200, 200);
  TH1F* hfilRes2 = new TH1F("hfilRes2", "", 40, -200, 200);

  TH1F* hhalRes1 = new TH1F("hhalRes1", "", 40, -200, 200);
  TH1F* hhalRes2 = new TH1F("hhalRes2", "", 40, -200, 200);

  TTree* tEvts = getTree(fIn, sample, "tEvts_MET"); assert(tEvts);
  int n = tEvts->Draw("ZLep1PtGen:ZLep1EtaGen:ZLep1PhiGen:ZLep2PtGen:ZLep2EtaGen:ZLep2PhiGen:ZLep1Pt:ZLep1PtErr:ZLep1Eta:ZLep1Phi:ZLep2Pt:ZLep2PtErr:ZLep2Eta:ZLep2Phi", Form("weight*(max(ZLep1PtGen,ZLep2PtGen)>800)*(EvtType>=2)"), "para goff");
  //int n = tEvts->Draw("ZLep1PtGen:ZLep1EtaGen:ZLep1PhiGen:ZLep2PtGen:ZLep2EtaGen:ZLep2PhiGen:ZLep1Pt:ZLep1Eta:ZLep1Phi:ZLep2Pt:ZLep2Eta:ZLep2Phi", Form("weight*(ZLep1PtGen>10)*(abs(max(ZLep1PtGen,ZLep2PtGen)-ZLep1Pt)>200)*(abs(min(ZLep1PtGen,ZLep2PtGen)-ZLep2Pt)>10)"), "para goff");
  for(int ievt=0; ievt<min(n,200); ++ievt){  //loop over events
    printf("--------------------\n");
    int treeIdx = 0;
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
    const float weight       = tEvts->GetW()[ievt];
    
    TLorentzVector l1_gen, l2_gen, sum_gen;
    if(ZLep1PtGen > ZLep2PtGen){
      l1_gen.SetPtEtaPhiM(   ZLep1PtGen,ZLep1EtaGen,ZLep1PhiGen,0.); 
      l2_gen.SetPtEtaPhiM(   ZLep2PtGen,ZLep2EtaGen,ZLep2PhiGen,0.); 
    }else{
      l2_gen.SetPtEtaPhiM(   ZLep1PtGen,ZLep1EtaGen,ZLep1PhiGen,0.); 
      l1_gen.SetPtEtaPhiM(   ZLep2PtGen,ZLep2EtaGen,ZLep2PhiGen,0.); 
    }
    sum_gen = l1_gen + l2_gen;
    float mZ_gen = sum_gen.M();

    TLorentzVector l1_old, l2_old, sum_old;
    l1_old.SetPtEtaPhiM(   ZLep1Pt,ZLep1Eta,ZLep1Phi,0.); 
    l2_old.SetPtEtaPhiM(   ZLep2Pt,ZLep2Eta,ZLep2Phi,0.); 
    sum_old = l1_old + l2_old;
    float mZ_old = sum_old.M();

    if(l1_gen.DeltaR(l1_old) > 0.1 || l2_gen.DeltaR(l2_old) > 0.1) continue;

    printf("  GEN: mz %.2f zpt %.2f zpz %.2f zp %.2f ze %.2f zeta %.2f zphi %.2f dr %.2f l1pt %.2f l2pt %.2f l1eta %.2f l2eta %.2f l1phi %.2f l2phi %.2f cosdphi %.4f\n", mZ_gen, sum_gen.Pt(),sum_gen.Pz(),sum_gen.P(),sum_gen.E(),sum_gen.Eta(), sum_gen.Phi(), l1_gen.DeltaR(l2_gen), l1_gen.Pt(),l2_gen.Pt(),l1_gen.Eta(),l2_gen.Eta(),l1_gen.Phi(),l2_gen.Phi(), (sum_gen.Perp2() - l1_gen.Perp2() - l2_gen.Perp2())/(2.*l1_gen.Perp()*l2_gen.Perp()));
    l1_gen.Print();
    l2_gen.Print();
    sum_gen.Print();
    printf("  OLD: mz %.2f zpt %.2f zpz %.2f zp %.2f ze %.2f zeta %.2f zphi %.2f dr %.2f l1pt %.2f l2pt %.2f l1eta %.2f l2eta %.2f l1phi %.2f l2phi %.2f l1pterr %.2f l2pterr %.2f cosdphi %.4f\n", mZ_old, sum_old.Pt(),sum_old.Pz(),sum_old.P(),sum_old.E(),sum_old.Eta(),sum_old.Phi(),l1_old.DeltaR(l2_old), l1_old.Pt(),l2_old.Pt(),l1_old.Eta(),l2_old.Eta(),l1_old.Phi(),l2_old.Phi(), ZLep1PtErr, ZLep2PtErr, (sum_old.Perp2() - l1_old.Perp2() - l2_old.Perp2())/(2.*l1_old.Perp()*l2_old.Perp()));
    l1_old.Print();
    l2_old.Print();
    sum_old.Print();
    printf("Want: k1: %.2f k2: %.2f arctan(k1/k2) = %.2f \n", l1_gen.Pt() / l1_old.Pt(), l2_gen.Pt() / l2_old.Pt(), TMath::ATan( (l1_gen.Pt()/l1_old.Pt()) / (l2_gen.Pt()/l2_old.Pt()) ) );

    printf("prod: gen Z+l1 %.4f gen z+l2 %.4f old z+l1 %.4f old z+l2 %.4f \n", 
           (sum_gen - l1_gen).Mag2(), (sum_gen - l2_gen).Mag2(), (sum_old - l1_old).Mag2(), (sum_old - l2_old).Mag2());

    float bestProd = 0.;
    float bestDiff = 999.;
    for(float k = 0.01; k<10.; k+=0.01){
      TLorentzVector temp1,temp2;
      temp1.SetPtEtaPhiM(k*ZLep1Pt,ZLep1Eta,ZLep1Phi,0.); 
      temp2.SetPtEtaPhiM(  ZLep2Pt,ZLep2Eta,ZLep2Phi,0.); 
      TLorentzVector temp = temp1 + temp2;
      float diff = abs(temp.M() - 91.188);
      if(diff < bestDiff){
        bestProd = k;
        bestDiff = diff;
      }
    }
    printf("Best k1*k2 = %.2f \n", bestProd);

    float bestResult = 9e9;
    float bestk1     = 1.;
    float bestk2     = 1.;
    for(float k1 = 0.1; k1<10.; k1+=0.01){
      for(float k2 = 0.1; k2<10.; k2+=0.01){
        //float k2 =  bestProd / k1;
        
        TLorentzVector l1_new, l2_new, sum_new;
        l1_new.SetPtEtaPhiM(k1*ZLep1Pt,ZLep1Eta,ZLep1Phi,0.); 
        l2_new.SetPtEtaPhiM(k2*ZLep2Pt,ZLep2Eta,ZLep2Phi,0.); 
        sum_new = l1_new + l2_new;

        if(l1_new.Pt() < l2_new.Pt()) continue;
        
        float Result = Analyze(k1,k2,l1_gen, l2_gen, sum_gen, l1_old, l2_old, sum_old, l1_new, l2_new, sum_new);

        if(Result < bestResult){
          //cout<<"Best is now "<<Result<<" beating out "<<bestResult<<endl;
          bestResult = Result;
          bestk1 = k1;
          bestk2 = k2;
        }
      }
    }
    cout<<"Best was "<<bestk1<<" and "<<bestk2<<endl;
    TLorentzVector l1_new, l2_new, sum_new;
    l1_new.SetPtEtaPhiM(bestk1*ZLep1Pt,ZLep1Eta,ZLep1Phi,0.); 
    l2_new.SetPtEtaPhiM(bestk2*ZLep2Pt,ZLep2Eta,ZLep2Phi,0.); 
    sum_new = l1_new + l2_new;
    Analyze(bestk1,bestk2,l1_gen, l2_gen, sum_gen, l1_old, l2_old, sum_old, l1_new, l2_new, sum_new, true);

    hnewRes1->Fill(l1_new.Pt() - l1_gen.Pt());
    hnewRes2->Fill(l2_new.Pt() - l2_gen.Pt());
    cout<<"New Res1 is now "<<hnewRes1->GetRMS()<<" +- "<<hnewRes1->GetRMSError()<<endl;
    cout<<"New Res2 is now "<<hnewRes2->GetRMS()<<" +- "<<hnewRes2->GetRMSError()<<endl;

    float wantk1 = l1_gen.Pt() / l1_old.Pt();
    float wantk2 = l2_gen.Pt() / l2_old.Pt();
    cout<<"Want was "<<wantk1<<" and "<<wantk2<<endl;
    l1_new.SetPtEtaPhiM(wantk1*ZLep1Pt,ZLep1Eta,ZLep1Phi,0.); 
    l2_new.SetPtEtaPhiM(wantk2*ZLep2Pt,ZLep2Eta,ZLep2Phi,0.); 
    sum_new = l1_new + l2_new;
    Analyze(wantk1,wantk2,l1_gen, l2_gen, sum_gen, l1_old, l2_old, sum_old, l1_new, l2_new, sum_new, true);

    //Try Me!
    float tryk1 = pow(bestProd, ZLep1Pt / (ZLep1Pt + ZLep2Pt)) ; //Weight by pt
    float tryk2 = pow(bestProd, ZLep2Pt / (ZLep1Pt + ZLep2Pt)) ; //Weight by pt
    //Assume softer muon is better measured
    //float tryk1 = bestProd;
    //float tryk2 = 1.;


    cout<<"Try was "<<tryk1<<" and "<<tryk2<<endl;
    l1_new.SetPtEtaPhiM(tryk1*ZLep1Pt,ZLep1Eta,ZLep1Phi,0.); 
    l2_new.SetPtEtaPhiM(tryk2*ZLep2Pt,ZLep2Eta,ZLep2Phi,0.); 
    sum_new = l1_new + l2_new;
    Analyze(tryk1,tryk2,l1_gen, l2_gen, sum_gen, l1_old, l2_old, sum_old, l1_new, l2_new, sum_new, true);
    htryRes1->Fill(l1_new.Pt() - l1_gen.Pt());
    htryRes2->Fill(l2_new.Pt() - l2_gen.Pt());
    cout<<"Try Res1 is now "<<htryRes1->GetRMS()<<" +- "<<htryRes1->GetRMSError()<<endl;
    cout<<"Try Res2 is now "<<htryRes2->GetRMS()<<" +- "<<htryRes2->GetRMSError()<<endl;

    //Fil Me!
    float filk1 = pow(bestProd, ZLep1Pt / (ZLep1Pt + ZLep2Pt)) ; //Weight by sigma
    float filk2 = pow(bestProd, ZLep2Pt / (ZLep1Pt + ZLep2Pt)) ; //Weight by sigma
    float A = (ZLep1Pt / ZLep1PtErr) / (ZLep2Pt/ZLep2PtErr);
    float bestFil = 9e9;
    for(float k1 = 0.1; k1<10.; k1+=0.01){
      float k2 = bestProd / k1;
      float B = fabs((1-k2) / (1-k1));
      float C = fabs( 1. - k1 ) * ZLep1Pt / ZLep1PtErr + fabs( 1. - k2 ) * ZLep2Pt / ZLep2PtErr;
      //if( C > 10) continue;
      if(C < bestFil){
        //if(fabs( A - B ) < bestFil){
        bestFil = C;//fabs( A - B );
        printf("fil: k1 %.2f k2 %.2f A %.2f B: %.2f C: %.2f diff %.4f\n", k1, k2, A, B, C, fabs( A - B ));
        filk1 = k1;
        filk2 = k2;
      }
    }

    cout<<"Fil was "<<filk1<<" and "<<filk2<<endl;
    l1_new.SetPtEtaPhiM(filk1*ZLep1Pt,ZLep1Eta,ZLep1Phi,0.); 
    l2_new.SetPtEtaPhiM(filk2*ZLep2Pt,ZLep2Eta,ZLep2Phi,0.); 
    sum_new = l1_new + l2_new;
    Analyze(filk1,filk2,l1_gen, l2_gen, sum_gen, l1_old, l2_old, sum_old, l1_new, l2_new, sum_new, true);
    hfilRes1->Fill(l1_new.Pt() - l1_gen.Pt());
    hfilRes2->Fill(l2_new.Pt() - l2_gen.Pt());
    cout<<"Fil Res1 is now "<<hfilRes1->GetRMS()<<" +- "<<hfilRes1->GetRMSError()<<endl;
    cout<<"Fil Res2 is now "<<hfilRes2->GetRMS()<<" +- "<<hfilRes2->GetRMSError()<<endl;


    //Hal Me!, weight
    float halk1 = 1.;
    float halk2 = 1.;
    float bestProb = -1;
    TF1 zprob ( "zprob","gaus",60, 120);  zprob.SetParameters(1,91.188, 1.);
    TF1 k1prob("k1prob","gaus",0,10); k1prob.SetParameters(1,1., ZLep1PtErr/ZLep1Pt); 
    TF1 k2prob("k2prob","gaus",0,10); k2prob.SetParameters(1,1., ZLep2PtErr/ZLep2Pt); 
    for(float k1 = 0.1; k1<5.; k1+=0.001){
      for(float k2 = 0.1; k2<5.; k2+=0.001){
        TLorentzVector temp1,temp2;
        temp1.SetPtEtaPhiM(k1*ZLep1Pt,ZLep1Eta,ZLep1Phi,0.); 
        temp2.SetPtEtaPhiM(k2*ZLep2Pt,ZLep2Eta,ZLep2Phi,0.); 
        TLorentzVector temp = temp1 + temp2;
        float prob = zprob(temp.M()) * k1prob(k1) * k2prob(k2);
        if(prob > bestProb){
          printf("prob is %.4f k1: %.2f k2: %.2f\n", bestProb, k1, k2);
          bestProb = prob;
          halk1 = k1;
          halk2 = k2;
        }
      }
    }

    cout<<"Hal was "<<halk1<<" and "<<halk2<<endl;
    l1_new.SetPtEtaPhiM(halk1*ZLep1Pt,ZLep1Eta,ZLep1Phi,0.); 
    l2_new.SetPtEtaPhiM(halk2*ZLep2Pt,ZLep2Eta,ZLep2Phi,0.); 
    sum_new = l1_new + l2_new;
    Analyze(halk1,halk2,l1_gen, l2_gen, sum_gen, l1_old, l2_old, sum_old, l1_new, l2_new, sum_new, true);
    hhalRes1->Fill(l1_new.Pt() - l1_gen.Pt());
    hhalRes2->Fill(l2_new.Pt() - l2_gen.Pt());
    cout<<"Hal Res1 is now "<<hhalRes1->GetRMS()<<" +- "<<hhalRes1->GetRMSError()<<endl;
    cout<<"Hal Res2 is now "<<hhalRes2->GetRMS()<<" +- "<<hhalRes2->GetRMSError()<<endl;

    //which boost gives half eneregy?

    //Old
    holdRes1->Fill(l1_old.Pt() - l1_gen.Pt());
    holdRes2->Fill(l2_old.Pt() - l2_gen.Pt());
    cout<<"Old Res1 is now "<<holdRes1->GetRMS()<<" +- "<<holdRes1->GetRMSError()<<" after iteration "<<holdRes1->GetEntries()<<endl;
    cout<<"Old Res2 is now "<<holdRes2->GetRMS()<<" +- "<<holdRes2->GetRMSError()<<endl;

  }//Evt loop

  //hnewRes1->Draw();
  //hnewRes2->Draw("same");

  htryRes1->Draw();
  htryRes2->SetLineColor(kRed);
  htryRes2->Draw("sames");
}

float
Analyze(float k1, float k2, 
        TLorentzVector & l1_gen, TLorentzVector & l2_gen, TLorentzVector &sum_gen,
        TLorentzVector & l1_old, TLorentzVector & l2_old, TLorentzVector &sum_old,
        TLorentzVector & l1_new, TLorentzVector & l2_new, TLorentzVector &sum_new,
        bool print){
      float mZ_new = sum_new.M();

      //if(fabs(sum_new.P() - pt(l1_new.DeltaR(l2_new))) > 10 ||
      //fabs(mZ_new - 91.188                        ) > 1.) return 9e9;

      
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
 
      return pow(sum_new.Pt() - fn_pt,2) + pow(50*(mZ_new - 91.188),2) + 100*(k1*k1 + k2*k2);
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
