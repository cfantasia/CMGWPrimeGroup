//Usage: root -b -q 'MatrixCalc.C+(useEWK)'
#include <iostream>
#include <iomanip>

#include "../Limits/consts.h"

using namespace std;

const float eElec = 0.9359; const float SeElec = 0.0004;//done
const float pElec = 0.07/**(1. - 0.3381)*/; const float SpElec = 0.01;//done

const float eMuon = 0.9711; const float SeMuon = 0.0002;
const float pMuon = 0.06/**(1. - -1*0.3021)*/; const float SpMuon = 0.01;//done

struct MatrixResult{
  float Nt_lep;
  float DeltaNt_lep;
  float Nt_jet;
  float DeltaNt_jet;
};

void MatrixMethod(const TH1F *hLoose, const TH1F *hTight, bool allChannel=true);
MatrixResult CalcMatrix(const float eTight, const  float Delta_eTight,
                        const float pFake, const float Delta_pFake, 
                        const float Nl, const float Delta_Nl,
                        const float Nt, const float Delta_Nt);


void MatrixCalc(bool useEWK=true){
  cout<<"useEWK: "<<useEWK<<endl;

  string tightfilename;
  string loosefilename;
  
  if(useEWK){
    tightfilename = "EWKWZ.root";
    loosefilename = "EWKWZMatrix.root";
  }else{
    tightfilename = "../../../WprimeWZ.root";
    loosefilename = "../../../WprimeWZMatrix.root";
  }

  TFile *fTight = TFile::Open(tightfilename.c_str(), "read"); assert(fTight);
  TFile *fLoose = TFile::Open(loosefilename.c_str(), "read"); assert(fLoose);

  vector<string> vBkg, vMC, samples;

  vBkg.push_back("DYJetsToLL");
  vBkg.push_back("TTJets");
  vBkg.push_back("ZZ");
  vBkg.push_back("GVJets");  
  vBkg.push_back("WJetsToLNu");
  vBkg.push_back("WWTo2L2Nu");

  if(!useEWK)  vBkg.push_back("WZJetsTo3LNu");

  vMC = vBkg;

  if(useEWK) vMC.push_back("WZJetsTo3LNu");

  samples = vMC;
  samples.push_back("BKG");
  samples.push_back("MC");
  if(!useEWK) samples.push_back("EXO");
  samples.push_back("data");

  cout << setiosflags(ios::fixed) << setprecision(2);    
  if(useEWK){
    for(unsigned i=0; i<samples.size(); ++i){
      vector<string> names;
      if     (samples[i] == "MC"){
        names = vMC;
      }else if(samples[i] == "BKG"){
        names = vBkg;
      }else{
        names.push_back(samples[i]);
      }
      TH1F* hTight = get_sum_of_hists(fTight, names, "hEvtType_AllCuts", 0);
      TH1F* hLoose = get_sum_of_hists(fLoose, names, "hEvtType_AllCuts", 0);
      MatrixMethod(hLoose, hTight, samples[i]=="data");
    }
  }else{
    TTree* tEvts = new TTree("tEvts", "Cut Values per sample");
    tEvts->ReadFile("../Limits/cutValues.wz.dat");
    tEvts->Draw("SignalCode:Mass:minWindow:maxWindow:HtCut:ZptCut:WptCut","", "para goff");
    float ntEvts = tEvts->GetSelectedRows(); cout<<"Found "<<ntEvts<<endl;

    bool useHt(1), useZpt(0), useWpt(0), useWindow(1);

    for(int isignal=0; isignal<ntEvts; ++isignal){
      const int SignalCode = tEvts->GetVal(0)[isignal];
      const string SignalName = SampleName(SignalCode)[0];
      const float mass = tEvts->GetVal(1)[isignal];
      const float minWindow = useWindow ? tEvts->GetVal(2)[isignal] : 0;
      const float maxWindow = useWindow ? tEvts->GetVal(3)[isignal] : 99999;
      const float minHt     = useHt     ? tEvts->GetVal(4)[isignal] : 0;
      const float minZpt    = useZpt    ? tEvts->GetVal(5)[isignal] : 0;
      const float minWpt    = useWpt    ? tEvts->GetVal(6)[isignal] : 0;
      
      
      const string cuts = Form("(WZMass > %.0f && WZMass < %.0f && Ht > %.0f && Zpt > %.0f && Wpt > %.0f)*weight",
                               minWindow, maxWindow, minHt, minZpt, minWpt); 
      cout<<"\n-------\nSignal: "<<SignalName<<" with cuts "<<cuts<<endl;
      
      for(unsigned i=0; i<samples.size(); ++i){
        vector<string> names;
        if     (samples[i] == "MC"){
          names = vMC;
          names.push_back(SignalName);
        }else if(samples[i] == "BKG"){
          names = vBkg;
        }else if(samples[i] == "EXO"){
          names.push_back(SignalName);
        }else{
          //if(samples[i] != "data") continue;//Cory: Remove
          names.push_back(samples[i]);
        }
//        cout<<"---------------------------------\n";
        cout<<"sample now is "<<samples[i]<<endl;
        
        TH1F hTight("hTight", "hTight", 4, 0, 4);
        get_sum_of_hists(fTight, names, "tEvts_ValidWZCand", "EvtType", cuts, hTight);

        TH1F hLoose("hLoose", "hLoose", 4, 0, 4);
        get_sum_of_hists(fLoose, names, "tEvts_ValidWZCand", "EvtType", cuts, hLoose);
        
        MatrixMethod(&hLoose, &hTight, samples[i]=="data");
      }//samples loop
    }//signal loop
  }
}

void
MatrixMethod(const TH1F *hLoose, const TH1F *hTight, bool allChannels){
  float e3Tight = hTight->GetBinContent(1);
  float e2Tight = hTight->GetBinContent(2);
  float e1Tight = hTight->GetBinContent(3);
  float e0Tight = hTight->GetBinContent(4);  
  float allTight = e3Tight + e2Tight + e1Tight + e0Tight;
  
  float Se3Tight = sqrt(e3Tight);
  float Se2Tight = sqrt(e2Tight);
  float Se1Tight = sqrt(e1Tight);
  float Se0Tight = sqrt(e0Tight);
  float SallTight = sqrt(allTight);
  
  cout<<"Tight Sample: "<<" N total: "<<Value(allTight,SallTight)<<" +/- "<<Value(SallTight)<<endl;

  if(0) return;

  if(allChannels)
    cout<<"N 3e: "<<e3Tight<<" +/- "<<Se3Tight<<endl//" frac: "<<e3/total<<endl
        <<"N 2e: "<<e2Tight<<" +/- "<<Se2Tight<<endl//" frac: "<<e2/total<<endl
        <<"N 1e: "<<e1Tight<<" +/- "<<Se1Tight<<endl//" frac: "<<e1/total<<endl
        <<"N 0e: "<<e0Tight<<" +/- "<<Se0Tight<<endl//" frac: "<<e0/total<<endl
        <<endl;

  float e3Loose = hLoose->GetBinContent(1);
  float e2Loose = hLoose->GetBinContent(2);
  float e1Loose = hLoose->GetBinContent(3);
  float e0Loose = hLoose->GetBinContent(4);    
  float allLoose = e3Loose + e2Loose + e1Loose + e0Loose;
  
  float Se3Loose = sqrt(e3Loose);
  float Se2Loose = sqrt(e2Loose);
  float Se1Loose = sqrt(e1Loose);
  float Se0Loose = sqrt(e0Loose);
  float SallLoose = sqrt(allLoose);
  
  cout<<"Loose Sample: "<<" N total: "<<Value(allLoose,SallLoose)<<" +/- "<<Value(SallLoose)<<endl;
  if(allChannels)
    cout<<"N 3e: "<<e3Loose<<" +/- "<<Se3Loose<<endl//" frac: "<<e3/total<<endl
        <<"N 2e: "<<e2Loose<<" +/- "<<Se2Loose<<endl//" frac: "<<e2/total<<endl
        <<"N 1e: "<<e1Loose<<" +/- "<<Se1Loose<<endl//" frac: "<<e1/total<<endl
        <<"N 0e: "<<e0Loose<<" +/- "<<Se0Loose<<endl//" frac: "<<e0/total<<endl
        <<endl;
  
/////////////////////////////
  cout<<" For All Channels\n";
  CalcMatrix(eElec, SeElec,//Cory: This is wrong.  Have to weight elec & muon eff/fake rates
             pElec, SpElec,
             allLoose, SallLoose,
             allTight, SallTight);

  if(!allChannels) return;

  cout<<" For 3e\n";
  MatrixResult result3e = CalcMatrix(eElec, SeElec,
                                     pElec, SpElec, 
                                     e3Loose, Se3Loose,
                                     e3Tight, Se3Tight);
  cout<<" For 2e\n";
  MatrixResult result2e = CalcMatrix(eMuon, SeMuon,
                                     pMuon, SpMuon, 
                                     e2Loose, Se2Loose,
                                     e2Tight, Se2Tight);
  cout<<" For 1e\n";
  MatrixResult result1e = CalcMatrix(eElec, SeElec,
                                     pElec, SpElec, 
                                     e1Loose, Se1Loose,
                                     e1Tight, Se1Tight);
  cout<<" For 0e\n";
  MatrixResult result0e = CalcMatrix(eMuon, SeMuon,
                                     pMuon, SpMuon, 
                                     e0Loose, Se0Loose,
                                     e0Tight, Se0Tight);
  
  MatrixResult resultAll;
  resultAll.Nt_lep = result3e.Nt_lep
    + result2e.Nt_lep
    + result1e.Nt_lep 
    + result0e.Nt_lep;

  resultAll.Nt_jet = result3e.Nt_jet +
    result2e.Nt_jet +
    result1e.Nt_jet +
    result0e.Nt_jet;

  resultAll.DeltaNt_lep = sqrt(result3e.DeltaNt_lep*result3e.DeltaNt_lep +
                               result2e.DeltaNt_lep*result2e.DeltaNt_lep +
                               result1e.DeltaNt_lep*result1e.DeltaNt_lep +
                               result0e.DeltaNt_lep*result0e.DeltaNt_lep);

  resultAll.DeltaNt_jet = sqrt(result3e.DeltaNt_jet*result3e.DeltaNt_jet +
                               result2e.DeltaNt_jet*result2e.DeltaNt_jet +
                               result1e.DeltaNt_jet*result1e.DeltaNt_jet +
                               result0e.DeltaNt_jet*result0e.DeltaNt_jet);

  cout<<" & Tight nLep +- error & Tight nJet +- error "<<endl
      <<" & "<<Value(resultAll.Nt_lep,resultAll.DeltaNt_lep)
      <<" & "<<Value(resultAll.DeltaNt_lep)
      <<" & "<<Value(resultAll.Nt_jet,resultAll.DeltaNt_jet)
      <<" & "<<Value(resultAll.DeltaNt_jet)
      <<" &  & \\\\ \\hline"<<endl;
}

MatrixResult
CalcMatrix(const float eTight, const  float Delta_eTight, 
           const float pFake, const float Delta_pFake, 
           const float Nl, const float Delta_Nl,
           const float Nt, const float Delta_Nt){
    MatrixResult result;

    float N1 = Nl - Nt;
    float N2 = Nt;

    float dN1 = TMath::Sqrt(N1);//Vuko's improvement
    //float dN1 = TMath::Sqrt(Delta_Nl*Delta_Nl + Delta_Nt*Delta_Nt);  
    float dN2 = Delta_Nt;
    

    float Nlep = (Nt - pFake*Nl)/(eTight - pFake);
    float Njet = (eTight*Nl - Nt)/(eTight - pFake);
    //cout << "Loose: N_lep: " << Nlep << " and: N_jet: " << Njet << endl;
  
    result.Nt_lep = eTight*Nlep;
    result.Nt_jet = pFake*Njet;
    
    //cout << "Tight: N_lep: " << Nt_lep << " and: N_jet: " << Nt_jet << endl;  
    
    //errors
    
    float dNlep_desig = (pFake*(pFake*(N1+N2)-N2))/((eTight-pFake)*(eTight-pFake));
    
    float dNlep_deqcd = (eTight*(N2-eTight*(N1+N2)))/((eTight-pFake)*(eTight-pFake));   
    
    float dNlep_dN1 = (-eTight*pFake)/(eTight-pFake);
    
    float dNlep_dN2 = (eTight*(1-pFake))/(eTight-pFake);


    float dNjet_desig =  (pFake*(N2-pFake*(N1+N2)))/((eTight-pFake)*(eTight-pFake)); 

    float dNjet_deqcd = (eTight*(eTight*(N1+N2)-N2))/((eTight-pFake)*(eTight-pFake));

    float dNjet_dN1 = (pFake*eTight)/(eTight-pFake);

    float dNjet_dN2 = (pFake*(eTight-1))/(eTight-pFake);


    result.DeltaNt_lep = TMath::Sqrt(((dNlep_desig)*(dNlep_desig)*Delta_eTight*Delta_eTight)+((dNlep_deqcd*dNlep_deqcd)*Delta_pFake*Delta_pFake) + ((dNlep_dN1*dNlep_dN1)*dN1*dN1) + ((dNlep_dN2*dNlep_dN2)*dN2*dN2)); 

    result.DeltaNt_jet = TMath::Sqrt(((dNjet_desig)*(dNjet_desig)*Delta_eTight*Delta_eTight)+((dNjet_deqcd*dNjet_deqcd)*Delta_pFake*Delta_pFake) + ((dNjet_dN1*dNjet_dN1)*dN1*dN1) + ((dNjet_dN2*dNjet_dN2)*dN2*dN2));
  
    //cout << "Error on tNlep: " << DeltaNt_lep << " Total: " << Nt_lep << " +- " << DeltaNt_lep << endl;

    //cout << "Error on tNjet: " << DeltaNt_jet << " Total: " << Nt_jet << " +- " << DeltaNt_jet << endl;
/*
    cout<<"Debug Results are: \n"
        <<"e: "<<eTight*100<<" +/- "<<Delta_eTight*100<<"%" 
        <<" p: "<<pFake*100<<" +/- "<<Delta_pFake*100<<"%"
        <<endl
        <<"nLoose: "<<Nl<<" +/- "<<Delta_Nl
        <<" nTight: "<<Nt<<" +/- "<<Delta_Nt
        <<" ratio: "<<Nt/Nl*100<<"%"
        <<endl
        <<"Loose nLep: "<<Nlep
        <<" Loose nJet: "<<Njet
        <<endl
        <<"Tight nLep: "<<result.Nt_lep<<" +/- "<<result.DeltaNt_lep
        <<" Tight nJet: "<<result.Nt_jet<<" +/- "<<result.DeltaNt_jet
        <<endl;
*/
    return result;
}
