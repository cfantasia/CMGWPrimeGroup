//Usage: root -b -l -q 'MatrixCalc.C+(useEWK)'
#include <iostream>
#include <iomanip>

#include "../Limits/consts.h"
#include "THStack.h"

using namespace std;

const double eElec = 0.9359; const double SeElec = 0.0004;//done
const double pElec = 0.07/**(1. - 0.3381)*/; const double SpElec = 0.01;//done

const double eMuon = 0.9711; const double SeMuon = 0.0002;
const double pMuon = 0.06/**(1. - -1*0.3021)*/; const double SpMuon = 0.01;//done

struct Rates{
  double Reff;
  double Delta_Reff;
  double Rfake;
  double Delta_Rfake;
  Rates():Reff(0),Delta_Reff(0),Rfake(0),Delta_Rfake(0){}
  Rates(double e, double se, double f, double sf):Reff(e),Delta_Reff(se),Rfake(f),Delta_Rfake(sf){}
  Rates(Value e, Value p):Reff(e.val),Delta_Reff(e.err),Rfake(p.val),Delta_Rfake(p.err){}
};

struct MatrixResult{
  double Nt_lep;
  double DeltaNt_lep;
  double Nt_jet;
  double DeltaNt_jet;
  
  MatrixResult():Nt_lep(0),DeltaNt_lep(0),Nt_jet(0),DeltaNt_jet(0){}
  MatrixResult(double nl, double snl, double nj, double snj):Nt_lep(nl),DeltaNt_lep(snl),Nt_jet(nj),DeltaNt_jet(snj){}
  
  MatrixResult operator+( const MatrixResult &rhs ) const{
    MatrixResult result = *this;     // Make a copy of myself.  Same as MatrixResult result(*this);
    result += rhs;            // Use += to add other to the copy.
    return result;
  }
  MatrixResult & operator+=(const MatrixResult &rhs) {
    Nt_lep += rhs.Nt_lep;
    DeltaNt_lep = sqrt(pow(DeltaNt_lep,2) + pow(rhs.DeltaNt_lep,2));
    Nt_jet += rhs.Nt_jet;
    DeltaNt_jet = sqrt(pow(DeltaNt_jet,2) + pow(rhs.DeltaNt_jet,2));
    return *this;
  }

  friend std::ostream& operator << (std::ostream &o, const MatrixResult & m){
    o<<"NT_Lep: "<<m.Nt_lep<<" +/- "<<m.DeltaNt_lep<<endl;
    o<<"NT_Jet: "<<m.Nt_jet<<" +/- "<<m.DeltaNt_jet<<endl;
    //o<<"NT_Lep: "<<Value(m.Nt_lep,m.DeltaNt_lep)<<" +/- "<<Value(m.DeltaNt_lep)<<endl;
    //o<<"NT_Jet: "<<Value(m.Nt_jet,m.DeltaNt_jet)<<" +/- "<<Value(m.DeltaNt_jet)<<endl;
    return o;
  }
  
};

struct MatrixSomething{
  double NGood;
  double DeltaNGood;
  double NBad;
  double DeltaNBad;
  
  MatrixSomething():NGood(0),DeltaNGood(0),NBad(0),DeltaNBad(0){}
  MatrixSomething(double nl, double snl, double nj, double snj):NGood(nl),DeltaNGood(snl),NBad(nj),DeltaNBad(snj){}
  
  MatrixSomething operator+( const MatrixSomething &rhs ) const{
    MatrixSomething result = *this;     // Make a copy of myself.  Same as MatrixSomething result(*this);
    result += rhs;            // Use += to add other to the copy.
    return result;
  }
  MatrixSomething & operator+=(const MatrixSomething &rhs) {
    NGood += rhs.NGood;
    DeltaNGood = sqrt(pow(DeltaNGood,2) + pow(rhs.DeltaNGood,2));
    NBad += rhs.NBad;
    DeltaNBad = sqrt(pow(DeltaNBad,2) + pow(rhs.DeltaNBad,2));
    return *this;
  }

  friend std::ostream& operator << (std::ostream &o, const MatrixSomething & m){
    o<<"NT_Good: "<<m.NGood<<" +/- "<<m.DeltaNGood<<endl;
    o<<"NT_Bad: "<<m.NBad<<" +/- "<<m.DeltaNBad<<endl;
    //o<<"NT_Lep: "<<Value(m.NGood,m.DeltaNGood)<<" +/- "<<Value(m.DeltaNGood)<<endl;
    //o<<"NT_Jet: "<<Value(m.NBad,m.DeltaNBad)<<" +/- "<<Value(m.DeltaNBad)<<endl;
    return o;
  }
  
};

double FindNEvts(TFile* f, vector<string>& samples, string hName, double pt1, double eta1, double pt2, double eta2, float &err);
double GetNEvts(TH2F* h, double xmin, double ymin, double xmax, double ymax, float& err);

Value FindE(double pt, double eta, bool isElec, bool isW);
Value FindP(double pt, double eta, bool isElec, bool isW);

void MatrixMethod(const TH1F *hLoose, const TH1F *hTight, bool allChannel=true);
MatrixResult CalcMatrix(const double eTight, const  double Delta_eTight,
                        const double pFake, const double Delta_pFake, 
                        const double Nl, const double Delta_Nl,
                        const double Nt, const double Delta_Nt);

MatrixSomething
CalcSomething(const Rates & wRates, const Rates & z1Rates, const Rates & z2Rates,
              const Value & vLoose, const Value & vTight, const Value & vTT, const Value & vTTT);


void MatrixCalc(bool useEWK=true){
  cout<<"useEWK: "<<useEWK<<endl;

  string tightfilename;
  string loosefilename;
  string ttfilename, tttfilename;
  
  if(useEWK){
    tightfilename = "../../../EWKWZ.root";
    loosefilename = "../../../EWKWZMatrix.root";
    ttfilename = "../../../EWKWZ-TT.root";
    tttfilename = "../../../EWKWZ-TT.root";
  }else{
    tightfilename = "../../../WprimeWZ.root";
    loosefilename = "../../../WprimeWZMatrix.root";
  }

  TFile *fTight = TFile::Open(tightfilename.c_str(), "read"); assert(fTight);
  TFile *fLoose = TFile::Open(loosefilename.c_str(), "read"); assert(fLoose);
  //TFile *fTT    = TFile::Open(   ttfilename.c_str(), "read"); assert(fTT);
  //TFile *fTTT    = TFile::Open( tttfilename.c_str(), "read"); assert(fTTT);

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
        if(samples[i] != "data") continue;//Cory: Remove
        names.push_back(samples[i]);
      }
      cout<<" Doing MM for sample "<<samples[i]<<endl;

      //TH1F* hTight = get_sum_of_hists(fTight, names, "hEvtType_AllCuts", 0);
      //TH1F* hLoose = get_sum_of_hists(fLoose, names, "hEvtType_AllCuts", 0);
      //MatrixMethod(hLoose, hTight, samples[i]=="data");

      MatrixResult resultCombined;
      MatrixResult resultAll[4];
      MatrixSomething somethingCombined;
      MatrixSomething somethingAll[4];
      Value vLooseCombined, vTightCombined, vTTCombined, vTTTCombined;
      Value vLooseAll[4], vTightAll[4], vTTAll[4], vTTTAll[4];

      if(samples[i] != "data") continue;//Cory: Remove this to do MC

      TTree* tLoose = (TTree*) fLoose->Get("data/tEvts_MET");
      TTree* tTight = (TTree*) fTight->Get("data/tEvts_MET");
      TTree* tTT  = (TTree*) fTight->Get("data/tEvts_ZLep1Tight");
      TTree* tTTT = (TTree*) fTight->Get("data/tEvts_ZLep2Tight");

      const int nchannel = 4;
      for(int channel=0; channel<nchannel; ++channel){
        double pt [] = {0., 20., 25., 30., 50., 1000}; const int npt  = 6;
        

        string histName = Form("hEtaVsPt%ie%im_AllCuts", nchannel-1-channel, channel);
        bool wElec = channel%2 == false;
        bool zElec = channel<=1;
        double Eleceta[] = {0., 1.5, 2.5}; int nEleceta = 3;
        double Muoneta[] = {0., 2.4}; int nMuoneta = 2;
        double *Weta, *Zeta;
        int nWeta, nZeta;
        if(wElec){ 
          Weta = Eleceta;
          nWeta = nEleceta;
        }else{
          Weta = Muoneta;
          nWeta = nMuoneta;
        }
        if(zElec){ 
          Zeta = Eleceta;
          nZeta = nEleceta;
        }else{
          Zeta = Muoneta;
          nZeta = nMuoneta;
        }
        printf("channel=%i, wElec=%i, zElec=%i, npt=%i, nWeta=%i\n", channel, wElec, zElec, npt, nWeta);
        
        //W Loops
        for(int iWeta=0; iWeta<nWeta-1; ++iWeta){
          for(int iWpt=0; iWpt<npt-1; ++iWpt){
            string wCuts = Form("(%.0f <= WLepPt && WLepPt < %.0f)*(%.1f <= abs(WLepEta) && abs(WLepEta) < %.1f)",
                                pt[iWpt], pt[iWpt+1], Weta[iWeta],Weta[iWeta+1]); 
            //Z1 Loops
            for(int iZ1eta=0; iZ1eta<nZeta-1; ++iZ1eta){
              for(int iZ1pt=0; iZ1pt<npt-1; ++iZ1pt){
                string z1Cuts = Form("(%.0f <= ZLep1Pt && ZLep1Pt < %.0f)*(%.1f <= abs(ZLep1Eta) && abs(ZLep1Eta) < %.1f)",
                                     pt[iZ1pt], pt[iZ1pt+1], Zeta[iZ1eta],Zeta[iZ1eta+1]); 
                //Z2 Loops
                for(int iZ2eta=0; iZ2eta<nZeta-1; ++iZ2eta){
                  for(int iZ2pt=0; iZ2pt<npt-1; ++iZ2pt){
                    string z2Cuts = Form("(%.0f <= ZLep2Pt && ZLep2Pt < %.0f)*(%.1f <= abs(ZLep2Eta) && abs(ZLep2Eta) < %.1f)",
                                         pt[iZ2pt], pt[iZ2pt+1], Zeta[iZ2eta],Zeta[iZ2eta+1]); 
                    string cuts = Form("weight*(EvtType==%i)*(%s)*(%s)*(%s)",channel, wCuts.c_str(), z1Cuts.c_str(), z2Cuts.c_str());
                    cout<<"Cuts are "<<cuts<<endl;

                    double weta = (Weta[iWeta]+Weta[iWeta+1])/2;
                    double wpt = (pt[iWpt]+pt[iWpt+1])/2;
                    double z1eta = (Zeta[iZ1eta]+Zeta[iZ1eta+1])/2;
                    double z1pt = (pt[iZ1pt]+pt[iZ1pt+1])/2;
                    double z2eta = (Zeta[iZ2eta]+Zeta[iZ2eta+1])/2;
                    double z2pt = (pt[iZ2pt]+pt[iZ2pt+1])/2;

                    Rates  wRates(FindE(wpt, weta, wElec, true), FindP(wpt, weta, wElec, true));
                    Rates z1Rates(FindE(z1pt, z1eta, zElec, false), FindP(z1pt, z1eta, zElec, false));
                    Rates z2Rates(FindE(z2pt, z2eta, zElec, false), FindP(z2pt, z2eta, zElec, false));

                    //Calculate Number of 3 real lepton events
                    Value vLoose, vTight, vTT, vTTT;
                    //vLoose.val = FindNEvts(fLoose, names, histName, pt[ipt], eta[ieta], pt[ipt+1], eta[ieta+1], vLoose.err);
                    //vTight.val = FindNEvts(fTight, names, histName, pt[ipt], eta[ieta], pt[ipt+1], eta[ieta+1], vTight.err);
                    //vTT.val  = FindNEvts(fTT, names, histName, pt[ipt], eta[ieta], pt[ipt+1], eta[ieta+1], vTT.err);
                    //vTTT.val = FindNEvts(fTTT, names, histName, pt[ipt], eta[ieta], pt[ipt+1], eta[ieta+1], vTTT.err);

                    tLoose->Draw("WLepPt:WLepEta:ZLep1Pt:ZLep1Eta:ZLep2Pt:ZLep2Eta", cuts.c_str(), "para goff");
                    vLoose.val = tLoose->GetSelectedRows();
                    vLoose.err = sqrt(vLoose.val);
                    vLooseAll[channel] += vLoose;
                    vLooseCombined += vLoose;
                    cout<<"Loose: "<<vLoose<<endl;
                    tTight->Draw("WLepPt:WLepEta:ZLep1Pt:ZLep1Eta:ZLep2Pt:ZLep2Eta", cuts.c_str(), "para goff");
                    vTight.val = tTight->GetSelectedRows();
                    vTight.err = sqrt(vTight.val);
                    vTightAll[channel] += vTight;
                    vTightCombined += vTight;
                    cout<<"Tight: "<<vTight<<endl;
                    tTT->Draw("WLepPt:WLepEta:ZLep1Pt:ZLep1Eta:ZLep2Pt:ZLep2Eta", cuts.c_str(), "para goff");
                    vTT.val = tTT->GetSelectedRows();
                    vTT.err = sqrt(vTT.val);
                    vTTAll[channel] += vTT;
                    vTTCombined += vTT;
                    cout<<"TT: "<<vTT<<endl;
                    tTTT->Draw("WLepPt:WLepEta:ZLep1Pt:ZLep1Eta:ZLep2Pt:ZLep2Eta", cuts.c_str(), "para goff");
                    vTTT.val = tTTT->GetSelectedRows();
                    vTTT.err = sqrt(vTTT.val);
                    vTTTAll[channel] += vTTT;
                    vTTTCombined += vTTT;
                    cout<<"TTT: "<<vTTT<<endl;

                    

                    //These asserts fail when W in Loose sample doesn't pass tight selection but another lepton would
                    //skip for now
                    //if(vLoose.val < vTight.val || vTight.val < vTT.val || vTTT.val < vTT.val) continue;
                    if(vLoose.val < vTight.val) vLoose.val = vTight.val;
                    assert(vLoose.val >= vTight.val);
                    assert(vTight.val >= vTT.val);
                    assert(vTT.val >= vTTT.val);
                    
                    MatrixResult tmpResult = CalcMatrix(wRates.Reff, wRates.Delta_Reff,
                                                        wRates.Rfake, wRates.Delta_Rfake,
                                                        vLoose.val, vLoose.err,
                                                        vTight.val, vTight.err);
                    
                    cout<<" The MM results from channel="<<channel<<" pt="<<wpt<<" eta="<<weta<<" are\n"<<tmpResult
                        <<"and NLoose="<<vLoose.val<<" NTight="<<vTight.val<<" ( ratio ="<<vTight.val/vLoose.val<<")"<<endl;
                    printf("Calculating with W : e=%.2f +- %.2f, p=%.2f +- %.2f\n", wRates.Reff, wRates.Delta_Reff, wRates.Rfake, wRates.Delta_Rfake);
                    printf("Calculating with Z1: e=%.2f +- %.2f, p=%.2f +- %.2f\n", z1Rates.Reff, z1Rates.Delta_Reff, z1Rates.Rfake, z1Rates.Delta_Rfake);
                    printf("Calculating with Z2: e=%.2f +- %.2f, p=%.2f +- %.2f\n", z2Rates.Reff, z2Rates.Delta_Reff, z2Rates.Rfake, z2Rates.Delta_Rfake);
                    
                    resultAll[channel] += tmpResult;
                    
                    MatrixSomething tmpSomething = CalcSomething( wRates, z1Rates, z2Rates,
                                                                  vLoose, vTight, vTT, vTTT);
                    
                    cout<<" The MMM results from channel="<<channel<<" pt="<<wpt<<" eta="<<weta<<" are\n"<<tmpSomething
                        <<"and NTight="<<vTight.val<<" NTT="<<vTT.val<<" ( ratio ="<<vTT.val/vTight.val<<")"
                        <<"and NTTT= "<<vTTT.val<<endl;
                    
                    somethingAll[channel] += tmpSomething;
                  }//z2pt
                }//z2eta
              }//z1pt
            }//Z1eta
          }//Wpt
        }//Weta

        cout<<"Loose: "<<vLooseAll[channel]<<" Tight: "<<vTightAll[channel]<<" TT: "<<vTTAll[channel]<<" TTT: "<<vTTTAll[channel]<<endl;
        cout<<" MM Channel "<<channel<<" Total is \n"<<resultAll[channel]<<endl;
        cout<<" MMM Channel "<<channel<<" Total is \n"<<somethingAll[channel]<<endl;
        resultCombined    += resultAll[channel];
        somethingCombined += somethingAll[channel];
        
      }//channel
      cout<<"Loose: "<<vLooseCombined<<" Tight: "<<vTightCombined<<" TT: "<<vTTCombined<<" TTT: "<<vTTTCombined<<endl;
      cout<<" MM All Channel Total is \n"<<resultCombined<<endl;
      cout<<" MMM All Channel Total is \n"<<somethingCombined<<endl;

      //Do Matrix Method to predict MET Dist
      //Cory: Needs to take into account the new results

      cout<<"Now do it bin by bin in MET\n";
      TCanvas c1;
      c1.Divide(2,2);
     
      TFile* fMM = new TFile("DataDrivenMET.root", "recreate");
      fMM->mkdir("DYJetsToLL-DataDriven");
        
      MatrixResult resultMETCombined;
      MatrixResult resultMETAll[4];
      TH1F* hResultMET[4];
      THStack* hMCMET[4];
      for(int channel=0; channel<nchannel; ++channel){
        string histName = Form("hMET%ie%im_ValidW", nchannel-1-channel, channel);
        bool wElec = channel%2 == false;
        printf("channel=%i, wElec=%i\n", channel, wElec);
        
        double pFake(0), SpFake(0), eTight(0), SeTight(0);
        double pFakeElec(0.0631), SpFakeElec(0.0086), eTightElec(0.9217), SeTightElec(0.0003);
        double pFakeMuon(0.0423), SpFakeMuon(0.0075), eTightMuon(0.9644), SeTightMuon(0.0002);

        if(wElec){
          pFake =  pFakeElec;
          SpFake = SpFakeElec;
          eTight = eTightElec;
          SeTight = SeTightElec;
        }else{
          pFake =  pFakeMuon;
          SpFake = SpFakeMuon;
          eTight = eTightMuon;
          SeTight = SeTightMuon;
        }
        
        TH1F* hLoose = get_sum_of_hists(fLoose, names, histName, 0);
        TH1F* hTight = get_sum_of_hists(fTight, names, histName, 0);
        
        fMM->cd("DYJetsToLL-DataDriven");
        hResultMET[channel] = (TH1F*) hLoose->Clone(Form("hMET%ie%im_ValidW", nchannel-1-channel, channel));

        hMCMET[channel] = new THStack(Form("hMCMET%ie%im_ValidW", nchannel-1-channel, channel), "");
        TH1F* hzjets = (TH1F*) fTight->Get(Form("DYJetsToLL/%s", histName.c_str())); assert(hzjets); 
        TH1F* httbar = (TH1F*) fTight->Get(Form("TTJets/%s", histName.c_str()));  assert(httbar);
        hzjets->SetFillColor(kOrange+7);
        httbar->SetFillColor(kViolet+2);

        hMCMET[channel]->Add(hzjets);
        hMCMET[channel]->Add(httbar);
        //; assert(hMCMET[channel]);

        int nbins = hLoose->GetNbinsX();
        for(int bin=1; bin<=nbins; ++bin){//loop over met bins
          int minbin = bin;
          int maxbin = bin;
          double min = hLoose->GetBinLowEdge(bin);
          double max = hLoose->GetBinLowEdge(bin+1);
        
          double NLoose(0), SNLoose(0);
          NLoose = hLoose->IntegralAndError(minbin, maxbin, SNLoose);
      
          double NTight(0), SNTight(0);      
          NTight = hTight->IntegralAndError(minbin, maxbin, SNTight);

        
          MatrixResult tmpResult = CalcMatrix(eTight, SeTight,
                                              pFake,  SpFake, 
                                              NLoose, SNLoose,
                                              NTight, SNTight);
        
          cout<<" The results from channel="<<channel<<" "<<min<<" < met <"<<max<<" are\n"<<tmpResult
              <<"and NLoose="<<NLoose<<" NTight="<<NTight<<endl;
          printf("Calculating with e=%.2f +- %.2f, p=%.2f +- %.2f\n", eTight, SeTight,pFake,  SpFake);

          resultMETAll[channel] = resultMETAll[channel] + tmpResult;
          hResultMET[channel]->SetBinContent(bin, tmpResult.Nt_jet);
          hResultMET[channel]->SetBinError  (bin, tmpResult.DeltaNt_jet);

        }//met bins
        cout<<" Channel "<<channel<<" Total is \n"<<resultMETAll[channel]<<endl;
        resultMETCombined = resultMETCombined + resultMETAll[channel];
        c1.cd(channel+1);
        hMCMET[channel]->Draw("hist");
        hResultMET[channel]->Draw("e1 same");
        hMCMET[channel]->SetMaximum(1.05* max(hMCMET[channel]->GetMaximum(), hResultMET[channel]->GetMaximum()));
        hMCMET[channel]->GetXaxis()->SetRangeUser(0,300);

        TLegend *legend = new TLegend(0.55, 0.65, 0.91,0.92,"");      
        legend->AddEntry(hResultMET[channel], "Data-Driven Est", "PE");
        legend->AddEntry(hzjets, "Z+Jets MC", "F");
        legend->AddEntry(httbar, "TTbar MC", "F");
        
        legend->SetTextSize(0.05);
        legend->SetTextFont(42);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->Draw();
        
        hResultMET[channel]->Write();

      }//channel
      cout<<" All Channel Total is \n"<<resultMETCombined<<endl;
      c1.Print("FakeRateMET.pdf");
      fMM->Save();
      
    }//sample loop
  }else{//not EWK
    TTree* tEvts = new TTree("tEvts", "Cut Values per sample");
    tEvts->ReadFile("../Limits/cutValues.wz.dat");
    tEvts->Draw("SignalCode:Mass:minWindow:maxWindow:HtCut:ZptCut:WptCut","", "para goff");
    double ntEvts = tEvts->GetSelectedRows(); cout<<"Found "<<ntEvts<<endl;

    bool useHt(1), useZpt(0), useWpt(0), useWindow(1);

    for(int isignal=0; isignal<ntEvts; ++isignal){
      const int SignalCode = tEvts->GetVal(0)[isignal];
      const string SignalName = SampleName(SignalCode)[0];
      const double mass = tEvts->GetVal(1)[isignal];
      const double minWindow = useWindow ? tEvts->GetVal(2)[isignal] : 0;
      const double maxWindow = useWindow ? tEvts->GetVal(3)[isignal] : 99999;
      const double minHt     = useHt     ? tEvts->GetVal(4)[isignal] : 0;
      const double minZpt    = useZpt    ? tEvts->GetVal(5)[isignal] : 0;
      const double minWpt    = useWpt    ? tEvts->GetVal(6)[isignal] : 0;
      
      
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
  double e3Tight = hTight->GetBinContent(1);
  double e2Tight = hTight->GetBinContent(2);
  double e1Tight = hTight->GetBinContent(3);
  double e0Tight = hTight->GetBinContent(4);  
  double allTight = e3Tight + e2Tight + e1Tight + e0Tight;
  
  double Se3Tight = sqrt(e3Tight);
  double Se2Tight = sqrt(e2Tight);
  double Se1Tight = sqrt(e1Tight);
  double Se0Tight = sqrt(e0Tight);
  double SallTight = sqrt(allTight);
  
  cout<<"Tight Sample: "<<" N total: "<<Value(allTight,SallTight)<<" +/- "<<Value(SallTight)<<endl;

  if(0) return;

  if(allChannels)
    cout<<"N 3e: "<<e3Tight<<" +/- "<<Se3Tight<<endl//" frac: "<<e3/total<<endl
        <<"N 2e: "<<e2Tight<<" +/- "<<Se2Tight<<endl//" frac: "<<e2/total<<endl
        <<"N 1e: "<<e1Tight<<" +/- "<<Se1Tight<<endl//" frac: "<<e1/total<<endl
        <<"N 0e: "<<e0Tight<<" +/- "<<Se0Tight<<endl//" frac: "<<e0/total<<endl
        <<endl;

  double e3Loose = hLoose->GetBinContent(1);
  double e2Loose = hLoose->GetBinContent(2);
  double e1Loose = hLoose->GetBinContent(3);
  double e0Loose = hLoose->GetBinContent(4);    
  double allLoose = e3Loose + e2Loose + e1Loose + e0Loose;
  
  double Se3Loose = sqrt(e3Loose);
  double Se2Loose = sqrt(e2Loose);
  double Se1Loose = sqrt(e1Loose);
  double Se0Loose = sqrt(e0Loose);
  double SallLoose = sqrt(allLoose);
  
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
CalcMatrix(const double eTight, const  double Delta_eTight, 
           const double pFake, const double Delta_pFake, 
           const double Nl, const double Delta_Nl,
           const double Nt, const double Delta_Nt){
    MatrixResult result;

    double N1 = Nl - Nt;//number of failed events
    double N2 = Nt;

    double dN1 = TMath::Sqrt(N1);//Vuko's improvement
    //double dN1 = TMath::Sqrt(Delta_Nl*Delta_Nl + Delta_Nt*Delta_Nt);  
    double dN2 = Delta_Nt;
    

    double Nlep = (Nt - pFake*Nl)/(eTight - pFake);
    double Njet = (eTight*Nl - Nt)/(eTight - pFake);
    //cout << "Loose: N_lep: " << Nlep << " and: N_jet: " << Njet << endl;
  
    result.Nt_lep = eTight*Nlep;
    result.Nt_jet = pFake*Njet;
    
    //cout << "Tight: N_lep: " << Nt_lep << " and: N_jet: " << Nt_jet << endl;  
    
    //errors
    
    double dNlep_desig = (pFake*(pFake*(N1+N2)-N2))/((eTight-pFake)*(eTight-pFake));
    
    double dNlep_deqcd = (eTight*(N2-eTight*(N1+N2)))/((eTight-pFake)*(eTight-pFake));   
    
    double dNlep_dN1 = (-eTight*pFake)/(eTight-pFake);
    
    double dNlep_dN2 = (eTight*(1-pFake))/(eTight-pFake);


    double dNjet_desig =  (pFake*(N2-pFake*(N1+N2)))/((eTight-pFake)*(eTight-pFake)); 

    double dNjet_deqcd = (eTight*(eTight*(N1+N2)-N2))/((eTight-pFake)*(eTight-pFake));

    double dNjet_dN1 = (pFake*eTight)/(eTight-pFake);

    double dNjet_dN2 = (pFake*(eTight-1))/(eTight-pFake);


    result.DeltaNt_lep = TMath::Sqrt((dNlep_desig*dNlep_desig * Delta_eTight*Delta_eTight) + 
                                     (dNlep_deqcd*dNlep_deqcd * Delta_pFake *Delta_pFake ) + 
                                     (dNlep_dN1  *dNlep_dN1   * dN1         *dN1         ) + 
                                     (dNlep_dN2  *dNlep_dN2   * dN2         *dN2         )); 

    result.DeltaNt_jet = TMath::Sqrt((dNjet_desig*dNjet_desig * Delta_eTight*Delta_eTight) + 
                                     (dNjet_deqcd*dNjet_deqcd * Delta_pFake *Delta_pFake ) +
                                     (dNjet_dN1  *dNjet_dN1   * dN1         *dN1         ) +
                                     (dNjet_dN2  *dNjet_dN2   * dN2         *dN2         ));
  
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

MatrixSomething
CalcSomething(const Rates & wRates, const Rates & z1Rates, const Rates & z2Rates,
              const Value & vLoose, const Value & vTight, const Value & vTT, const Value & vTTT){
  MatrixSomething result;
  
  //Form Independent Variables
  double Nf = vLoose.val - vTight.val;
  double Delta_Nf = sqrt(Nf);
  double Ntf = vTight.val - vTT.val;
  double Delta_Ntf = sqrt(Ntf);
  double Nttf = vTT.val - vTTT.val;
  double Delta_Nttf = sqrt(Nttf);
  
  const double & e1 = wRates.Reff;
  const double & e2 = z1Rates.Reff;
  const double & e3 = z2Rates.Reff;
  const double & p1 = wRates.Rfake;
  const double & p2 = z1Rates.Rfake;
  const double & p3 = z2Rates.Rfake;

  const double & Nttt = vTTT.val;
  const double & Delta_Nttt = vTTT.err;

  //Determine Value
  result.NGood = 
    p1 * (Ntf + Nttf + Nttt - e1*(Nf+Ntf+Nttf+Nttt))/(e1 - p1) + 
    (Nttf + Nttt - p2*(Ntf + Nttf + Nttt))/(e2 - p2) + 
    (Nttt - e3*(Nttf + Nttt))/(e2*(e3-p3));
  result.NBad  = vTight.val - result.NGood;
  
  //Determine Error (NTTT, NTTF, NTF, NF, e1, e2, e3, p1, p2, p3)
  double ddNttt = (p1 - e1*p1)/(e1-p1) + (p2-1)/(p2-e2) + (e3-1)/(e2*(p3-e3));
  double ddNttf = (p1 - e1*p1)/(e1-p1) + (p2-1)/(p2-e2) + (e3  )/(e2*(p3-e3));
  double ddNtf = e1*(1-p1)/(e1-p1) + (e2)/(p2-e2);
  double ddNf = e1 * p1 / (e1-p1);
  double dde1 = p1 * (Ntf+Nttf+Nttt - p1*(Nf+Ntf+Nttf+Nttt)) / pow(e1-p1,2);
  double dde2 = (Nttf+Nttt - p2*(Ntf+Nttf+Nttt))/pow(e2-p2,2) + (Nttt-e3*(Nttf+Nttt))/(e2*e2*(e3-p3));
  double dde3 = (Nttt*(p3-1)+p3*Nttf)/(e2*pow(e3-p3,2));
  double ddp1 = e1*(Ntf+Nttf+Nttt-e1*(Nf+Ntf+Nttf+Nttt))/pow(e1-p1,2);
  double ddp2 = (Nttf+Nttt-e2*(Ntf+Nttf+Nttt))/pow(e2-p2,2);
  double ddp3 = (Nttt-e3*(Nttf+Nttt))/(e2*pow(e3-p3,2));
                
  result.DeltaNGood = sqrt( pow(ddNttt*Delta_Nttt, 2) +
                            pow(ddNttf*Delta_Nttf, 2) +
                            pow(ddNtf*Delta_Ntf, 2) +
                            pow(ddNf*Delta_Nf, 2) +
                            pow(dde1*wRates.Delta_Reff, 2) +
                            pow(dde2*z1Rates.Delta_Reff, 2) +
                            pow(dde3*z2Rates.Delta_Reff, 2) +
                            pow(ddp1*wRates.Delta_Rfake, 2) +
                            pow(ddp2*z1Rates.Delta_Rfake, 2) +
                            pow(ddp3*z2Rates.Delta_Rfake, 2) );
  result.DeltaNGood = std::min(result.DeltaNGood, (double)vTight.val);//Cory: Temp hack since we know it can't be bigger than this
  
  printf("Nttt=%.1f Nttf=%.1f Ntf=%.1f Nf=%.1f e1=%.1f e2=%.1f e3=%.1f p1=%.1f p2=%.1f p3=%.1f\n", 
         ddNttt, ddNttf, ddNtf, ddNf, dde1, dde2, dde3, ddp1, ddp2, ddp3);
  printf("Nttt=%.1f Nttf=%.1f Ntf=%.1f Nf=%.1f e1=%.1f e2=%.1f e3=%.1f p1=%.1f p2=%.1f p3=%.1f\n", 
         Delta_Nttt, Delta_Nttf, Delta_Ntf, Delta_Nf, 
         wRates.Delta_Reff, z1Rates.Delta_Reff, z2Rates.Delta_Reff,
         wRates.Delta_Rfake, z1Rates.Delta_Rfake, z2Rates.Delta_Rfake);
         
  result.DeltaNBad = 0.;//Cory: Fix

  return result;
                }

/*
MatrixSomething
CalcSomething(const double Nt, const double Delta_Nt,
              const double Ntt, const double Delta_Ntt,
              const Rates tRates, const Rates ttRates,
              const double Nf, const double Delta_Nf,
              const MatrixResult & matrix){
  MatrixSomething result;
  double Ntf = Nt - Ntt;//number of failed events
  double Delta_Ntf = sqrt(Ntf);
  double Nl = Nf + Ntt + Ntf;

  const double &      NtJet = matrix.Nt_jet;
  const double & DeltaNtJet = matrix.DeltaNt_jet;

  result.NGood = (Ntt - ttRates.Reff*ttRates.Rfake*(Ntf+Ntt) - ttRates.Reff*(ttRates.Reff-ttRates.Rfake)*NtJet) / (ttRates.Reff*(ttRates.Reff-ttRates.Rfake));
  result.NBad  = ( ttRates.Reff*ttRates.Reff*(Ntf+Ntt) - Ntt )  / (ttRates.Reff*(ttRates.Reff-ttRates.Rfake));

  double Tdiff = tRates.Rfake - tRates.Reff;
  double TTdiff = ttRates.Reff - ttRates.Rfake;
  double term1 = (1/ttRates.Reff - ttRates.Rfake)/TTdiff + (tRates.Rfake * (tRates.Reff-1))/Tdiff;
  double term2 = ttRates.Rfake/TTdiff + (tRates.Rfake * (1-tRates.Reff))/Tdiff;
  double term3 = tRates.Rfake*tRates.Reff/Tdiff;
  double term4 = (-2*ttRates.Reff*Ntt + ttRates.Rfake*Ntt+ttRates.Reff*ttRates.Reff*ttRates.Rfake*Nt)/pow(ttRates.Reff*TTdiff,2);
  double term5 = (Ntt - ttRates.Reff*ttRates.Reff*Nt)/(ttRates.Reff*TTdiff*TTdiff);
  double term6 = (tRates.Rfake*(Nt - Nl*tRates.Rfake))/pow(Tdiff,2);
  double term7 = (tRates.Reff*(Nt-Nl*tRates.Reff))/pow(Tdiff,2);
  result.DeltaNGood = sqrt( pow(term1*Delta_Ntt, 2) +
                            pow(term2*Delta_Ntf, 2) +
                            pow(term3*Delta_Nf , 2) +
                            pow(term4*ttRates.Delta_Reff, 2) +
                            pow(term5*ttRates.Delta_Rfake, 2) +
                            pow(term6*tRates.Delta_Reff, 2) +
                            pow(term7*tRates.Delta_Rfake, 2) );
  
  result.DeltaNBad = 0.;//Cory: Fix

  return result;
}
*/

double 
FindNEvts(TFile* f, vector<string>& samples, string hName, double pt1, double eta1, double pt2, double eta2, float& err){
  TH2F* hist = get_sum_of_hists2D(f, samples, hName, 0);
  //now get the integral of pt1->pt2, eta1->eta2

  return GetNEvts(hist, pt1, eta1, pt2, eta2, err);
}

double 
Closest(double value, double* array, int size){
  for(int i=0; i<size-1; ++i){
    if(array[i] <= value  && value < array[i+1]) return i;
  }
  cout<<" Value "<<value<<" is not found in the array!!!!\n";
  abort();
}



double
GetNEvts(TH2F* h, double xmin, double ymin, double xmax, double ymax, float& err){
  double total(0);
  double error(0);

  int xminbin,xmaxbin,yminbin,ymaxbin;
  xminbin = h->GetXaxis()->FindBin(xmin);
  xmaxbin = h->GetXaxis()->FindBin(xmax)-1;
  yminbin = h->GetYaxis()->FindBin(ymin);
  ymaxbin = h->GetYaxis()->FindBin(ymax);
  
  total  = h->IntegralAndError(xminbin, xmaxbin, yminbin, ymaxbin, error);
  err = error;
  
  //Cory: Deal with negative eta (Y axis)
  xminbin = h->GetXaxis()->FindBin(xmin);
  xmaxbin = h->GetXaxis()->FindBin(xmax)-1;
  yminbin = h->GetYaxis()->FindBin(-1*ymax);
  ymaxbin = h->GetYaxis()->FindBin(-1*ymin)-1;//avoid overlap
  
  error = 0;
  total  += h  ->IntegralAndError(xminbin, xmaxbin, yminbin, ymaxbin, error);
  err = sqrt(error*error + err*err);
  
  return total;
}


Value
  FindE(double pt, double eta, bool isElec, bool isW){
  if(isW){//T Values
    if(isElec){
      if(eta > 1.5){//EndCap
        if(pt > 50) return Value(0.9633,0.0007);
        if(pt > 30) return Value(0.9312,0.0003);
        if(pt > 25) return Value(0.8826,0.0017);
        if(pt > 20) return Value(0.8633,0.0027);
        else        return Value(0.8633,0.0027);
      }else{//Barrel
        if(pt > 50) return Value(0.9344,0.0018);
        if(pt > 30) return Value(0.8787,0.0008);
        if(pt > 25) return Value(0.8144,0.0032);
        if(pt > 20) return Value(0.7804,0.0042);
        else        return Value(0.7804,0.0042);
      }
    }else{//Muons
      if(pt > 50) return Value(0.9870,0.0003);
      if(pt > 30) return Value(0.9674,0.0002);
      if(pt > 25) return Value(0.9352,0.0010);
      if(pt > 20) return Value(0.9161,0.0015);
      else        return Value(0.9161,0.0015);
    }
  }else{//Z
    if(isElec){
      if(eta > 1.5){//EndCap
        if(pt > 50) return Value(0.8781, 0.0024 );
        if(pt > 30) return Value(0.8460, 0.0009 );
        if(pt > 25) return Value(0.7968, 0.0030 );
        if(pt > 20) return Value(0.7726, 0.0040 );
        else        return Value(0.7726, 0.0040 );
      }else{//Barrel
        if(pt > 50) return Value(0.7989, 0.0010 );
        if(pt > 30) return Value(0.7398, 0.0004 );
        if(pt > 25) return Value(0.6909, 0.0018 );
        if(pt > 20) return Value(0.6569, 0.0028 );
        else        return Value(0.6569, 0.0028 );
      }
    }else{//Muons
      if(pt > 50) return Value(0.9832, 0.0004 );
      if(pt > 30) return Value(0.9588, 0.0002 );
      if(pt > 25) return Value(0.9155, 0.0012 );
      if(pt > 20) return Value(0.9013, 0.0017 );
      else        return Value(0.9013, 0.0017 );
    }
  }
}

Value
  FindP(double pt, double eta, bool isElec, bool isW){
  if(isW){//T Values
    if(isElec){
      if(eta > 1.5){//EndCap
        if(pt > 50) return Value(0.0213,0.0281);
        if(pt > 30) return Value(0.0759,0.0792);
        if(pt > 25) return Value(0.2297,0.2338);
        if(pt > 20) return Value(0.1172,0.1205);
        else        return Value(0.0909,0.1066);
      }else{//Barrel
        if(pt > 50) return Value(0.0476,0.0528);
        if(pt > 30) return Value(0.0400,0.0424);
        if(pt > 25) return Value(0.0748,0.0781);
        if(pt > 20) return Value(0.0588,0.0609);
        else        return Value(0.0580,0.0659);
      }
    }else{//Muons
      if(pt > 50) return Value(0.1522,0.1611);
      if(pt > 30) return Value(0.0851,0.0885);
      if(pt > 25) return Value(0.0460,0.0524);
      if(pt > 20) return Value(0.0169,0.0204);
      else        return Value(0.0000,0.0255);
    }
  }else{//Z
    if(isElec){
      if(eta > 1.5){//EndCap
        if(pt > 50) return Value(0.0323, 0.0423);
        if(pt > 30) return Value(0.1038, 0.1080);
        if(pt > 25) return Value(0.2113, 0.2158);
        if(pt > 20) return Value(0.1546, 0.1586);
        else        return Value(0.0047, 0.0057);
      }else{//Barrel
        if(pt > 50) return Value(0.0471, 0.0536);
        if(pt > 30) return Value(0.0643, 0.0679);
        if(pt > 25) return Value(0.1571, 0.1628);
        if(pt > 20) return Value(0.1333, 0.1373);
        else        return Value(0.0037, 0.0042);
      }
    }else{//Muons
      if(pt > 50) return Value(0.9400, 0.8976);
      if(pt > 30) return Value(0.8235, 0.7997);
      if(pt > 25) return Value(0.8333, 0.7592);
      if(pt > 20) return Value(0.8000, 0.7213);
      else        return Value(0.0000, 0.0184);
    }
  }
 }
