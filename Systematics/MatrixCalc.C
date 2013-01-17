//Usage: root -b -l -q 'MatrixCalc.C+(useEWK, doSystematics)'
#include <iostream>
#include <iomanip>

#include "../Limits/consts.h"
#include "THStack.h"

/*
  struct Rates: Eff and Fake Rate plus errors
  struct MatrixResult: Output of MM
  struct MatrixYields: Number of TTT, TTF, ... events
  struct MatrixRates: Rates for W and Z leptons
  struct MatrixContainer: MatrixYields + MatrixRates
  struct MatrixVariable: How to bin a variable for MM
*/

const bool useMyEff  = true;
const bool useMyFake = false;//true;

using namespace std;

const double eElec = 0.9359; const double SeElec = 0.0004;//done
const double pElec = 0.07/**(1. - 0.3381)*/; const double SpElec = 0.01;//done

const double eMuon = 0.9711; const double SeMuon = 0.0002;
const double pMuon = 0.06/**(1. - -1*0.3021)*/; const double SpMuon = 0.01;//done

const bool debug = false;
const int NChannels = 4;

struct Rates{
  double Reff;
  double Delta_Reff;
  double Rfake;
  double Delta_Rfake;
  Rates():Reff(0),Delta_Reff(0),Rfake(0),Delta_Rfake(0){}
  Rates(double e, double se, double f, double sf):Reff(e),Delta_Reff(se),Rfake(f),Delta_Rfake(sf){}
  Rates(Value e, Value p):Reff(e.val),Delta_Reff(e.err),Rfake(p.val),Delta_Rfake(p.err){}
  friend std::ostream& operator << (std::ostream &o, const Rates & r){
    o<<"e="<<r.Reff<<" +/-"<< r.Delta_Reff<<"p="<< r.Rfake<<" +/- "<<r.Delta_Rfake;
    return o;
  }
};

struct MatrixResult{
  double NGood;
  double DeltaNGood;
  double NBad;
  double DeltaNBad;
  
  MatrixResult():NGood(0),DeltaNGood(0),NBad(0),DeltaNBad(0){}
  MatrixResult(double ng, double sng, double nb, double snb):NGood(ng),DeltaNGood(sng),NBad(nb),DeltaNBad(snb){}
  
  MatrixResult operator+( const MatrixResult &rhs ) const{
    MatrixResult result = *this;     // Make a copy of myself.  Same as MatrixResult result(*this);
    result += rhs;            // Use += to add other to the copy.
    return result;
  }
  MatrixResult & operator+=(const MatrixResult &rhs) {
    NGood += rhs.NGood;
    DeltaNGood = sqrt(pow(DeltaNGood,2) + pow(rhs.DeltaNGood,2));
    NBad += rhs.NBad;
    DeltaNBad = sqrt(pow(DeltaNBad,2) + pow(rhs.DeltaNBad,2));
    return *this;
  }

  friend std::ostream& operator << (std::ostream &o, const MatrixResult & m){
    o<< setiosflags(ios::fixed) << setprecision(2);
    o<<"NT_Good: "<<m.NGood<<" +/- "<<m.DeltaNGood;
    o<<"\tNT_Bad : "<<m.NBad<<" +/- "<<m.DeltaNBad;
    o<<"\tNT_SUM : "<<m.NBad+m.NGood;
    //o<<"NT_Lep: "<<Value(m.NGood,m.DeltaNGood)<<" +/- "<<Value(m.DeltaNGood)<<endl;
    //o<<"NT_Jet: "<<Value(m.NBad,m.DeltaNBad)<<" +/- "<<Value(m.DeltaNBad)<<endl;
    return o;
  }
  
};

struct MatrixYields{
  Value vTTT, vTTF, vTFT, vTFF, vFTT, vFTF, vFFT, vFFF;
  MatrixYields():vTTT(), vTTF(), vTFT(), vTFF(), vFTT(), vFTF(), vFFT(), vFFF(){}
  MatrixYields(Value ttt, Value ttf, Value tft, Value tff,
               Value ftt, Value ftf, Value fft, Value fff):
    vTTT(ttt), vTTF(ttf), vTFT(tft), vTFF(tff), 
    vFTT(ftt), vFTF(ftf), vFFT(fft), vFFF(fff){}
  MatrixYields operator+( const MatrixYields &rhs ) const{
    MatrixYields result = *this;     // Make a copy of myself.  Same as MatrixYields result(*this);
    result += rhs;            // Use += to add other to the copy.
    return result;
  }
  MatrixYields & operator+=(const MatrixYields &rhs) {
    vTTT   += rhs.vTTT;
    vTTF   += rhs.vTTF;
    vTFT   += rhs.vTFT;
    vTFF   += rhs.vTFF;
    vFTT   += rhs.vFTT;
    vFTF   += rhs.vFTF;
    vFFT   += rhs.vFFT;
    vFFF   += rhs.vFFF;
    return *this;
  }

  friend std::ostream& operator << (std::ostream &o, const MatrixYields & m){
    o<< setiosflags(ios::fixed) << setprecision(2);
    o<<"TTT: "<<m.vTTT<<" TTF: "<<m.vTTF<<" TFT: "<<m.vTFT<<" TFF: "<<m.vTFF
     <<" FTT: "<<m.vFTT<<" FTF: "<<m.vFTF<<" FFT: "<<m.vFFT<<" FFF: "<<m.vFFF;
    return o;
  }
};

struct MatrixRates{
  Rates wRates, z1Rates, z2Rates;
  MatrixRates():wRates(), z1Rates(), z2Rates(){}
  MatrixRates(Rates w, Rates z1, Rates z2):wRates(w), z1Rates(z1), z2Rates(z2){}
};

struct MatrixContainer{
  MatrixYields yield;
  MatrixResult result;

  MatrixContainer operator+( const MatrixContainer &rhs ) const{
    MatrixContainer out = *this;     // Make a copy of myself.  Same as MatrixContainer result(*this);
    out += rhs;            // Use += to add other to the copy.
    return out;
  }
  MatrixContainer & operator+=(const MatrixContainer &rhs) {
    yield  += rhs.yield;
    result += rhs.result;
    return *this;
  }

  friend std::ostream& operator << (std::ostream &o, const MatrixContainer & m){
    o<< setiosflags(ios::fixed) << setprecision(2);
    o<<"Yields : "<<m.yield<<endl;
    o<<"Results: "<<m.result;
    return o;
  }
};
  
struct MatrixVariable{
  string variable;
  int nbins;
  float min, max;
  vector<float> bins;
  TH1F* hist;
  MatrixVariable(string v, int nb, float m):
    variable(v), nbins(nb), min(0), max(m), hist(NULL){
    hist = new TH1F(v.c_str(), v.c_str(), nbins, min, max);
    setBins();
  }
  MatrixVariable(string v, int nb, float minimum, float maximum):
    variable(v), nbins(nb), min(minimum), max(maximum), hist(NULL){
    hist = new TH1F(v.c_str(), v.c_str(), nbins, min, max);
    setBins();
  }
  MatrixVariable(string v, string file, string histo):
    variable(v){
    TFile* f = new TFile(file.c_str(), "read"); assert(f);
    hist = (TH1F*)f->Get(histo.c_str());
    //cout<<"getting hist "<<histo<<" from file "<<file<<endl;
    assert(hist);
    nbins = hist->GetNbinsX();
    min = hist->GetXaxis()->GetXmin();
    max = hist->GetXaxis()->GetXmax();
    setBins();
  }
  void setBins(){
    bins.assign(nbins+1, 0.);
    float step = (max-min) / nbins;
    float lowEdge = min;
    for(int i=0; i<nbins; ++i,lowEdge+=step){
      bins[i] = lowEdge;
    }
    bins[nbins] = 9999;
    nbins++;//Account for overflow
  }
};

double FindNEvts(TFile* f, vector<string>& samples, string hName, double pt1, double eta1, double pt2, double eta2, float &err);
double GetNEvts(TH2F* h, double xmin, double ymin, double xmax, double ymax, float& err);

Value FindE(double pt, double eta, bool isElec, bool isW, bool isMC);
Value FindP(double pt, double eta, bool isElec, bool isW, bool isMC, bool isSystematics=false);

void MatrixMethod(const TH1F *hLoose, const TH1F *hTight, bool allChannel=true);
MatrixResult CalcMatrix2x2(const double eTight, const  double Delta_eTight,
                           const double pFake, const double Delta_pFake, 
                           const double Nl, const double Delta_Nl,
                           const double Nt, const double Delta_Nt);
MatrixResult CalcMatrix2x2(const Rates & wRates, const Value & vl, const Value & vt);

MatrixResult CalcMatrix4x4(const Rates & wRates, const Rates & z1Rates, const Rates & z2Rates,
                           const Value & vLoose, const Value & vTight, const Value & vTT, const Value & vTTT);

MatrixResult CalcMatrix8x8(const Rates & wRates, const Rates & z1Rates, const Rates & z2Rates,
                           const Value & vTTT, const Value & vTTF, const Value & vTFT, const Value & vTFF,
                           const Value & vFTT, const Value & vFTF, const Value & vFFT, const Value & vFFF);
MatrixResult CalcMatrix8x8(const Rates & wRates, const Rates & z1Rates, const Rates & z2Rates,
                           const MatrixYields & yields);


void
PrintLatexTable(const Value vTTT[], const Value vTTF[], const Value vTFT[], const Value vTFF[],
                const Value vFTT[], const Value vFTF[], const Value vFFT[], const Value vFFF[],
                const MatrixResult result[], const string & sample);
void
PrintLatexTable(const MatrixYields vByChan[], const MatrixResult result[], const string & sample);
void PrintLatexTable(const MatrixContainer container[], const string & sample);


void PrintRateTable(bool printEff, bool isElec, bool isMC);

MatrixRates
GetRates(int channel, float wpt, float weta, float z1pt, float z1eta, float z2pt, float z2eta, bool isMC, bool doSystematics);
MatrixYields
GetYields(int WTightCode, int ZTightCode, float weight);
MatrixContainer
AnalyzeTree(TTree* tree, const string & cuts, bool isMC, bool doSystematics);
void Analyze(TTree* tree, const MatrixVariable & variable, const string & sample, const bool isMC, const bool doSystematics, const int & verbose=0, const bool & printNBad=true);
void Analyze(TTree* tree, const vector<MatrixVariable> & variables, const string & sample, const bool isMC, const bool doSystematics, const int & verbose=0, const bool & printNBad=true);
void PrintLatexTable(const MatrixVariable & var, const MatrixContainer containerBin[][NChannels], const string & sample, const bool & printNBad);
void WriteHisto(const MatrixVariable & var, const MatrixContainer containerBin[][NChannels]);

void MatrixCalc(bool useEWK=true, bool doSystematics=false){
  cout<<"useEWK: "<<useEWK<<" doSystematics: "<<doSystematics<<endl;

  string tightfilename, loosefilename;
  
  if(useEWK){
    tightfilename = "../../../EWKWZ.root";
    loosefilename = "../../../EWKWZMatrix.root";
  }else{
    tightfilename = "../../../WprimeWZ.root";
    loosefilename = "../../../WprimeWZMatrix.root";
  }

  TFile *fMatrix = TFile::Open(loosefilename.c_str(), "read"); assert(fMatrix);

  vector<string> vBkg, vMC, samples;

  vBkg.push_back("DYJetsToLL");
  vBkg.push_back("TTJets");
  vBkg.push_back("ZZ");
  vBkg.push_back("GVJets");  
  //vBkg.push_back("WJetsToLNu");
  vBkg.push_back("WWTo2L2Nu");

  if(!useEWK)  vBkg.push_back("WZJetsTo3LNu");

  vMC = vBkg;

  if(useEWK) vMC.push_back("WZJetsTo3LNu");

//Cory: commenting to only show data
//  samples = vMC;
//  if(!useEWK) samples.push_back("BKG");
//  samples.push_back("MC");
//  if(!useEWK) samples.push_back("EXO");
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
        //if(!(samples[i] == "data" || samples[i] == "WZJetsTo3LNu")) continue;//Cory: Remove
        names.push_back(samples[i]);
      }

      bool isMC = (samples[i] != "data");
      cout<<" Doing MM for sample "<<samples[i]<<"( MC= "<<isMC<<" )"<<endl;

      TTree* tValidW = getTree(fMatrix, names, "tEvts_ValidW"); assert(tValidW);
      TTree* tF = getTree(fMatrix, names, "tEvts_MET"); assert(tF);
    
      //By Bin
      Analyze(tF, MatrixVariable("", 1, 1), samples[i], isMC, doSystematics, 50);
      //Analyze(tF, MatrixVariable("MET", 3, 80., 350.), samples[i], isMC, doSystematics, 0);
      //Analyze(tF, MatrixVariable("Zpt", 3, 80., 350.), samples[i], isMC, doSystematics, 0);
      Analyze(tF, MatrixVariable("WCharge", 2, -2, 2), samples[i], isMC, doSystematics, 0);
      Analyze(tValidW, MatrixVariable("MET", "../../../EWKWZMatrix.root", "data/hMET_ValidW"), samples[i], isMC, doSystematics, 0);
      Analyze(tF, MatrixVariable("Zpt", "../../../EWKWZMatrix.root", "data/hZpt_MET"), samples[i], isMC, doSystematics, 0);
      Analyze(tF, MatrixVariable("Wpt", "../../../EWKWZMatrix.root", "data/hWpt_MET"), samples[i], isMC, doSystematics, 0);
      Analyze(tF, MatrixVariable("ZMass", "../../../EWKWZMatrix.root", "data/hZMass_MET"), samples[i], isMC, doSystematics, 0);
      Analyze(tF, MatrixVariable("WTransMass", "../../../EWKWZMatrix.root", "data/hWTransMass_MET"), samples[i], isMC, doSystematics, 0);
      Analyze(tF, MatrixVariable("WLepPt", 30, 0, 300), samples[i], isMC, doSystematics, 0);
      Analyze(tF, MatrixVariable("NJets", "../../../EWKWZMatrix.root", "data/hNJets_MET"), samples[i], isMC, doSystematics, 0);

      if(0){
        for(int ch=3; ch<4; ++ch){
          for(int wcode=0; wcode<2; ++wcode){
            for(int zcode=0; zcode<4; ++zcode){
              cout<<"=======Ch: "<<ch<<" W Tight: "<<wcode<<" Z Tight: "<<zcode<<" ======="<<endl;
              tF->Draw("Run:Lumi:Event", Form("WTightCode==%i && ZTightCode==%i && EvtType==%i",wcode,zcode,ch), "para goff");
              int n = tF->GetSelectedRows();
              //loop over tree
              for(int ientry=0; ientry<n; ++ientry){
                int idx = 0;
                uint run    = round(tF->GetVal(idx++)[ientry]);
                uint lumi   = round(tF->GetVal(idx++)[ientry]);
                uint event  = round(tF->GetVal(idx++)[ientry]);
                cout<<"run: "<<run<<" lumi: "<<lumi<<" event: "<<event<<endl;
              }
            }
          }
        }
      }//print event lists?
      //tF->Draw("Run:Lumi:Event", Form("Event==1238207696"), "para goff");
      //int n = tF->GetSelectedRows();
      //cout<<"Found "<<n<<" events"<<endl;
      
    }//sample loop
  }else{//not EWK
    TFile *fTight = TFile::Open(tightfilename.c_str(), "read"); assert(fTight);
    TFile *fLoose = TFile::Open(loosefilename.c_str(), "read"); assert(fLoose);


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
  }//EWK IF
}//MatrixCalc main fn

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
  cout<<" For ByChan Channels\n";
  CalcMatrix2x2(eElec, SeElec,//Cory: This is wrong.  Have to weight elec & muon eff/fake rates
                pElec, SpElec,
                allLoose, SallLoose,
                allTight, SallTight);

  if(!allChannels) return;

  cout<<" For 3e\n";
  MatrixResult result2x23e = CalcMatrix2x2(eElec, SeElec,
                                           pElec, SpElec, 
                                           e3Loose, Se3Loose,
                                           e3Tight, Se3Tight);
  cout<<" For 2e\n";
  MatrixResult result2x22e = CalcMatrix2x2(eMuon, SeMuon,
                                           pMuon, SpMuon, 
                                           e2Loose, Se2Loose,
                                           e2Tight, Se2Tight);
  cout<<" For 1e\n";
  MatrixResult result2x21e = CalcMatrix2x2(eElec, SeElec,
                                           pElec, SpElec, 
                                           e1Loose, Se1Loose,
                                           e1Tight, Se1Tight);
  cout<<" For 0e\n";
  MatrixResult result2x20e = CalcMatrix2x2(eMuon, SeMuon,
                                           pMuon, SpMuon, 
                                           e0Loose, Se0Loose,
                                           e0Tight, Se0Tight);

  MatrixResult result2x2ByChan;
  result2x2ByChan.NGood = result2x23e.NGood
    + result2x22e.NGood
    + result2x21e.NGood 
    + result2x20e.NGood;

  result2x2ByChan.NBad = result2x23e.NBad +
    result2x22e.NBad +
    result2x21e.NBad +
    result2x20e.NBad;

  result2x2ByChan.DeltaNGood = sqrt(result2x23e.DeltaNGood*result2x23e.DeltaNGood +
                                    result2x22e.DeltaNGood*result2x22e.DeltaNGood +
                                    result2x21e.DeltaNGood*result2x21e.DeltaNGood +
                                    result2x20e.DeltaNGood*result2x20e.DeltaNGood);

  result2x2ByChan.DeltaNBad = sqrt(result2x23e.DeltaNBad*result2x23e.DeltaNBad +
                                   result2x22e.DeltaNBad*result2x22e.DeltaNBad +
                                   result2x21e.DeltaNBad*result2x21e.DeltaNBad +
                                   result2x20e.DeltaNBad*result2x20e.DeltaNBad);

  cout<<" & Tight nLep +- error & Tight nJet +- error "<<endl
      <<" & "<<Value(result2x2ByChan.NGood,result2x2ByChan.DeltaNGood)
      <<" & "<<Value(result2x2ByChan.DeltaNGood)
      <<" & "<<Value(result2x2ByChan.NBad,result2x2ByChan.DeltaNBad)
      <<" & "<<Value(result2x2ByChan.DeltaNBad)
      <<" &  & \\\\ \\hline"<<endl;
}

MatrixResult
CalcMatrix2x2(const Rates & wRates, const Value & vl, const Value & vt){
  return CalcMatrix2x2(wRates.Reff, wRates.Delta_Reff, 
                       wRates.Rfake, wRates.Delta_Rfake,
                       vl.val, vl.err, vt.val, vt.err);
}

MatrixResult
CalcMatrix2x2(const double eTight, const  double Delta_eTight, 
              const double pFake, const double Delta_pFake, 
              const double Nl, const double Delta_Nl,
              const double Nt, const double Delta_Nt){
  MatrixResult result2x2;

  double N1 = Nl - Nt;//number of failed events
  double N2 = Nt;

  double dN1 = TMath::Sqrt(N1);//Vuko's improvement
  //double dN1 = TMath::Sqrt(Delta_Nl*Delta_Nl + Delta_Nt*Delta_Nt);  
  double dN2 = Delta_Nt;
    

  double Nlep = (Nt - pFake*Nl)/(eTight - pFake);
  double Njet = (eTight*Nl - Nt)/(eTight - pFake);
  //cout << "Loose: N_lep: " << Nlep << " and: N_jet: " << Njet << endl;
  
  result2x2.NGood = eTight*Nlep;
  result2x2.NBad = pFake*Njet;
    
  //cout << "Tight: N_lep: " << NGood << " and: N_jet: " << NBad << endl;  
    
  //errors
    
  double dNlep_desig = (pFake*(pFake*(N1+N2)-N2))/((eTight-pFake)*(eTight-pFake));
    
  double dNlep_deqcd = (eTight*(N2-eTight*(N1+N2)))/((eTight-pFake)*(eTight-pFake));   
    
  double dNlep_dN1 = (-eTight*pFake)/(eTight-pFake);
    
  double dNlep_dN2 = (eTight*(1-pFake))/(eTight-pFake);


  double dNjet_desig =  (pFake*(N2-pFake*(N1+N2)))/((eTight-pFake)*(eTight-pFake)); 

  double dNjet_deqcd = (eTight*(eTight*(N1+N2)-N2))/((eTight-pFake)*(eTight-pFake));

  double dNjet_dN1 = (pFake*eTight)/(eTight-pFake);

  double dNjet_dN2 = (pFake*(eTight-1))/(eTight-pFake);


  result2x2.DeltaNGood = TMath::Sqrt((dNlep_desig*dNlep_desig * Delta_eTight*Delta_eTight) + 
                                     (dNlep_deqcd*dNlep_deqcd * Delta_pFake *Delta_pFake ) + 
                                     (dNlep_dN1  *dNlep_dN1   * dN1         *dN1         ) + 
                                     (dNlep_dN2  *dNlep_dN2   * dN2         *dN2         )); 

  result2x2.DeltaNBad = TMath::Sqrt((dNjet_desig*dNjet_desig * Delta_eTight*Delta_eTight) + 
                                    (dNjet_deqcd*dNjet_deqcd * Delta_pFake *Delta_pFake ) +
                                    (dNjet_dN1  *dNjet_dN1   * dN1         *dN1         ) +
                                    (dNjet_dN2  *dNjet_dN2   * dN2         *dN2         ));
  
  //cout << "Error on tNlep: " << DeltaNGood << " Total: " << NGood << " +- " << DeltaNGood << endl;

  //cout << "Error on tNjet: " << DeltaNBad << " Total: " << NBad << " +- " << DeltaNBad << endl;
/*
  cout<<"Debug Result2x2s are: \n"
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
  <<"Tight nLep: "<<result2x2.NGood<<" +/- "<<result2x2.DeltaNGood
  <<" Tight nJet: "<<result2x2.NBad<<" +/- "<<result2x2.DeltaNBad
  <<endl;
*/
  return result2x2;
}

MatrixResult
CalcMatrix4x4(const Rates & wRates, const Rates & z1Rates, const Rates & z2Rates,
              const Value & vLoose, const Value & vTight, const Value & vTT, const Value & vTTT){
  MatrixResult result;
  
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

  //result.NGood = min((double)vTight.val, result.NGood);//Cory: Hack to test
  //result.NGood = max(0., result.NGood);//Cory: Hack to test
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
  if(result.DeltaNGood > vLoose.val) printf("Error too large %.2f vs %.2f\n", result.DeltaNGood, vLoose.val);
  //result.DeltaNGood = std::min(result.DeltaNGood, (double)vTight.val);//Cory: Temp hack since we know it can't be bigger than this
  
  //printf("Partial Deriv: Nttt=% .1f Nttf=% .1f Ntf=% .1f Nf=% .1f e1=% .1f e2=% .1f e3=% .1f p1=% .1f p2=% .1f p3=% .1f\n", 
  //ddNttt, ddNttf, ddNtf, ddNf, dde1, dde2, dde3, ddp1, ddp2, ddp3);
  //printf("Errors:        Nttt=% .1f Nttf=% .1f Ntf=% .1f Nf=% .1f e1=% .1f e2=% .1f e3=% .1f p1=% .1f p2=% .1f p3=% .1f\n", 
  //       Delta_Nttt, Delta_Nttf, Delta_Ntf, Delta_Nf, 
  //       wRates.Delta_Reff, z1Rates.Delta_Reff, z2Rates.Delta_Reff,
  //       wRates.Delta_Rfake, z1Rates.Delta_Rfake, z2Rates.Delta_Rfake);
         
  result.DeltaNBad = 0.;//Cory: Fix

  return result;
}


MatrixResult CalcMatrix8x8(const Rates & wRates, const Rates & z1Rates, const Rates & z2Rates,
                           const Value & vTTT, const Value & vTTF, const Value & vTFT, const Value & vTFF,
                           const Value & vFTT, const Value & vFTF, const Value & vFFT, const Value & vFFF){
  MatrixResult result;
  const double & e1 = wRates.Reff;
  const double & e2 = z1Rates.Reff;
  const double & e3 = z2Rates.Reff;
  const double & p1 = wRates.Rfake;
  const double & p2 = z1Rates.Rfake;
  const double & p3 = z2Rates.Rfake;
  
  //Form Independent Variable`s
  const double & Nttt = vTTT.val;
  const double & Nttf = vTTF.val;
  const double & Ntft = vTFT.val;
  const double & Ntff = vTFF.val;
  const double & Nftt = vFTT.val;
  const double & Nftf = vFTF.val;
  const double & Nfft = vFFT.val;
  const double & Nfff = vFFF.val;

  const double Nt = Nttt+Nttf+Ntft+Ntff;
  const double Nl = Nftt+Nftf+Nfft+Nfff + Nt;
  //Determine Value
  double DenomInv = 1. / ((e1 - p1)*(e2 - p2)*(e3 - p3));
  double preFactor = e1 * e2 * e3 * DenomInv;//Now include z lepton eff
/*
  double base = (   Nttt*(1 - p1)*(1 - p2)*(1 - p3)
  - Nttf*(1 - p1)*(1 - p2)*(    p3)
  - Ntft*(1 - p1)*(    p2)*(1 - p3)
  + Ntff*(1 - p1)*(    p2)*(    p3)    
  - Nftt*(    p1)*(1 - p2)*(1 - p3)
  + Nftf*(    p1)*(1 - p2)*(    p3)
  + Nfft*(    p1)*(    p2)*(1 - p3)
  - Nfff*(    p1)*(    p2)*(    p3)
  );
*/
  //The pattern is l->factor of p
  //               j->factor of e
  //               t->factor of (1 - X)
  //               f->factor of    X
  //               sign of term is (-1)^(#T + #L)?

  double base_Nlll = ( 0//Just to make the next lines line up
                       + Nttt*(1 - p1)*(1 - p2)*(1 - p3)
                       - Nttf*(1 - p1)*(1 - p2)*(    p3)
                       - Ntft*(1 - p1)*(    p2)*(1 - p3)
                       + Ntff*(1 - p1)*(    p2)*(    p3)    
                       - Nftt*(    p1)*(1 - p2)*(1 - p3)
                       + Nftf*(    p1)*(1 - p2)*(    p3)
                       + Nfft*(    p1)*(    p2)*(1 - p3)
                       - Nfff*(    p1)*(    p2)*(    p3)     
    );
  double base_Nllj = ( 0//Just to make the next lines line up
                       - Nttt*(1 - p1)*(1 - p2)*(1 - e3)
                       + Nttf*(1 - p1)*(1 - p2)*(    e3)
                       + Ntft*(1 - p1)*(    p2)*(1 - e3)
                       - Ntff*(1 - p1)*(    p2)*(    e3)    
                       + Nftt*(    p1)*(1 - p2)*(1 - e3)
                       - Nftf*(    p1)*(1 - p2)*(    e3)
                       - Nfft*(    p1)*(    p2)*(1 - e3)
                       + Nfff*(    p1)*(    p2)*(    e3)     );
  double base_Nljl = ( 0//Just to make the next lines line up
                       - Nttt*(1 - p1)*(1 - e2)*(1 - p3)
                       + Nttf*(1 - p1)*(1 - e2)*(    p3)
                       + Ntft*(1 - p1)*(    e2)*(1 - p3)
                       - Ntff*(1 - p1)*(    e2)*(    p3)    
                       + Nftt*(    p1)*(1 - e2)*(1 - p3)
                       - Nftf*(    p1)*(1 - e2)*(    p3)
                       - Nfft*(    p1)*(    e2)*(1 - p3)
                       + Nfff*(    p1)*(    e2)*(    p3)     );
  double base_Nljj = ( 0//Just to make the next lines line up
                       + Nttt*(1 - p1)*(1 - e2)*(1 - e3)
                       - Nttf*(1 - p1)*(1 - e2)*(    e3)
                       - Ntft*(1 - p1)*(    e2)*(1 - e3)
                       + Ntff*(1 - p1)*(    e2)*(    e3)    
                       - Nftt*(    p1)*(1 - e2)*(1 - e3)
                       + Nftf*(    p1)*(1 - e2)*(    e3)
                       + Nfft*(    p1)*(    e2)*(1 - e3)
                       - Nfff*(    p1)*(    e2)*(    e3)     );
  double base_Njll = ( 0//Just to make the next lines line up
                       - Nttt*(1 - e1)*(1 - p2)*(1 - p3)
                       + Nttf*(1 - e1)*(1 - p2)*(    p3)
                       + Ntft*(1 - e1)*(    p2)*(1 - p3)
                       - Ntff*(1 - e1)*(    p2)*(    p3)    
                       + Nftt*(    e1)*(1 - p2)*(1 - p3)
                       - Nftf*(    e1)*(1 - p2)*(    p3)
                       - Nfft*(    e1)*(    p2)*(1 - p3)
                       + Nfff*(    e1)*(    p2)*(    p3)     );
  double base_Njlj = ( 0//Just to make the next lines line up
                       + Nttt*(1 - e1)*(1 - p2)*(1 - e3)
                       - Nttf*(1 - e1)*(1 - p2)*(    e3)
                       - Ntft*(1 - e1)*(    p2)*(1 - e3)
                       + Ntff*(1 - e1)*(    p2)*(    e3)    
                       - Nftt*(    e1)*(1 - p2)*(1 - e3)
                       + Nftf*(    e1)*(1 - p2)*(    e3)
                       + Nfft*(    e1)*(    p2)*(1 - e3)
                       - Nfff*(    e1)*(    p2)*(    e3)     );
  double base_Njjl = ( 0//Just to make the next lines line up
                       + Nttt*(1 - e1)*(1 - e2)*(1 - p3)
                       - Nttf*(1 - e1)*(1 - e2)*(    p3)
                       - Ntft*(1 - e1)*(    e2)*(1 - p3)
                       + Ntff*(1 - e1)*(    e2)*(    p3)    
                       - Nftt*(    e1)*(1 - e2)*(1 - p3)
                       + Nftf*(    e1)*(1 - e2)*(    p3)
                       + Nfft*(    e1)*(    e2)*(1 - p3)
                       - Nfff*(    e1)*(    e2)*(    p3)     );
  double base_Njjj = ( 0//Just to make the next lines line up
                       - Nttt*(1 - e1)*(1 - e2)*(1 - e3)
                       + Nttf*(1 - e1)*(1 - e2)*(    e3)
                       + Ntft*(1 - e1)*(    e2)*(1 - e3)
                       - Ntff*(1 - e1)*(    e2)*(    e3)    
                       + Nftt*(    e1)*(1 - e2)*(1 - e3)
                       - Nftf*(    e1)*(1 - e2)*(    e3)
                       - Nfft*(    e1)*(    e2)*(1 - e3)
                       + Nfff*(    e1)*(    e2)*(    e3)     );

  double Nl_test = base_Nlll + base_Nllj + base_Nljl + base_Nljj +
    base_Njll + base_Njlj + base_Njjl + base_Njjj;
  Nl_test *= DenomInv;
  //assert(fabs(Nl_test - Nl) < 0.001);//Cory

  if(debug) printf("preFactor: %.2f base: %.2f\n", preFactor, base_Nlll);
  result.NGood = preFactor*base_Nlll;
  result.NBad  = Nttt - result.NGood;

  //Determine Error on Ngood
  double dGdNttt = preFactor*(1 - p1)*(1 - p2)*(1 - p3);
  double dGdNttf = preFactor*(1 - p1)*(1 - p2)*(    p3);
  double dGdNtft = preFactor*(1 - p1)*(    p2)*(1 - p3);
  double dGdNtff = preFactor*(1 - p1)*(    p2)*(    p3);
  double dGdNftt = preFactor*(    p1)*(1 - p2)*(1 - p3);
  double dGdNftf = preFactor*(    p1)*(1 - p2)*(    p3);
  double dGdNfft = preFactor*(    p1)*(    p2)*(1 - p3);
  double dGdNfff = preFactor*(    p1)*(    p2)*(    p3);

  double dGde1 = result.NGood * (1/e1 - 1/(e1 - p1)); 
  double dGde2 = result.NGood * (1/e2 - 1/(e2 - p2));
  double dGde3 = result.NGood * (1/e3 - 1/(e3 - p3));

  double dBasedp1 = ( - Nttt*(1 - p2)*(1 - p3)
                      + Nttf*(1 - p2)*(    p3)
                      + Ntft*(    p2)*(1 - p3)
                      - Ntff*(    p2)*(    p3)
                      - Nftt*(1 - p2)*(1 - p3)
                      + Nftf*(1 - p2)*(    p3)
                      + Nfft*(    p2)*(1 - p3)
                      - Nfff*(    p2)*(    p3)
    );
  double dGdp1 = preFactor*(1/(e1-p1) * base_Nlll + dBasedp1);
  double dBasedp2 = (   - Nttt*(1 - p1)*(1 - p3)
                        + Nttf*(1 - p1)*(    p3)
                        - Ntft*(1 - p1)*(1 - p3)
                        + Ntff*(1 - p1)*(    p3)    
                        + Nftt*(    p1)*(1 - p3)
                        - Nftf*(    p1)*(    p3)
                        + Nfft*(    p1)*(1 - p3)
                        - Nfff*(    p1)*(    p3)
    );
  double dGdp2 = preFactor*(1/(e2-p2) * base_Nlll + dBasedp2);
  double dBasedp3 = ( - Nttt*(1 - p1)*(1 - p2)
                      - Nttf*(1 - p1)*(1 - p2)
                      + Ntft*(1 - p1)*(    p2)
                      + Ntff*(1 - p1)*(    p2)
                      + Nftt*(    p1)*(1 - p2)
                      + Nftf*(    p1)*(1 - p2)
                      - Nfft*(    p1)*(    p2)
                      - Nfff*(    p1)*(    p2)
    );
  double dGdp3 = preFactor*(1/(e3-p3) * base_Nlll + dBasedp3);

  result.DeltaNGood = sqrt( pow(dGdNttt*vTTT.err, 2) +
                            pow(dGdNttf*vTTF.err, 2) +
                            pow(dGdNtft*vTFT.err, 2) +
                            pow(dGdNtff*vTFF.err, 2) +
                            pow(dGdNftt*vFTT.err, 2) +
                            pow(dGdNftf*vFTF.err, 2) +
                            pow(dGdNfft*vFFT.err, 2) +
                            pow(dGdNfff*vFFF.err, 2) +
                            pow(dGde1* wRates.Delta_Reff, 2) +
                            pow(dGde2*z1Rates.Delta_Reff, 2) +
                            pow(dGde3*z2Rates.Delta_Reff, 2) +
                            pow(dGdp1* wRates.Delta_Rfake, 2) +
                            pow(dGdp2*z1Rates.Delta_Rfake, 2) +
                            pow(dGdp3*z2Rates.Delta_Rfake, 2) );

  //if(debug) if(result.DeltaNGood > Nl) printf("Error too large %.2f vs %.2f\n", result.DeltaNGood, Nl);
  //result.DeltaNGood = std::min(result.DeltaNGood, Nttt);//Cory: Temp hack since we know it can't be bigger than this

  //result.DeltaNBad = sqrt(pow(vTTT.err, 2) + pow(result.DeltaNGood,2));//Cory: Not right, but upper limit on error
 

  //Determine Error on NBad
  double dBdNttt = 1 - preFactor*(1 - p1)*(1 - p2)*(1 - p3);
  double dBdNttf = preFactor*(1 - p1)*(1 - p2)*(    p3);
  double dBdNtft = preFactor*(1 - p1)*(    p2)*(1 - p3);
  double dBdNtff = preFactor*(1 - p1)*(    p2)*(    p3);
  double dBdNftt = preFactor*(    p1)*(1 - p2)*(1 - p3);
  double dBdNftf = preFactor*(    p1)*(1 - p2)*(    p3);
  double dBdNfft = preFactor*(    p1)*(    p2)*(1 - p3);
  double dBdNfff = preFactor*(    p1)*(    p2)*(    p3);

  double dBde1 = result.NBad * (1/e1 - 1/(e1 - p1)); 
  double dBde2 = result.NBad * (1/e2 - 1/(e2 - p2));
  double dBde3 = result.NBad * (1/e3 - 1/(e3 - p3));

  double dBdp1 = 0 - dGdp1;
  double dBdp2 = 0 - dGdp2;
  double dBdp3 = 0 - dGdp3;

  result.DeltaNBad = sqrt( pow(dBdNttt*vTTT.err, 2) +
                           pow(dBdNttf*vTTF.err, 2) +
                           pow(dBdNtft*vTFT.err, 2) +
                           pow(dBdNtff*vTFF.err, 2) +
                           pow(dBdNftt*vFTT.err, 2) +
                           pow(dBdNftf*vFTF.err, 2) +
                           pow(dBdNfft*vFFT.err, 2) +
                           pow(dBdNfff*vFFF.err, 2) +
                           pow(dBde1* wRates.Delta_Reff, 2) +
                           pow(dBde2*z1Rates.Delta_Reff, 2) +
                           pow(dBde3*z2Rates.Delta_Reff, 2) +
                           pow(dBdp1* wRates.Delta_Rfake, 2) +
                           pow(dBdp2*z1Rates.Delta_Rfake, 2) +
                           pow(dBdp3*z2Rates.Delta_Rfake, 2) );


  return result;
}

MatrixResult CalcMatrix8x8(const MatrixRates & rates, const MatrixYields & yields){
  return CalcMatrix8x8(rates.wRates, rates.z1Rates, rates.z2Rates,
                       yields.vTTT, yields.vTTF, yields.vTFT, yields.vTFF,
                       yields.vFTT, yields.vFTF, yields.vFFT, yields.vFFF);
}

MatrixResult CalcMatrix8x8(const Rates & wRates, const Rates & z1Rates, const Rates & z2Rates,
                           const MatrixYields & yields){
  return CalcMatrix8x8(wRates, z1Rates, z2Rates,
                       yields.vTTT, yields.vTTF, yields.vTFT, yields.vTFF,
                       yields.vFTT, yields.vFTF, yields.vFFT, yields.vFFF);
}

void
PrintLatexTable(const MatrixContainer container[], const string & sample){
  cout<<"\\begin{table}[tbh] \\begin{center}"<<endl;
  cout<<"\\begin{tabular}{c|cccc} \\hline"<<endl;
  cout<<"             &  $3e$     & $2e1\\mu$   &  $ 2\\mu 1e$     &    $3\\mu$   \\\\  \\hline"<<endl;
  cout<<"$\\epsilon_1 \\epsilon_2 \\epsilon_3 n_{\\ell \\ell \\ell} $          &  "
      <<container[0].result.NGood<<" $\\pm$ "<<container[0].result.DeltaNGood<<" &  "
      <<container[1].result.NGood<<" $\\pm$ "<<container[1].result.DeltaNGood<<" &  "
      <<container[2].result.NGood<<" $\\pm$ "<<container[2].result.DeltaNGood<<" &  "
      <<container[3].result.NGood<<" $\\pm$ "<<container[3].result.DeltaNGood<<"   \\\\ \\hline  "<<endl
      <<"$N_{TTT} - \\epsilon_1 \\epsilon_2 \\epsilon_3 n_{\\ell \\ell \\ell} $          &  "
      <<container[0].result.NBad<<" $\\pm$ "<<container[0].result.DeltaNBad<<" &  "
      <<container[1].result.NBad<<" $\\pm$ "<<container[1].result.DeltaNBad<<" &  "
      <<container[2].result.NBad<<" $\\pm$ "<<container[2].result.DeltaNBad<<" &  "
      <<container[3].result.NBad<<" $\\pm$ "<<container[3].result.DeltaNBad<<"   \\\\ \\hline  "<<endl;
  cout<<"$N_{TTT} $  &   "<<container[0].yield.vTTT<<"  &   "<<container[1].yield.vTTT<<"  &   "<<container[2].yield.vTTT<<"   &   "<<container[3].yield.vTTT<<"   \\\\ "<<endl;
  cout<<"$N_{TTF} $  &   "<<container[0].yield.vTTF<<"  &   "<<container[1].yield.vTTF<<"  &   "<<container[2].yield.vTTF<<"   &   "<<container[3].yield.vTTF<<"   \\\\ "<<endl;
  cout<<"$N_{TFT} $  &   "<<container[0].yield.vTFT<<"  &   "<<container[1].yield.vTFT<<"  &   "<<container[2].yield.vTFT<<"   &   "<<container[3].yield.vTFT<<"   \\\\ "<<endl;
  cout<<"$N_{TFF} $  &   "<<container[0].yield.vTFF<<"  &   "<<container[1].yield.vTFF<<"  &   "<<container[2].yield.vTFF<<"   &   "<<container[3].yield.vTFF<<"   \\\\ "<<endl;
  cout<<"$N_{FTT} $  &   "<<container[0].yield.vFTT<<"  &   "<<container[1].yield.vFTT<<"  &   "<<container[2].yield.vFTT<<"   &   "<<container[3].yield.vFTT<<"   \\\\ "<<endl;
  cout<<"$N_{FTF} $  &   "<<container[0].yield.vFTF<<"  &   "<<container[1].yield.vFTF<<"  &   "<<container[2].yield.vFTF<<"   &   "<<container[3].yield.vFTF<<"   \\\\ "<<endl;
  cout<<"$N_{FFT} $  &   "<<container[0].yield.vFFT<<"  &   "<<container[1].yield.vFFT<<"  &   "<<container[2].yield.vFFT<<"   &   "<<container[3].yield.vFFT<<"   \\\\ "<<endl;
  cout<<"$N_{FFF} $  &   "<<container[0].yield.vFFF<<"  &   "<<container[1].yield.vFFF<<"  &   "<<container[2].yield.vFFF<<"   &   "<<container[3].yield.vFFF<<"   \\\\ "<<endl;
  cout<<"\\hline \n \\end{tabular}"<<endl;
  cout<<"\\caption{Yields from "<<sample<<"}"<<endl;
  cout<<"\\label{tab:MM-"<<sample<<"}"<<endl;
  cout<<"\\end{center} \\end{table}"<<endl;
}

void
PrintLatexTable(const MatrixYields vByChan[], const MatrixResult result[], const string & sample){
  cout<<"\\begin{table}[tbh] \\begin{center}"<<endl;
  cout<<"\\begin{tabular}{c|cccc} \\hline \\hline"<<endl;
  cout<<"             &  $3e$     & $2e1\\mu$   &  $ 2\\mu 1e$     &    $3\\mu$   \\\\  \\hline"<<endl;
  cout<<"$\\epsilon_1 \\epsilon_2 \\epsilon_3 n_{\\ell \\ell \\ell} $          &  "
      <<result[0].NGood<<" $\\pm$ "<<result[0].DeltaNGood<<" &  "
      <<result[1].NGood<<" $\\pm$ "<<result[1].DeltaNGood<<" &  "
      <<result[2].NGood<<" $\\pm$ "<<result[2].DeltaNGood<<" &  "
      <<result[3].NGood<<" $\\pm$ "<<result[3].DeltaNGood<<"   \\\\ \\hline  "<<endl
      <<result[0].NBad<<" $\\pm$ "<<result[0].DeltaNBad<<" &  "
      <<result[1].NBad<<" $\\pm$ "<<result[1].DeltaNBad<<" &  "
      <<result[2].NBad<<" $\\pm$ "<<result[2].DeltaNBad<<" &  "
      <<result[3].NBad<<" $\\pm$ "<<result[3].DeltaNBad<<"   \\\\ \\hline  "<<endl;
  cout<<"$N_{TTT} $  &   "<<vByChan[0].vTTT<<"  &   "<<vByChan[1].vTTT<<"  &   "<<vByChan[2].vTTT<<"   &   "<<vByChan[3].vTTT<<"   \\\\ \\hline"<<endl;
  cout<<"$N_{TTF} $  &   "<<vByChan[0].vTTF<<"  &   "<<vByChan[1].vTTF<<"  &   "<<vByChan[2].vTTF<<"   &   "<<vByChan[3].vTTF<<"   \\\\ \\hline"<<endl;
  cout<<"$N_{TFT} $  &   "<<vByChan[0].vTFT<<"  &   "<<vByChan[1].vTFT<<"  &   "<<vByChan[2].vTFT<<"   &   "<<vByChan[3].vTFT<<"   \\\\ \\hline"<<endl;
  cout<<"$N_{TFF} $  &   "<<vByChan[0].vTFF<<"  &   "<<vByChan[1].vTFF<<"  &   "<<vByChan[2].vTFF<<"   &   "<<vByChan[3].vTFF<<"   \\\\ \\hline"<<endl;
  cout<<"$N_{FTT} $  &   "<<vByChan[0].vFTT<<"  &   "<<vByChan[1].vFTT<<"  &   "<<vByChan[2].vFTT<<"   &   "<<vByChan[3].vFTT<<"   \\\\ \\hline"<<endl;
  cout<<"$N_{FTF} $  &   "<<vByChan[0].vFTF<<"  &   "<<vByChan[1].vFTF<<"  &   "<<vByChan[2].vFTF<<"   &   "<<vByChan[3].vFTF<<"   \\\\ \\hline"<<endl;
  cout<<"$N_{FFT} $  &   "<<vByChan[0].vFFT<<"  &   "<<vByChan[1].vFFT<<"  &   "<<vByChan[2].vFFT<<"   &   "<<vByChan[3].vFFT<<"   \\\\ \\hline"<<endl;
  cout<<"$N_{FFF} $  &   "<<vByChan[0].vFFF<<"  &   "<<vByChan[1].vFFF<<"  &   "<<vByChan[2].vFFF<<"   &   "<<vByChan[3].vFFF<<"   \\\\ \\hline"<<endl;
  cout<<"\\hline \n \\end{tabular}"<<endl;
  cout<<"\\caption{Yields from "<<sample<<"}"<<endl;
  cout<<"\\label{tab:MM-"<<sample<<"}"<<endl;
  cout<<"\\end{center} \\end{table}"<<endl;
}

void
PrintLatexTable(const Value vTTT[], const Value vTTF[], const Value vTFT[], const Value vTFF[],
                const Value vFTT[], const Value vFTF[], const Value vFFT[], const Value vFFF[],
                const MatrixResult result[], const string & sample){
  cout<<"\\begin{table}[tbh]\n\\begin{center}"<<endl;
  cout<<"\\begin{tabular}{c|cccc} \\hline \\hline"<<endl;
  cout<<"             &  $3e$     & $2e1\\mu$   &  $ 2\\mu 1e$     &    $3\\mu$   \\\\  \\hline"<<endl;
  cout<<"$\\epsilon_1 \\epsilon_2 \\epsilon_3 n_{\\ell \\ell \\ell} $          &  "
      <<result[0].NGood<<" $\\pm$ "<<result[0].DeltaNGood<<" &  "
      <<result[1].NGood<<" $\\pm$ "<<result[1].DeltaNGood<<" &  "
      <<result[2].NGood<<" $\\pm$ "<<result[2].DeltaNGood<<" &  "
      <<result[3].NGood<<" $\\pm$ "<<result[3].DeltaNGood<<"   \\\\ \\hline  "<<endl;
  cout<<"$N_{TTT} $  &   "<<vTTT[0].val<<"  &   "<<vTTT[1].val<<"  &   "<<vTTT[2].val<<"   &   "<<vTTT[3].val<<"   \\\\ \\hline"<<endl;
  cout<<"$N_{TTF} $  &   "<<vTTF[0].val<<"  &   "<<vTTF[1].val<<"  &   "<<vTTF[2].val<<"   &   "<<vTTF[3].val<<"   \\\\ \\hline"<<endl;
  cout<<"$N_{TFT} $  &   "<<vTFT[0].val<<"  &   "<<vTFT[1].val<<"  &   "<<vTFT[2].val<<"   &   "<<vTFT[3].val<<"   \\\\ \\hline"<<endl;
  cout<<"$N_{TFF} $  &   "<<vTFF[0].val<<"  &   "<<vTFF[1].val<<"  &   "<<vTFF[2].val<<"   &   "<<vTFF[3].val<<"   \\\\ \\hline"<<endl;
  cout<<"$N_{FTT} $  &   "<<vFTT[0].val<<"  &   "<<vFTT[1].val<<"  &   "<<vFTT[2].val<<"   &   "<<vFTT[3].val<<"   \\\\ \\hline"<<endl;
  cout<<"$N_{FTF} $  &   "<<vFTF[0].val<<"  &   "<<vFTF[1].val<<"  &   "<<vFTF[2].val<<"   &   "<<vFTF[3].val<<"   \\\\ \\hline"<<endl;
  cout<<"$N_{FFT} $  &   "<<vFFT[0].val<<"  &   "<<vFFT[1].val<<"  &   "<<vFFT[2].val<<"   &   "<<vFFT[3].val<<"   \\\\ \\hline"<<endl;
  cout<<"$N_{FFF} $  &   "<<vFFF[0].val<<"  &   "<<vFFF[1].val<<"  &   "<<vFFF[2].val<<"   &   "<<vFFF[3].val<<"   \\\\ \\hline"<<endl;
  cout<<"\\hline \n \\end{tabular}"<<endl;
  cout<<"\\caption{Yields from "<<sample<<"}"<<endl;
  cout<<"\\label{tab:mmm-"<<sample<<"}"<<endl;
  cout<<"\\end{center}"<<endl;
  cout<<"\\end{table}"<<endl;
}


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
  
  //Deal with negative eta (Y axis)
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
FindE(double pt, double eta, bool isElec, bool isW, bool isMC){
  if(isMC){
    if(isW){//For W use T Values
      if(isElec){
        if(fabs(eta) > 1.479){//EndCap
          if(pt > 50) return Value(0.9651,0.0007);
          if(pt > 30) return Value(0.9321,0.0003);
          if(pt > 25) return Value(0.8913,0.0014);
          else        return Value(0.8660,0.0020);
        }else{//Barrel
          if(pt > 50) return Value(0.9833,0.0002);
          if(pt > 30) return Value(0.9676,0.0001);
          if(pt > 25) return Value(0.9408,0.0007);
          else        return Value(0.9262,0.0011);
        }
      }else{//Muons
        if(pt > 50) return Value(0.9751,0.0003);
        if(pt > 30) return Value(0.9331,0.0002);
        if(pt > 25) return Value(0.8619,0.0008);
        else        return Value(0.8074,0.0016);
      }
    }else{//For Z use TT values
      if(isElec){
        if(fabs(eta) > 1.479){//EndCap
          if(pt > 50) return Value(0.9914,0.0004);
          if(pt > 30) return Value(0.9784,0.0002);
          if(pt > 25) return Value(0.9599,0.0008);  
          if(pt > 20) return Value(0.9446,0.0011);
          else        return Value(0.9206,0.0020);
        }else{//Barrel
          if(pt > 50) return Value(0.9979,0.0001);
          if(pt > 30) return Value(0.9950,0.0001);   
          if(pt > 25) return Value(0.9885,0.0003);
          if(pt > 20) return Value(0.9850,0.0005);
          else        return Value(0.9800,0.0009);
        }
      }else{//Muons
        if(pt > 50) return Value(0.9938,0.0002);
        if(pt > 30) return Value(0.9814,0.0001);
        if(pt > 25) return Value(0.9622,0.0008);
        if(pt > 20) return Value(0.9361,0.0012);
        else        return Value(0.8856,0.0023);
      }
    }
  }else{//Data
    if(isElec){
      if(useMyEff){//My values
        if(fabs(eta) > 1.479){//EndCap 
          if(pt > 50) return Value(0.9387,0.0011);
          if(pt > 30) return Value(0.8980,0.0005);
          if(pt > 25) return Value(0.8406,0.0020);
          if(pt > 20) return Value(0.8058,0.0026);
          else        return Value(0.6048,0.0053);
        }else{
          if(pt > 50) return Value(0.9736,0.0004);
          if(pt > 30) return Value(0.9572,0.0002);
          if(pt > 25) return Value(0.9244,0.0009);
          if(pt > 20) return Value(0.9031,0.0014);
          else        return Value(0.8145,0.0033);
        }
      }else{//Group II
        return Value(1.,0.); //Cory: This is to test other groups numbers
      }
    }else{//Muons
      if(useMyEff){//My Values
        if(pt > 50) return Value(0.9792,0.0004);
        if(pt > 30) return Value(0.9489,0.0002);
        if(pt > 25) return Value(0.9059,0.0009);
        if(pt > 20) return Value(0.8670,0.0015);
        else        return Value(0.7778,0.0024);
      }else{//Group II
        return Value(1.,0.); //Cory: This is to test other groups numbers
      }
    }
  }
}


double
FindPRatio(double pt, double eta, bool isElec, bool isW){
  Value pMC   = FindP(pt, eta, isElec, isW, true, false);
  Value pData = FindP(pt, eta, isElec, isW, false, false);
  return pData.val / pMC.val;
}

Value
FindP(double pt, double eta, bool isElec, bool isW, bool isMC, bool doSystematic){
  if(doSystematic){//Use Wjets Fake Rates
    double R = FindPRatio(pt,eta,isElec,isW);
    if      (isElec && fabs(eta) >  1.479){//Elec Endcap
      if(pt > 15) return Value(R*0.1501, 0.0390); //channel 2, pt  15 eta 1.5, (mean, err)= (0.1501, 0.0390) (15.65 / 104.21) +0.0369 -0.0410
      if(pt > 10) return Value(R*0.1886, 0.0274); //channel 2, pt  10 eta 1.5, (mean, err)= (0.1886, 0.0274) (44.32 / 235.04) +0.0275 -0.0273
    }else if(isElec && fabs(eta) <= 1.479){//Elec Barrel
      if(pt > 15) return Value(R*0.0797, 0.0144); //channel 2, pt  15 eta 0.0, (mean, err)= (0.0797, 0.0144) (33.15 / 415.91) +0.0153 -0.0135
      if(pt > 10) return Value(R*0.1116, 0.0098); //channel 2, pt  10 eta 0.0, (mean, err)= (0.1116, 0.0098) (124.66 / 1116.67) +0.0096 -0.0100
    }else{//Muons
      if(pt > 15) return Value(R*0.0236, 0.0155); //channel 1, pt  15 eta 0.0, (mean, err)= (0.0236, 0.0155) (2.97 / 125.97) +0.0130 -0.0179
      if(pt > 10) return Value(R*0.0585, 0.0124); //channel 1, pt  10 eta 0.0, (mean, err)= (0.0585, 0.0124) (24.69 / 421.79) +0.0121 -0.0128
    }
  }//end ifSystematic
  if(isMC){
    if      (isElec && fabs(eta) >  1.479){//Elec Endcap
      if(pt > 15) return Value(0.1734, 0.0402); //channel 2, pt  15 eta 1.5, (mean, err)= (0.1734, 0.0402) (18.86 / 108.77) +0.0369 -0.0435
      if(pt > 10) return Value(0.1870, 0.0271); //channel 2, pt  10 eta 1.5, (mean, err)= (0.1870, 0.0271) (44.54 / 238.22) +0.0264 -0.0278
    }else if(isElec && fabs(eta) <= 1.479){//Elec Barrel
      if(pt > 15) return Value(0.0825, 0.0145); //channel 2, pt  15 eta 0.0, (mean, err)= (0.0825, 0.0145) (34.58 / 419.15) +0.0142 -0.0147
      if(pt > 10) return Value(0.1141, 0.0098); //channel 2, pt  10 eta 0.0, (mean, err)= (0.1141, 0.0098) (128.67 / 1127.30) +0.0096 -0.0100
    }else{//Muons
      if(pt > 15) return Value(0.0300, 0.0174); //channel 1, pt  15 eta 0.0, (mean, err)= (0.0300, 0.0174) (3.86 / 128.54) +0.0156 -0.0193
      if(pt > 10) return Value(0.0602, 0.0125); //channel 1, pt  10 eta 0.0, (mean, err)= (0.0602, 0.0125) (25.80 / 428.93) +0.0118 -0.0131
    }
  }else{//Data
    if(isElec){//These cam from min jet cut of 35 GeV
      if(useMyFake){
        if(fabs(eta) < 1.479){//barrel
          return Value(0.1379, 0.0500);//Remove pt binning bc of low stats
          if(pt > 25) return Value(0.2174, 0.1050); //channel 2, pt  25 eta 0.0, (mean, err)= (0.2174, 0.1050) (5.00 / 23.00) +0.1198 -0.0901
          if(pt > 15) return Value(0.1207, 0.0509); //channel 2, pt  15 eta 0.0, (mean, err)= (0.1207, 0.0509) (7.00 / 58.00) +0.0587 -0.0432
          if(pt > 10) return Value(0.1290, 0.0398); //channel 2, pt  10 eta 0.0, (mean, err)= (0.1290, 0.0398) (12.00 / 93.00) +0.0444 -0.0353
        }else{
          return Value(0.0810, 0.0500);//Remove pt binning bc of low stats
          if(pt > 25) return Value(0.1538, 0.1329); //channel 2, pt  25 eta 1.5, (mean, err)= (0.1538, 0.1329) (2.00 / 13.00) +0.1673 -0.0984
          if(pt > 15) return Value(0.0323, 0.0483); //channel 2, pt  15 eta 1.5, (mean, err)= (0.0323, 0.0483) (1.00 / 31.00) +0.0700 -0.0266
          if(pt > 10) return Value(0.1000, 0.0706); //channel 2, pt  10 eta 1.5, (mean, err)= (0.1000, 0.0706) (3.00 / 30.00) +0.0874 -0.0537
        }
      }else{//Group II
        if(fabs(eta) < 1.0){
          if(pt > 30) return Value(0.085,0.017);
          if(pt > 25) return Value(0.085,0.012);
          if(pt > 20) return Value(0.085,0.009);
          if(pt > 15) return Value(0.057,0.009);
          else        return Value(0.066,0.015);
        }else if(fabs(eta) < 1.479){
          if(pt > 30) return Value(0.059,0.019);
          if(pt > 25) return Value(0.067,0.015);
          if(pt > 20) return Value(0.064,0.011);
          if(pt > 15) return Value(0.054,0.010);
          else        return Value(0.040,0.012);
        }else if(fabs(eta) < 2.0){
          if(pt > 30) return Value(0.072,0.019);
          if(pt > 25) return Value(0.052,0.012);
          if(pt > 20) return Value(0.057,0.010);
          if(pt > 15) return Value(0.018,0.007);
          else        return Value(0.016,0.009);
        }else if(fabs(eta) < 2.5){
          if(pt > 30) return Value(0.068,0.021);
          if(pt > 25) return Value(0.070,0.015);
          if(pt > 20) return Value(0.046,0.009);
          if(pt > 15) return Value(0.041,0.013);
          else        return Value(0.023,0.013);
        }
      }
    }else{//Muons
      if(useMyFake){//Use my values
        if(pt > 25) return Value(0.0588, 0.0549); //channel 1, pt  25 eta 0.0, (mean, err)= (0.0588, 0.0549) (2.00 / 34.00) +0.0720 -0.0378
        if(pt > 15) return Value(0.0753, 0.0326); //channel 1, pt  15 eta 0.0, (mean, err)= (0.0753, 0.0326) (7.00 / 93.00) +0.0380 -0.0272
        if(pt > 10) return Value(0.1429, 0.0274); //channel 1, pt  10 eta 0.0, (mean, err)= (0.1429, 0.0274) (28.00 / 196.00) +0.0294 -0.0254
      }else{//Group II
        if(fabs(eta) < 1.0){
          if(pt > 30) return Value(0.034,0.007);
          if(pt > 25) return Value(0.008,0.004);
          if(pt > 20) return Value(0.011,0.001);
          if(pt > 15) return Value(0.010,0.001);
          else        return Value(0.016,0.001);
        }else if(fabs(eta) < 1.479){
          if(pt > 30) return Value(0.017,0.002);
          if(pt > 25) return Value(0.009,0.002);
          if(pt > 20) return Value(0.015,0.002);
          if(pt > 15) return Value(0.029,0.013);
          else        return Value(0.050,0.012);
        }else if(fabs(eta) < 2.0){
          if(pt > 30) return Value(0.026,0.002);
          if(pt > 25) return Value(0.016,0.003);
          if(pt > 20) return Value(0.025,0.003);
          if(pt > 15) return Value(0.062,0.022);
          else        return Value(0.057,0.015);
        }else if(fabs(eta) < 2.5){
          if(pt > 30) return Value(0.022,0.004);
          if(pt > 25) return Value(0.023,0.006);
          if(pt > 20) return Value(0.031,0.006);
          if(pt > 15) return Value(0.001,0.001);
          else        return Value(0.052,0.023);
        }
      }
    }//end muons
  }//end data
  abort();
  return Value();
}
  
MatrixRates
GetRates(int channel, float wpt, float weta, float z1pt, float z1eta, float z2pt, float z2eta, bool isMC, bool doSystematics){
  bool wElec = channel%2 == false;
  bool zElec = channel<=1;
  
  MatrixRates rates;
  rates.wRates = Rates(FindE( wpt,  weta, wElec, true , isMC), FindP( wpt,  weta, wElec,  true, isMC,doSystematics));
  rates.z1Rates= Rates(FindE(z1pt, z1eta, zElec, false, isMC), FindP(z1pt, z1eta, zElec, false, isMC,doSystematics));
  rates.z2Rates= Rates(FindE(z2pt, z2eta, zElec, false, isMC), FindP(z2pt, z2eta, zElec, false, isMC,doSystematics));
  
  if(debug){
    printf("Calculating with W : e=%.4f +- %.4f, p=%.4f +- %.4f\n", rates.wRates.Reff, rates.wRates.Delta_Reff, rates.wRates.Rfake, rates.wRates.Delta_Rfake);
    printf("Calculating with Z1: e=%.4f +- %.4f, p=%.4f +- %.4f\n", rates.z1Rates.Reff, rates.z1Rates.Delta_Reff, rates.z1Rates.Rfake, rates.z1Rates.Delta_Rfake);
    printf("Calculating with Z2: e=%.4f +- %.4f, p=%.4f +- %.4f\n", rates.z2Rates.Reff, rates.z2Rates.Delta_Reff, rates.z2Rates.Rfake, rates.z2Rates.Delta_Rfake);
  }

  return rates;
}

MatrixYields
GetYields(int WTightCode, int ZTightCode, float weight){
  float err = sqrt(weight);
  MatrixYields yields;
  yields.vTTT = Value(weight*(WTightCode==1 && ZTightCode==3), err*(WTightCode==1 && ZTightCode==3));//Final sample
  yields.vTTF = Value(weight*(WTightCode==1 && ZTightCode==2), err*(WTightCode==1 && ZTightCode==2));
  yields.vTFT = Value(weight*(WTightCode==1 && ZTightCode==1), err*(WTightCode==1 && ZTightCode==1));
  yields.vTFF = Value(weight*(WTightCode==1 && ZTightCode==0), err*(WTightCode==1 && ZTightCode==0));     
  yields.vFTT = Value(weight*(WTightCode==0 && ZTightCode==3), err*(WTightCode==0 && ZTightCode==3));
  yields.vFTF = Value(weight*(WTightCode==0 && ZTightCode==2), err*(WTightCode==0 && ZTightCode==2));
  yields.vFFT = Value(weight*(WTightCode==0 && ZTightCode==1), err*(WTightCode==0 && ZTightCode==1));
  yields.vFFF = Value(weight*(WTightCode==0 && ZTightCode==0), err*(WTightCode==0 && ZTightCode==0));

  if(debug) cout<<" Bin Event Total is ..."<<endl<<yields<<endl;          

  return yields;
}

MatrixContainer
AnalyzeTree(TTree* tree, const string & cuts, bool isMC, bool doSystematics){
  MatrixContainer container;
  tree->Draw("EvtType:WTightCode:WLepPt:WLepEta:ZTightCode:ZLep1Pt:ZLep1Eta:ZLep2Pt:ZLep2Eta", Form("weight*(%s)",cuts.c_str()), "para goff");
  int n = tree->GetSelectedRows();
  //loop over tree
  for(int ientry=0; ientry<n; ++ientry){
    float weight = tree->GetW()[ientry];
    int idx = 0;
    int channel    = round(tree->GetVal(idx++)[ientry]);
    int WTightCode = round(tree->GetVal(idx++)[ientry]);
    float wpt    = tree->GetVal(idx++)[ientry];
    float weta   = tree->GetVal(idx++)[ientry];
    int ZTightCode = round(tree->GetVal(idx++)[ientry]);
    float z1pt    = tree->GetVal(idx++)[ientry];
    float z1eta   = tree->GetVal(idx++)[ientry];
    float z2pt    = tree->GetVal(idx++)[ientry];
    float z2eta   = tree->GetVal(idx++)[ientry];
    
    MatrixRates rates = GetRates(channel, wpt, weta, z1pt, z1eta, z2pt, z2eta, isMC, doSystematics);
    MatrixYields yields = GetYields(WTightCode, ZTightCode, weight);
    
    MatrixResult tmpResult8x8 = CalcMatrix8x8(rates, yields);
    if(debug) cout<<" The event MMMM results are: "<<tmpResult8x8<<endl;

    container.yield += yields;
    container.result += tmpResult8x8;
    
  }
  if(debug) cout<<" The MMMM results are: "<<container.result<<endl;
  return container;
}
  
void
Analyze(TTree* tree, const MatrixVariable & variable, const string & sample, const bool isMC, const bool doSystematics, const int & verbose, const bool & printNBad){
  vector<MatrixVariable> variables(1, variable);
  Analyze(tree, variables, sample, isMC, doSystematics, verbose, printNBad);
}

void
Analyze(TTree* tree, const vector<MatrixVariable> & variables, const string & sample, const bool isMC, const bool doSystematics, const int & verbose, const bool & printNBad){
  if(verbose) cout<<"Verbose level "<<verbose<<endl;
  for(int iVar=0; iVar<(int)variables.size(); ++iVar){
    const MatrixVariable & var = variables[iVar];
    MatrixContainer containerBins[var.nbins][NChannels];
    MatrixContainer container[NChannels];
    //loop over bins
    for(int iBin=0; iBin<var.nbins-1; ++iBin){
      //loop over channels
      for(int channel=0; channel<NChannels; ++channel){
        string cuts = var.variable.empty() ? "(1)" : Form("(%s >= %.0f)*(%s < %.0f)", var.variable.c_str(), var.bins[iBin],  var.variable.c_str(), var.bins[iBin+1]);
        cuts += Form("*(EvtType==%i)", channel);
        //cout<<" Cuts are "<<cuts<<endl;
        MatrixContainer & chContainer = containerBins[iBin][channel];
        chContainer = AnalyzeTree(tree, cuts, isMC, doSystematics);
        container[channel] += chContainer;
        //cout<<"The bin results with cuts "<<cuts<<" are "<<chContainer.result<<endl;
      }//channel loop
    }//bin loop

    if(verbose > 10){
      PrintRateTable(true,  true,  isMC);
      PrintRateTable(false, true,  isMC);
      PrintRateTable(true,  false,  isMC);
      PrintRateTable(false, false,  isMC);
    }
    //if(verbose > 1) PrintLatexTable(var, containerBins, sample, printNBad);
    if(verbose > 1) PrintLatexTable(container, sample);

    if(var.hist) WriteHisto(var, containerBins);
    
    //if(var.variable.empty()) PrintLatexTable(vByChan, result8x8ByChan, samples[i]);


    if(verbose > 5){
      for(int channel=0; channel<NChannels; ++channel){
        //Print the channel sum by bin
        cout<<"Totals for channel "<<channel<<" are \n"<<container[channel]<<endl;
      }
    }
  }//var loop
}

void
PrintLatexTable(const MatrixVariable & var, const MatrixContainer containerBin[][NChannels], const string & sample, const bool & printNBad){
  cout<<"\\begin{table}[tbh] \\begin{center}"<<endl;
  cout<<"\\begin{tabular}{c|cccc} \\hline \\hline"<<endl;
  cout<<"   "<<var.variable<<"          &  $3e$     & $2e1\\mu$   &  $ 2\\mu 1e$     &    $3\\mu$   \\\\  \\hline"<<endl;
  for(int iBin=0; iBin<var.nbins-1; ++iBin){
    if(var.variable.empty()) printf(" Yield");
    else                     printf(" %.0f-%.0f",var.bins[iBin], var.bins[iBin+1]);
    for(int channel=0; channel<NChannels; ++channel){
      const MatrixResult & result = containerBin[iBin][channel].result;
      cout<<" & "<<result.NGood<<" $\\pm$ "<<result.DeltaNGood;
      if(printNBad) cout<<" ("<<result.NBad <<" $\\pm$ "<<result.DeltaNBad<<")";
    }//channel loop
    cout<<"  \\\\"<<endl;
  }//bin loop
  cout<<"\\hline \n \\end{tabular}"<<endl;
  cout<<"\\caption{"<<sample<<" signal (bkg) yields binned in "<<var.variable<<". NB: Last bin includes overflow}"<<endl;
  cout<<"\\label{tab:mmBinned-"<<var.variable<<"-"<<sample<<"}"<<endl;
  cout<<"\\end{center}"<<endl;
  cout<<"\\end{table}"<<endl;
}

void
WriteHisto(const MatrixVariable & var, const MatrixContainer containerBin[][NChannels]){
  //cout<<"writing histos for "<<var.hist->GetName()<<endl;
  if(var.variable.empty()) return;
  if(!var.hist){
    cout<<"aborting with variable = "<<var.variable
        <<"because no hist"<<endl;
    abort();
  }
  TFile* fMM = new TFile("DataDrivenHistos.root", "update");
  if(!fMM->GetKey("DYJetsToLL-DataDriven"))
    fMM->mkdir("DYJetsToLL-DataDriven");//change
  fMM->cd("DYJetsToLL-DataDriven");//change
  
  TH1F* hResult[NChannels] = {NULL};
  for(int channel=0; channel<NChannels; ++channel){
    string histName = var.hist->GetName();
    string chString = Form("%ie%im_", NChannels-1-channel, channel);
    if(histName.find("_") == string::npos) histName += chString;
    else histName.replace(histName.find("_"), 1, chString.c_str()); 
    hResult[channel] = (TH1F*) var.hist->Clone(histName.c_str());
    hResult[channel]->Reset();
    for(int iBin=0; iBin<var.nbins-1; ++iBin){
      const MatrixResult & result = containerBin[iBin][channel].result;
      int bin = hResult[channel]->FindBin(var.bins[iBin]);
      hResult[channel]->SetBinContent(bin, result.NBad);
      hResult[channel]->SetBinError  (bin, result.DeltaNBad);
    }
    hResult[channel]->Write();
  }
  fMM->Save();
  fMM->Close();
  //cout<<" done writing histos for "<<var.hist->GetName()<<endl;
}  


bool 
CheckDataEvent(float run, float event){
  bool shouldRemove = false;
  return shouldRemove;
}


void
PrintRateTable(bool printEff, bool isElec, bool isMC){
  float pt [] = {10., 15., 20., 25., 30., 9999}; const int npt  = 5;
  float eta[] = {0., 1., 1.479, 2., 2.5}; const int neta  = 4;

  cout<<"\\begin{table}[tbh] \\begin{center}"<<endl;
  cout<<"\\begin{tabular}{c|*{"<<neta<<"}{|c}} \\hline"<<endl;
  //.."\\begin{tabular}{|c||*{"<<levels.size()*evtTypes.size()<<"}{c|}}

  for(int iEta=0; iEta<neta; ++iEta){
    printf(" & %.3f $< |\\eta| \\leq$ %.3f",eta[iEta], eta[iEta+1]);
  }
  cout<<" \\\\  \\hline"<<endl;
  for(int iPt=0; iPt<npt; ++iPt){
    printf(" %.0f $< p_{T} \\leq$ %.0f",pt[iPt], pt[iPt+1]);
    for(int iEta=0; iEta<neta; ++iEta){
      const Rates rates(FindE( pt[iPt]+1.,  eta[iEta]+0.01, isElec, true , isMC),
                        FindP( pt[iPt]+1.,  eta[iEta]+0.01, isElec, true , isMC));
        
      if(printEff) printf(" & %.4f $\\pm$ %.4f", rates.Reff, rates.Delta_Reff);
      else         printf(" & %.4f $\\pm$ %.4f", rates.Rfake, rates.Delta_Rfake);
    }
    cout<<"  \\\\"<<endl;
  }//pt loop

  cout<<"\\hline \n \\end{tabular}"<<endl;
  printf("\\caption{%s %s Rate binned in $p_{T}$ and $\\eta$}\n", 
         isElec ? "Electron" : "Muon",
         printEff ? " Efficiency" : " Fake");
  cout<<"\\label{tab:fakeRates}"<<endl;
  cout<<"\\end{center} \\end{table}"<<endl;
}
