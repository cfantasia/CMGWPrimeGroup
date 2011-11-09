//Usage: root -b -q 'ExpectedEvts.C+("inName.root", window)'

#include "consts.h"
void AddDataDriven(double& Evts, double& sEvts, const double& ZJets, const double&sZJets);

bool debug_ = false;
bool useHists_ = false;

void
ExpectedEvts(string inName, string config, int windFracTenths=-1, string opt=""){
  if(opt.find("debug") != string::npos) debug_ = true;
  if(opt.find("useHists") != string::npos) useHists_ = true;

  gErrorIgnoreLevel = kWarning;
  double windFrac = windFracTenths/10.;

  TFile *f = TFile::Open(inName.c_str(), "read");
  TCanvas* c1 = new TCanvas("c1", "Number of Events");
  string outfile("nEvents.txt");
  ofstream out(outfile.c_str());
  if(!out) { 
    cout << "Cannot open file " << outfile << endl; 
    abort();
  } 

  out<<"SignalCode/I:"
     <<"Mass/F:"
     <<"Lumi/F:"
     <<"DataEvts/F:"
     <<"MCEvts/F:"
     <<"statMCEvts/F:"
     <<"sysMCEvts/F:"
     <<"sMCEvts/F:"
     <<"BkgEvts/F:"
     <<"statBkgEvts/F:"
     <<"sysBkgEvts/F:"
     <<"sBkgEvts/F:"
     <<"Eff/F:"
     <<"statEff/F:"
     <<"sysEff/F:"
     <<"sEff/F"
     <<endl;
  out  << setiosflags(ios::fixed) << setprecision(4);

  vector<string> BkgSamples; 
  vector<string> dataSamples;

  string paramString = "SignalCode:Mass:minWindow:maxWindow";
  string treeName, histName, varName;
  if(inName.find("WprimeWZ") != string::npos){
    BkgSamples.push_back("WWTo2L2Nu");
    BkgSamples.push_back("WJetsToLNu");
    BkgSamples.push_back("GVJets");  
    BkgSamples.push_back("ZZ"); 
    BkgSamples.push_back("TTJets"); 
    BkgSamples.push_back("DYJetsToLL"); 
    BkgSamples.push_back("WZJetsTo3LNu");

    dataSamples.push_back("data");

    paramString += ":HtCut:ZptCut:WptCut:ZJets:sZJets";
    treeName = "tWZCand";
    histName = "hWZMass_AllCuts";
    varName = "WZMass";
  }else if(inName.find("HadVZ") != string::npos){
    BkgSamples.push_back("Summer11_ZZJets_2l2q");
    BkgSamples.push_back("Summer11_VGamma");
    BkgSamples.push_back("Summer11_WW");
    BkgSamples.push_back("Summer11_WZ");
    BkgSamples.push_back("Summer11_TTJets");
    BkgSamples.push_back("Summer11_DYJetsToLL_PtZ100");

    dataSamples.push_back("data_DoubleMu-Run2011A-May10ReReco-v1");
    dataSamples.push_back("data_SingleMu-Run2011A-May10ReReco-v1");
    dataSamples.push_back("data_DoubleElectron-Run2011A-May10ReReco-v1");
    dataSamples.push_back("data_SingleElectron-Run2011A-May10ReReco-v1");

    dataSamples.push_back("data_DoubleMu-Run2011A-PromptReco-v4");
    dataSamples.push_back("data_SingleMu-Run2011A-PromptReco-v4");
    dataSamples.push_back("data_DoubleElectron-Run2011A-PromptReco-v4");
    dataSamples.push_back("data_SingleElectron-Run2011A-PromptReco-v4");

    dataSamples.push_back("data_DoubleMu-Run2011A-05Aug2011-v1");
    dataSamples.push_back("data_SingleMu-Run2011A-05Aug2011-v1");
    dataSamples.push_back("data_DoubleElectron-Run2011A-05Aug2011-v1");
    dataSamples.push_back("data_SingleElectron-Run2011A-05Aug2011-v1");

    dataSamples.push_back("data_DoubleMu-Run2011A-PromptReco-v6");
    dataSamples.push_back("data_SingleMu-Run2011A-PromptReco-v6");
    dataSamples.push_back("data_DoubleElectron-Run2011A-PromptReco-v6");
    dataSamples.push_back("data_SingleElectron-Run2011A-PromptReco-v6");

    dataSamples.push_back("data_DoubleMu-Run2011B-PromptReco-v1");
    dataSamples.push_back("data_SingleMu-Run2011B-PromptReco-v1");
    dataSamples.push_back("data_DoubleElectron-Run2011B-PromptReco-v1");
    dataSamples.push_back("data_SingleElectron-Run2011B-PromptReco-v1");

    paramString += ":ZptCut:VptCut";
    treeName = "tVZCand";
    histName = "hVZMass_AllCuts";
    varName = "VZMass";
  }else{
    cerr<<"  Unknown input, aborting\n";
    abort();
  }
  
  TTree* tEvts = new TTree("tEvts", "Cut Values per sample");
  tEvts->ReadFile(config.c_str());
  
  double lumi = GetLumiUsed(f);
       
  tEvts->Draw(paramString.c_str(), "", "para goff");
  double n = tEvts->GetSelectedRows(); 
  if( debug_) cout<<"Found "<<n<<" samples "<<endl;
  for(int isample=0; isample<n; ++isample){
    const int SignalCode = tEvts->GetVal(0)[isample];
    const string SignalName = SampleName(SignalCode);
    const double mass = tEvts->GetVal(1)[isample];
    double minWindow = tEvts->GetVal(2)[isample];
    double maxWindow = tEvts->GetVal(3)[isample];

    string cuts;
    if(inName.find("WprimeWZ") != string::npos){ 
      const double minHt  = tEvts->GetVal(4)[isample];
      const double minZpt = tEvts->GetVal(5)[isample];
      const double minWpt = tEvts->GetVal(6)[isample];
      
      cuts = Form("(WZMass > %.0f && WZMass < %.0f && Ht > %.0f && Zpt > %.0f && Wpt > %.0f)*weight",
                  minWindow, maxWindow, minHt, minZpt, minWpt);
    }else if(inName.find("HadVZ") != string::npos){
      const double minZpt = tEvts->GetVal(4)[isample];
      const double minVpt = tEvts->GetVal(5)[isample];
      cuts = Form("(VZMass > %.0f && VZMass < %.0f && Zpt > %.0f && Vpt > %.0f)*weight",
                  minWindow, maxWindow, minZpt, minVpt);
    }
  
    if(debug_){
      cout<<"signalName: "<<SignalName<<endl;
      if(!useHists_)
        cout<<"Cuts are "<<cuts<<endl;
    }
      
    //Get Histograms
    vector<string> allSamples(BkgSamples);  
    allSamples.push_back(SignalName);
    TH1F *bkghist, *datahist, *allhist;
    if(useHists_){
      bkghist = get_sum_of_hists(f, BkgSamples, histName, 0, 1.);
      datahist = get_sum_of_hists(f, dataSamples, histName, 0, 1.);
      allhist = get_sum_of_hists(f, allSamples, histName, 0, 1.);
    }else{
      bkghist = new TH1F("bkghist", "bkghist", 250, 0, 2500);
      get_sum_of_hists(f, BkgSamples, treeName, varName, cuts, *bkghist);

      datahist = new TH1F("datahist", "datahist", 250, 0, 2500);
      get_sum_of_hists(f, dataSamples, treeName, varName, cuts, *datahist);
        
      allhist = new TH1F("allhist", "allhist", 250, 0, 2500);
      get_sum_of_hists(f, allSamples, treeName, varName, cuts, *allhist);
      allhist->SetLineColor(kRed);
    }        
    //Determine Mass Window
    int minBin, maxBin;
    if(windFrac > 0){
      double fitWindowLow = mass - 25;
      double fitWindowMax = mass + 25;
      allhist->Fit("gaus", "", "", fitWindowLow, fitWindowMax);
      TF1 *fit = allhist->GetFunction("gaus");
        
      double mean     = fit->GetParameter(1);
      double gaus_sig = fit->GetParameter(2);
      cout<<"Peak at "<<mean<<" +/- "<<gaus_sig<<endl;
        
      minWindow = mean - windFrac*gaus_sig;
      maxWindow = mean + windFrac*gaus_sig;
        
      minBin = allhist->FindBin(minWindow);
      maxBin = allhist->FindBin(maxWindow);   
    }else if(useHists_){
      minBin = allhist->FindBin(minWindow);
      maxBin = allhist->FindBin(maxWindow);        
    }else{
      minBin = 0;
      maxBin = allhist->GetNbinsX()+1;
    }

    if(debug_){
      cout<<"Integral from mass "<<minWindow<<" to "<<maxWindow<<" GeV."<<endl;
      cout<<"Integral from Bins "<<minBin<<" to "<<maxBin<<endl;  
    }

    double statMCEvts = 0.;
    double nMCEvts = allhist->IntegralAndError(minBin, maxBin, statMCEvts);

    double statBkgEvts = 0.;
    double nBkgEvts = bkghist->IntegralAndError(minBin, maxBin,statBkgEvts);

    double nSigEvts = nMCEvts - nBkgEvts;
    
    double DataEvts = datahist->Integral(minBin, maxBin);

    if(inName.find("WprimeWZ") != string::npos){ 
      const double ZJets = tEvts->GetVal(7)[isample];
      const double sZJets = tEvts->GetVal(8)[isample];

      //Not using anymore
      //AddDataDriven(nMCEvts, statMCEvts, ZJets, sZJets);
      //AddDataDriven(nBkgEvts, statBkgEvts, ZJets, sZJets);
    }
    if(debug_){
      cout<<"# of All Evts in Mass Window is "<<nMCEvts<<" +/- "<<statMCEvts<<" per "<<lumi<<" inv pb "<<endl;
      //cout<<"Total # of All Evts is "<<allhist->Integral()<<endl;
      cout<<"# of Bkg Evts in Mass Window is "<<nBkgEvts<<" +/- "<<statBkgEvts<<" per "<<lumi<<" inv pb "<<endl;
      //cout<<"Total # of Bkg Evts is "<<bkghist->Integral()<<endl;
      cout<<"# of Sig Evts in Mass Window is "<<nSigEvts<<" per "<<lumi<<" inv pb "<<endl;
      cout<<"# of Data Evts in Mass Window is "<<DataEvts<<" per "<<lumi<<" inv pb "<<endl;
      //cout<<"Total # of Data Evts is "<<datahist->Integral()<<endl;
    }

    //Read in Xsec from sample file
    double xsec = XSec(SignalName, mass);
    double nGenWeighted = lumi*xsec;
    double     Eff = nSigEvts / nGenWeighted;
    double statEff = TMath::Sqrt(Eff * (1-Eff)/nGenerated(SignalName)); 

    //Add in systematic errors for signal only (bkg now below)
    double sysMCEvts  = nSigEvts*SysErr(SignalName);
    double sysBkgEvts = 0.;
    double sysEff     = Eff*SysErr(SignalName);

    for(unsigned iBkg=0; iBkg<BkgSamples.size(); ++iBkg){
      vector<string> persample(1, BkgSamples[iBkg]);
      TH1F *perhist;
      if(useHists_){
        perhist = get_sum_of_hists(f, persample, histName, 0, 1.);
      }else{
        perhist = new TH1F("perhist", "perhist", 250, 0, 2500);
        get_sum_of_hists(f, persample, treeName, varName, cuts, *perhist);
      }
      double nPerHist = perhist->Integral(minBin, maxBin);

      if(debug_) cout<<BkgSamples[iBkg]<<": # of Evts in Mass Window is "<<nPerHist<<" per "<<lumi<<" inv pb "<<endl;
      double sysSample = nPerHist*SysErr(BkgSamples[iBkg]);
      double sysSampleWindow = nPerHist*BkgSysErrBySignal(SignalName);
      sysSample = AddInQuad(sysSampleWindow,  sysSample);
          
      sysMCEvts  = AddInQuad(sysMCEvts,  sysSample);
      sysBkgEvts = AddInQuad(sysBkgEvts, sysSample);

      delete perhist;
    }

    double sMCEvts    = AddInQuad(  statMCEvts,   sysMCEvts);
    double sBkgEvts   = AddInQuad( statBkgEvts,  sysBkgEvts);
    double sEff       = AddInQuad(statEff,      sysEff);

    out<<setprecision(0)
       <<SignalCode<<"\t"
       <<setw(6)<<mass<<"\t"
       <<setw(6)<<lumi<<"\t"
       <<setw(8)<<    DataEvts<<"\t"
       <<setprecision(4)
       <<setw(8)<<     nMCEvts<<"  "
       <<setw(8)<<  statMCEvts<<"  "
       <<setw(8)<<   sysMCEvts<<"  "
       <<setw(8)<<     sMCEvts<<"  "
       <<setw(8)<<    nBkgEvts<<"  "
       <<setw(8)<< statBkgEvts<<"  "
       <<setw(8)<<  sysBkgEvts<<"  "
       <<setw(8)<<    sBkgEvts<<"  "
       <<setw(6)<<    Eff<<"  "    
       <<setw(6)<<statEff<<"  "
       <<setw(6)<< sysEff<<"  "
       <<setw(6)<<   sEff<<"  "
       <<endl;
    c1->SaveAs((SignalName + ".pdf").c_str());

    //Clean Up
    delete bkghist;
    delete datahist;
    delete allhist;

  }//signal sample loop
  return;
}

void
AddDataDriven(double& Evts, double& sEvts, const double& DataDriven, const double&sDataDriven){
  Evts += DataDriven;
  sEvts = AddInQuad(sEvts, sDataDriven);
}
