//Usage: root -b -q 'ExpectedEvts.C+("inName.root", window)'

#include "consts.h"
void AddDataDriven(double& Evts, double& sEvts, const double& ZJets, const double&sZJets);

bool debug_ = false;
bool useHists_ = false;
bool noWind_ = false;
bool printTbl_ = false;
bool scaleMC_ = false;

void
ExpectedEvts(string inName, string config, int windFracTenths=-1, string opt=""){
  if(opt.find("debug") != string::npos) debug_ = true;
  if(opt.find("useHists") != string::npos) useHists_ = true;
  if(opt.find("noWind") != string::npos) noWind_ = true;
  if(opt.find("tbl") != string::npos) printTbl_ = true;
  if(opt.find("scaleMC") != string::npos) scaleMC_ = true;

  gErrorIgnoreLevel = kWarning;
  double windFrac = windFracTenths/10.;

  TFile *f = TFile::Open(inName.c_str(), "read");  assert(f);
  TCanvas* c1 = new TCanvas("c1", "Number of Events");
  string outfile("nEvents.txt");
  ofstream out(outfile.c_str());
  if(!out) { 
    cout << "Cannot open file " << outfile << endl; 
    abort();
  } 

  out<<"SignalCode/F:"
     <<"Mass/F:"
     <<"Lumi/F:"
     <<"DataEvts/F:"
     <<"MCEvts/F:"
     <<"statMCEvts/F:"
     <<"sysMCEvts/F:"
     <<"sMCEvts/F:"
     <<"SigEvts/F:"
     <<"statSigEvts/F:"
     <<"sysSigEvts/F:"
     <<"sSigEvts/F:"
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

  vector<string> BkgSamples, tblSamples; 
  vector<string> dataSamples;

  vector<float> BkgWeights; 

  string paramString = "SignalCode:Mass:nGen:bkgSysErr:minWindow:maxWindow";
  string treeName, histName, varName;
  if(inName.find("WprimeWZ") != string::npos){
    BkgSamples.push_back("TTJets"); 
    BkgSamples.push_back("ZZ"); 
    BkgSamples.push_back("WZJetsTo3LNu"); 
    BkgSamples.push_back("DYJetsToLL"); 
    BkgSamples.push_back("GVJets");  
    //BkgSamples.push_back("WWTo2L2Nu"); 
    //BkgSamples.push_back("WJetsToLNu"); 

    tblSamples.push_back("TTJets");
    //tblSamples.push_back("GVJets");
    tblSamples.push_back("ZZ");
    tblSamples.push_back("WZJetsTo3LNu");

    dataSamples.push_back("data");

    paramString += ":LtCut:ZptCut:WptCut:WZKFactor:ZJets:sZJets";
    treeName = "tEvts_MET";
    histName = "hWZMass_AllCuts";
    varName = "WZMass";
  }else if(inName.find("HadVZ") != string::npos){
    BkgSamples.push_back("Summer11_ZZ");
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
    dataSamples.push_back("data_DoubleElectron-Run2011B-PromptReco-v1");

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
       
  if(printTbl_){
    cout<<" Mass ";
    if(!noWind_) cout<<"& Window";
    for(unsigned iBkg=0; iBkg<BkgSamples.size(); ++iBkg){
      for(unsigned iTbl=0; iTbl<tblSamples.size(); ++iTbl){
        if(BkgSamples[iBkg] == tblSamples[iTbl])
          cout<<" & "<<BkgSamples[iBkg];
      }
    }
    cout<<" & Total & Data & $N_{Sig}$ \\\\"<<endl;
  }

  tEvts->Draw(paramString.c_str(), "", "para goff");
  double n = tEvts->GetSelectedRows(); 
  if( debug_) cout<<"Found "<<n<<" samples "<<endl;
  for(int isample=0; isample<n; ++isample){
    int treeIdx = 0;
    const double SignalCode = tEvts->GetVal(treeIdx++)[isample];
    const vector<string> SignalNames = SampleName(SignalCode);
    const double mass = tEvts->GetVal(treeIdx++)[isample];
    const double nGen = tEvts->GetVal(treeIdx++)[isample];
    const double bkgSysErr = tEvts->GetVal(treeIdx++)[isample];
    double minWindow = tEvts->GetVal(treeIdx++)[isample];
    double maxWindow = tEvts->GetVal(treeIdx++)[isample];
    if(noWind_){
      minWindow = -1;
      maxWindow = 9e9;
    }

    string cuts;
    if(inName.find("WprimeWZ") != string::npos){ 
      const double minLt  = tEvts->GetVal(treeIdx++)[isample];
      const double minZpt = tEvts->GetVal(treeIdx++)[isample];
      const double minWpt = tEvts->GetVal(treeIdx++)[isample];
      cuts = Form("(WZMass > %.0f && WZMass < %.0f && Lt > %.0f && Zpt > %.0f && Wpt > %.0f)*weight",
                  minWindow, maxWindow, minLt, minZpt, minWpt);

      if(scaleMC_){
        float wZee = 0.960;
        float wZmm = 1.050;
        float wWen = 1.046;
        float wWmn = 1.021;
        float w3e = wZee * wWen;
        float w2e = wZee * wWmn;
        float w1e = wZmm * wWen;
        float w0e = wZmm * wWmn;
        cuts += Form("*( (weight==1.0) + (weight!=1.0)*((EvtType==0)*%f + (EvtType==1)*%f + (EvtType==2)*%f + (EvtType==3)*%f))", w0e,w1e,w2e,w3e);
      }

      BkgWeights.assign(BkgSamples.size(), 1.);
      const double wzKFac = tEvts->GetVal(treeIdx++)[isample];
      for(unsigned iBkg=0; iBkg<BkgSamples.size(); ++iBkg){
        if(BkgSamples[iBkg].find("WZ") != string::npos){
          BkgWeights[iBkg] = wzKFac;
        }
      }

    }else if(inName.find("HadVZ") != string::npos){
      const double minZpt = tEvts->GetVal(treeIdx++)[isample];
      const double minVpt = tEvts->GetVal(treeIdx++)[isample];
      cuts = Form("(VZMass > %.0f && VZMass < %.0f && Zpt > %.0f && Vpt > %.0f)*weight",
                  minWindow, maxWindow, minZpt, minVpt);
    }
  
    if(debug_){
      for(unsigned i=0; i<SignalNames.size(); ++i) cout<<"signalName: "<<SignalNames[i]<<endl;
      if(!useHists_)
        cout<<"Cuts are "<<cuts<<endl;
    }
      
    //Get Histograms
    vector<string> allSamples(BkgSamples);  
    vector<float> allWeights(BkgWeights);  
    
    for(unsigned i=0; i<SignalNames.size(); ++i){
      allSamples.push_back(SignalNames[i]);
      if(allWeights.size()) allWeights.push_back(1.0);
    }

    TH1F *bkghist, *datahist, *allhist;
    if(useHists_){
      bkghist = get_sum_of_hists(f, BkgSamples, histName, 0, 1.);
      datahist = get_sum_of_hists(f, dataSamples, histName, 0, 1.);
      allhist = get_sum_of_hists(f, allSamples, histName, 0, 1.);
    }else{
      bkghist = new TH1F("bkghist", "bkghist", 250, 0, 2500);
      get_sum_of_hists(f, BkgSamples, treeName, varName, cuts, *bkghist, BkgWeights);

      datahist = new TH1F("datahist", "datahist", 250, 0, 2500);
      get_sum_of_hists(f, dataSamples, treeName, varName, cuts, *datahist);
        
      allhist = new TH1F("allhist", "allhist", 250, 0, 2500);
      get_sum_of_hists(f, allSamples, treeName, varName, cuts, *allhist, allWeights);
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

    double statSigEvts = sqrt(statMCEvts*statMCEvts - statBkgEvts*statBkgEvts);
    double nSigEvts = nMCEvts - nBkgEvts;

    double DataEvts = datahist->Integral(minBin, maxBin);

    if(inName.find("WprimeWZ") != string::npos){ 
      const double  ZJets = tEvts->GetVal(treeIdx++)[isample];
      const double sZJets = tEvts->GetVal(treeIdx++)[isample];

      if(0){      //Not using anymore
        AddDataDriven(nMCEvts, statMCEvts, ZJets, sZJets);
        AddDataDriven(nBkgEvts, statBkgEvts, ZJets, sZJets);
      }
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
    //double xsec = GetFileInfo(hInfo, "#sigma (pb)");
    //double lumi = GetFileInfo(hInfo, "#intL dt (pb^{-1})")
    //double nGen = GetFileInfo(hInfo, "Number of Events Produced");

    double xsec = XSec(SignalNames[0], mass);//Cory: this is wrong, how should i combine xsecs???
    double nGenWeighted = lumi*xsec;
    double     Eff = nSigEvts / nGenWeighted;
    double statEff = TMath::Sqrt(Eff * (1-Eff)/nGen); 
    //cout<<" mass: "<<mass<<" nSig:"<<nSigEvts<<" nGenWeighted:"<<nGenWeighted<<" eff:"<<Eff<<" TESTnSigEvts:"<<TESTnSigEvts<<endl;
    //Add in systematic errors for signal only (bkg now below)
    double sysMCEvts  = nSigEvts*SysErr(SignalNames[0]);
    double sysSigEvts = nSigEvts*SysErr(SignalNames[0]);
    double sysBkgEvts = 0.;
    double sysEff     = Eff*SysErr(SignalNames[0]);

    if(printTbl_){
      cout<<"W' "<<setiosflags(ios::fixed)<<setprecision(0)<<mass;
      if(!noWind_) cout<<" & "<<minWindow<<"-"<<maxWindow;
    }

    for(unsigned iBkg=0; iBkg<BkgSamples.size(); ++iBkg){
      vector<string> persample(1, BkgSamples[iBkg]);
      TH1F *perhist;
      if(useHists_){
        perhist = get_sum_of_hists(f, persample, histName, 0, 1.);
      }else{
        perhist = new TH1F("perhist", "perhist", 250, 0, 2500);
        get_sum_of_hists(f, persample, treeName, varName, cuts, *perhist);
      }
      double sPerHist = 0;
      double nPerHist = perhist->IntegralAndError(minBin, maxBin, sPerHist);
      if(printTbl_){
        for(unsigned iTbl=0; iTbl<tblSamples.size(); ++iTbl){
          if(BkgSamples[iBkg] == tblSamples[iTbl])
            cout<<" & "<<Value(nPerHist, sPerHist)<<" $\\pm$ "<<Value(sPerHist);
          //cout<<" & "<<std::fixed << std::setprecision(1)<<nPerHist<<" $\\pm$ "<<sPerHist;
        }
      }

      if(debug_) cout<<BkgSamples[iBkg]<<": # of Evts in Mass Window is "<<nPerHist<<" per "<<lumi<<" inv pb "<<endl;
      double sysSample = nPerHist*SysErr(BkgSamples[iBkg]); //Bkg due to this bkg
      double sysSampleWindow = nPerHist*bkgSysErr; //Bkg error in this signal window
      sysSample = AddInQuad(sysSampleWindow,  sysSample);
          
      sysMCEvts  = AddInQuad(sysMCEvts,  sysSample);
      sysBkgEvts = AddInQuad(sysBkgEvts, sysSample);

      delete perhist;
    }

    double sMCEvts    = AddInQuad(  statMCEvts,   sysMCEvts);
    double sSigEvts   = AddInQuad( statSigEvts,  sysSigEvts);
    double sBkgEvts   = AddInQuad( statBkgEvts,  sysBkgEvts);
    double sEff       = AddInQuad(statEff,      sysEff);

    if(printTbl_){ 
      cout<<" & "<<Value(nBkgEvts,statBkgEvts)<<" $\\pm$ "<<Value(statBkgEvts);
      //cout<<" & "<<std::fixed << std::setprecision(1)<<nBkgEvts<<" $\\pm$ "<<statBkgEvts;
      cout<<" & "<<(int)DataEvts;
      cout<<" & "<<Value(nSigEvts,statSigEvts)<<" $\\pm$ "<<Value(statSigEvts);
      cout<<" \\\\ \\hline"<<endl;
    }

    out<<setprecision(1)
       <<SignalCode<<"\t"
       <<setprecision(0)
       <<setw(6)<<mass<<"\t"
       <<setw(6)<<lumi<<"\t"
       <<setw(8)<<    DataEvts<<"\t"
       <<setprecision(4)
       <<setw(8)<<     nMCEvts<<"  "
       <<setw(8)<<  statMCEvts<<"  "
       <<setw(8)<<   sysMCEvts<<"  "
       <<setw(8)<<     sMCEvts<<"  "
       <<setw(8)<<    nSigEvts<<"  "
       <<setw(8)<< statSigEvts<<"  "
       <<setw(8)<<  sysSigEvts<<"  "
       <<setw(8)<<    sSigEvts<<"  "
       <<setw(8)<<    nBkgEvts<<"  "
       <<setw(8)<< statBkgEvts<<"  "
       <<setw(8)<<  sysBkgEvts<<"  "
       <<setw(8)<<    sBkgEvts<<"  "
       <<setw(6)<<    Eff<<"  "    
       <<setw(6)<<statEff<<"  "
       <<setw(6)<< sysEff<<"  "
       <<setw(6)<<   sEff<<"  "
       <<endl;
    if(0) c1->SaveAs((SignalNames[0] + ".pdf").c_str());

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
