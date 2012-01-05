//Usage: root -b -q 'printNEvts.C+(file, evtType)'

#include <vector>
#include "THStack.h"
#include "TStyle.h"
#include "TLine.h"
#include "../Limits/consts.h"

void
printNEvts(string infile, int evtType=-1){  
  TFile *f = TFile::Open(infile.c_str(), "read"); assert(f);

  vector<pair<vector<string>, string> > Samples;

  if(infile.find("WprimeWZ") != string::npos){
    Samples.push_back(make_pair(vector<string>(1,"DYJetsToLL"),"$Z+jets$"));
    Samples.push_back(make_pair(vector<string>(1,"TTJets"), "$t\\bar{t}$"));
    Samples.push_back(make_pair(vector<string>(1,"ZZ"), "$ZZ$"));
    Samples.push_back(make_pair(vector<string>(1,"GVJets"), "$V\\gamma$"));
    Samples.push_back(make_pair(vector<string>(1,"WJetsToLNu"), "$W+jets$"));
    Samples.push_back(make_pair(vector<string>(1,"WWTo2L2Nu"), "$WW$"));
    Samples.push_back(make_pair(vector<string>(1,"WZJetsTo3LNu"), "$WZ$"));
  }else if(infile.find("HadVZ") != string::npos){
    Samples.push_back(make_pair(vector<string>(1,"Summer11_DYJetsToLL_PtZ100"),"$Z+jets$"));
    Samples.push_back(make_pair(vector<string>(1,"Summer11_TTJets"), "$t\\bar{t}$"));
    Samples.push_back(make_pair(vector<string>(1,"Summer11_ZZ"), "$ZZ$"));
    Samples.push_back(make_pair(vector<string>(1,"Summer11_VGamma"), "$V\\gamma$"));
    Samples.push_back(make_pair(vector<string>(1,"Summer11_WJets_PtW100"), "$W+jets$"));
    Samples.push_back(make_pair(vector<string>(1,"Summer11_WW"), "$WW$"));
    Samples.push_back(make_pair(vector<string>(1,"Summer11_WZ"), "$WZ$"));
  }
  cout << std::fixed << std::setprecision(1);
    
  //total bkg
  vector<string> bkg;
  for(unsigned i=0; i<Samples.size(); ++i){
    for(unsigned j=0; j<Samples[i].first.size(); ++j){
      bkg.push_back(Samples[i].first[j]);
    }
  }
  Samples.push_back(make_pair(vector<string>(), "BREAK"));//Print hline
  Samples.push_back(make_pair(bkg, "Total Background"));
  Samples.push_back(make_pair(vector<string>(), "BREAK"));//Print hline

  vector<string> Data;
  if(infile.find("WprimeWZ") != string::npos || infile.find("EWKWZ") != string::npos){
    Data.push_back("data");
  }else if(infile.find("HadVZ") != string::npos){
    Data.push_back("data_DoubleMu-Run2011A-May10ReReco-v1");
    Data.push_back("data_SingleMu-Run2011A-May10ReReco-v1");
    Data.push_back("data_DoubleElectron-Run2011A-May10ReReco-v1");

    Data.push_back("data_DoubleMu-Run2011A-PromptReco-v4");
    Data.push_back("data_SingleMu-Run2011A-PromptReco-v4");
    Data.push_back("data_DoubleElectron-Run2011A-PromptReco-v4");

    Data.push_back("data_DoubleMu-Run2011A-05Aug2011-v1");
    Data.push_back("data_SingleMu-Run2011A-05Aug2011-v1");
    Data.push_back("data_DoubleElectron-Run2011A-05Aug2011-v1");

    Data.push_back("data_DoubleMu-Run2011A-PromptReco-v6");
    Data.push_back("data_SingleMu-Run2011A-PromptReco-v6");
    Data.push_back("data_DoubleElectron-Run2011A-PromptReco-v6");
    
    Data.push_back("data_DoubleMu-Run2011B-PromptReco-v1");
    Data.push_back("data_DoubleElectron-Run2011B-PromptReco-v1");
  } 
  Samples.push_back(make_pair(Data, "Data"));
  Samples.push_back(make_pair(vector<string>(), "BREAK"));//Print hline

  //signal
  if(infile.find("WprimeWZ") != string::npos || infile.find("EWKWZ") != string::npos){
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-200"), "200"));
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-250"), "250"));
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-300"), "300"));
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-400"), "400"));
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-500"), "500"));
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-600"), "600"));
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-700"), "700"));
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-800"), "800"));
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-900"), "900"));
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-1000"), "1000"));
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-1100"), "1100"));
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-1200"), "1200"));
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-1300"), "1300"));
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-1400"), "1400"));
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-1500"), "1500"));
  }else if(infile.find("HadVZ") != string::npos){
    vector<string> RS750;
    RS750.push_back("Summer11_RSZZeejj_750");
    RS750.push_back("Summer11_RSZZmmjj_750");
    Samples.push_back(make_pair(RS750, "$G_{RS}$ (750 \\GeV)"));

    vector<string> RS1000;
    RS1000.push_back("Summer11_RSZZeejj_1000");
    RS1000.push_back("Summer11_RSZZmmjj_1000");
    Samples.push_back(make_pair(RS1000, "$G_{RS}$ (1000 \\GeV)"));

   vector<string> RS1250;
    RS1250.push_back("Summer11_RSZZeejj_1250");
    RS1250.push_back("Summer11_RSZZmmjj_1250");
    Samples.push_back(make_pair(RS1250, "$G_{RS}$ (1250 \\GeV)"));

    vector<string> RS1500;
    RS1500.push_back("Summer11_RSZZeejj_1500");
    RS1500.push_back("Summer11_RSZZmmjj_1500");
    Samples.push_back(make_pair(RS1500, "$G_{RS}$ (1500 \\GeV)"));

    vector<string> RS1750;
    RS1750.push_back("Summer11_RSZZeejj_1750");
    RS1750.push_back("Summer11_RSZZmmjj_1750");
    Samples.push_back(make_pair(RS1750, "$G_{RS}$ (1750 \\GeV)"));

    vector<string> RS2000;
    RS2000.push_back("Summer11_RSZZeejj_2000");
    RS2000.push_back("Summer11_RSZZmmjj_2000");
    Samples.push_back(make_pair(RS2000, "$G_{RS}$ (2000 \\GeV)"));
    
    Samples.push_back(make_pair(vector<string>(), "BREAK"));//Print hline
    Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_300"), "\\Wprime (300 \\GeV)"));
    Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_400"), "\\Wprime (400 \\GeV)"));
    //Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_500"), "\\Wprime (500 \\GeV)"));
    Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_600"), "\\Wprime (600 \\GeV)"));
    Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_700"), "\\Wprime (700 \\GeV)"));
    Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_800"), "\\Wprime (800 \\GeV)"));
    Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_900"), "\\Wprime (900 \\GeV)"));
    Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_1000"), "\\Wprime (1000 \\GeV)"));
    //Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_1100"), "\\Wprime (1100 \\GeV)"));
    Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_1200"), "\\Wprime (1200 \\GeV)"));
    //Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_1300"), "\\Wprime (1300 \\GeV)"));
    Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_1400"), "\\Wprime (1400 \\GeV)"));
    Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_1500"), "\\Wprime (1500 \\GeV)"));
  }

  //////////////////////////////

  vector<pair<string, string> > levels;
  if(infile.find("WprimeWZ") != string::npos || infile.find("EWKWZ") != string::npos){
    if(evtType == -1){
      levels.push_back(make_pair("MinNLeptons", "Preselection"));
      levels.push_back(make_pair("ValidZ", "Z Selection"));
    }
    levels.push_back(make_pair("ValidW", "W Selection"));
    levels.push_back(make_pair("MET", "\\MET"));
    levels.push_back(make_pair("ValidWZCand", "WZ"));
  }else if(infile.find("HadVZ") != string::npos){
    if(evtType == -1){
    }
    levels.push_back(make_pair("Zpt", "Z Selection"));
    levels.push_back(make_pair("VMass", "V Selection"));
  }  

  cout<<" Sample ";
  for(unsigned level=0; level<levels.size(); ++level) cout<<" & "<<levels[level].second;
  cout<<" \\\\ \\hline"<<endl;

  for(unsigned i=0; i<Samples.size(); ++i){
    if(Samples[i].second == "BREAK"){
      cout<<" \\hline"<<endl;
      continue;
    }

    cout<<Samples[i].second;

    for(unsigned level=0; level<levels.size(); ++level){
      string hName;
      int bin = -1;
      if(evtType == -1){
        hName = "hNumEvts";
      }else{
        hName = "hEvtType_" + levels[level].first;
      }
      TH1F* hist = get_sum_of_hists(f, Samples[i].first, hName.c_str());
      if(evtType == -1){
        bin = hist->GetXaxis()->FindBin(levels[level].first.c_str());
      }else{
        bin = hist->FindBin(evtType);
      }
      float tot = hist->GetBinContent(bin);
      float sigma = hist->GetBinError(bin);
      
      cout<<" & "<<Value(tot,sigma)<<" $\\pm$ "<<Value(sigma);
    }//loop over cuts
    cout<<" \\\\ \\hline"<<endl;
  }//loop over samples
}
