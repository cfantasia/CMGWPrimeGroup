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
    Data.push_back("data_SingleElectron-Run2011A-May10ReReco-v1");

    Data.push_back("data_DoubleMu-Run2011A-PromptReco-v4");
    Data.push_back("data_SingleMu-Run2011A-PromptReco-v4");
    Data.push_back("data_DoubleElectron-Run2011A-PromptReco-v4");
    Data.push_back("data_SingleElectron-Run2011A-PromptReco-v4");

    Data.push_back("data_DoubleMu-Run2011A-05Aug2011-v1");
    Data.push_back("data_SingleMu-Run2011A-05Aug2011-v1");
    Data.push_back("data_DoubleElectron-Run2011A-05Aug2011-v1");
    Data.push_back("data_SingleElectron-Run2011A-05Aug2011-v1");

    Data.push_back("data_DoubleMu-Run2011A-PromptReco-v6");
    Data.push_back("data_SingleMu-Run2011A-PromptReco-v6");
    Data.push_back("data_DoubleElectron-Run2011A-PromptReco-v6");
    Data.push_back("data_SingleElectron-Run2011A-PromptReco-v6");
    
    Data.push_back("data_DoubleMu-Run2011B-PromptReco-v1");
    Data.push_back("data_DoubleElectron-Run2011B-PromptReco-v1");
  } 
  Samples.push_back(make_pair(Data, "Data"));
  Samples.push_back(make_pair(vector<string>(), "BREAK"));//Print hline

  //signal
  if(infile.find("WprimeWZ") != string::npos || infile.find("EWKWZ") != string::npos){
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
    Samples.push_back(make_pair(RS750, "$G_{RS}$ 750"));

    vector<string> RS1000;
    RS1000.push_back("Summer11_RSZZeejj_1000");
    RS1000.push_back("Summer11_RSZZmmjj_1000");
    Samples.push_back(make_pair(RS1000, "$G_{RS}$ 1000"));

   vector<string> RS1250;
    RS1250.push_back("Summer11_RSZZeejj_1250");
    RS1250.push_back("Summer11_RSZZmmjj_1250");
    Samples.push_back(make_pair(RS1250, "$G_{RS}$ 1250"));

    vector<string> RS1500;
    RS1500.push_back("Summer11_RSZZeejj_1500");
    RS1500.push_back("Summer11_RSZZmmjj_1500");
    Samples.push_back(make_pair(RS1500, "$G_{RS}$ 1500"));

    vector<string> RS1750;
    RS1750.push_back("Summer11_RSZZeejj_1750");
    RS1750.push_back("Summer11_RSZZmmjj_1750");
    Samples.push_back(make_pair(RS1750, "$G_{RS}$ 1750"));

    vector<string> RS2000;
    RS2000.push_back("Summer11_RSZZeejj_2000");
    RS2000.push_back("Summer11_RSZZmmjj_2000");
    Samples.push_back(make_pair(RS2000, "$G_{RS}$ 2000"));

    Samples.push_back(make_pair(vector<string>(), "BREAK"));//Print hline

    Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_300"), "W' 300"));
    Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_400"), "W' 400"));
    //Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_500"), "W' 500"));
    Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_600"), "W' 600"));
    Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_700"), "W' 700"));
    Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_800"), "W' 800"));
    Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_900"), "W' 900"));
    Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_1000"), "W' 1000"));
    //Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_1100"), "W' 1100"));
    Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_1200"), "W' 1200"));
    //Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_1300"), "W' 1300"));
    Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_1400"), "W' 1400"));
    Samples.push_back(make_pair(vector<string>(1, "Summer11_WPrimeZZlljj_1500"), "W' 1500"));
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
      levels.push_back(make_pair("MinNJets", "Preselection"));
    }
    levels.push_back(make_pair("ZMass", "Z Selection"));
    levels.push_back(make_pair("VMass", "V Selection"));
    //levels.push_back(make_pair("ValidVZ", "VZ"));
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
      float tot = 0;
      float sigma2 = 0;
      for(unsigned subsam=0; subsam<Samples[i].first.size(); ++subsam){
        TH1F* hist = NULL;
        int bin = -1;
        if(evtType == -1){
          string hist_name = Samples[i].first[subsam] + "/hNumEvts";
          hist = (TH1F*) f->Get(hist_name.c_str()); 
          if(!hist){
            cout<<"\n\nFailed getting "<<hist_name<<endl;
            abort();
          }
          bin = hist->GetXaxis()->FindBin(levels[level].first.c_str());
        }else{
          string hist_name = Samples[i].first[subsam] + "/hEvtType_" + levels[level].first;
          hist = (TH1F*) f->Get(hist_name.c_str()); 
          if(!hist){
            cout<<"\n\nDidn't find histo "<<hist_name<<endl;
            abort();
          }
          bin = hist->FindBin(evtType);
        }
        tot += hist->GetBinContent(bin);//underflow
        sigma2 += hist->GetBinError(bin)*hist->GetBinError(bin);

      }//loop over subsamples
      cout<<" & "<<tot<<" $\\pm$ "<<sqrt(sigma2);
    }//loop over cuts
    cout<<" \\\\ \\hline"<<endl;
  }//loop over samples
}
