//Author: Cory Fantasia 2012
//Purpose: Print latex tables of number of events and errors
//Usage: root -b -l -q 'printNEvts.C+(file, mode)'
//eg: root -b -l -q 'printNEvts.C+("../../../EWKWZ.root", -1)'

#include <vector>
#include "THStack.h"
#include "TStyle.h"
#include "TLine.h"
#include "../Limits/consts.h"

struct Channel{
  string name;
  vector<int> sub;//subchannels
  Channel(){}
  Channel(string n, int s1=-9999, int s2=-999, int s3=-999, int s4=-999){
    name = n; 
    if(s1 != -999) sub.push_back(s1);
    if(s2 != -999) sub.push_back(s2);
    if(s3 != -999) sub.push_back(s3);
    if(s4 != -999) sub.push_back(s4);
  }
//  Channel(string n, int[] arr){
//    name = n; 
//    sub = vector<int> (arr, arr + sizeof(arr) / sizeof(arr[0]) );
//  }

};

void
printNEvts(string infile, int mode=-1){  
  TFile *f = TFile::Open(infile.c_str(), "read"); assert(f);

  vector<pair<vector<string>, string> > Samples;

  if(infile.find("WprimeWZ") != string::npos  || infile.find("EWKWZ") != string::npos){
    Samples.push_back(make_pair(vector<string>(1,"DYJetsToLL"),"$Z+jets$"));
    Samples.push_back(make_pair(vector<string>(1,"TTJets"), "$t\\bar{t}$"));
    Samples.push_back(make_pair(vector<string>(1,"ZZ"), "$ZZ$"));
    Samples.push_back(make_pair(vector<string>(1,"GVJets"), "$V\\gamma$"));
    //Samples.push_back(make_pair(vector<string>(1,"WJetsToLNu"), "$W+jets$"));
    //Samples.push_back(make_pair(vector<string>(1,"WWTo2L2Nu"), "$WW$"));
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
  Samples.push_back(make_pair(bkg, infile.find("EWKWZ") != string::npos ? "Total MC" : "Total Background"));
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

  //signal
  if(infile.find("WprimeWZ") != string::npos){
    Samples.push_back(make_pair(vector<string>(), "BREAK"));//Print hline
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-200"), "200"));
    //Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-250"), "250"));
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
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-1600"), "1600"));
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-1700"), "1700"));
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-1800"), "1800"));
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-1900"), "1900"));
    Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-2000"), "2000"));
  }else if(infile.find("HadVZ") != string::npos){
    Samples.push_back(make_pair(vector<string>(), "BREAK"));//Print hline
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

  vector<Channel> evtTypes;
  vector<pair<string, string> > levels;
  if(infile.find("WprimeWZ") != string::npos || infile.find("EWKWZ") != string::npos){
    if(mode == -1){
      evtTypes.push_back(Channel("All Channels", 0,1,2,3));
      levels.push_back(make_pair("MinNLeptons", "Preselection"));
      levels.push_back(make_pair("ValidZ", "Z Selection"));
      levels.push_back(make_pair("ValidW", "W Selection"));
      levels.push_back(make_pair("MET", "\\MET"));
    }else if(mode == -2){
      evtTypes.push_back(Channel("ee", 0,1));
      evtTypes.push_back(Channel("$\\mu\\mu$", 2,3));
      levels.push_back(make_pair("ValidZ", "Z Selection"));
    }else if(mode == -4){
      evtTypes.push_back(Channel("3e", 0));
      evtTypes.push_back(Channel("2e1$\\mu$", 1));
      evtTypes.push_back(Channel("1e2$\\mu$", 2));
      evtTypes.push_back(Channel("3$\\mu$", 3));
      levels.push_back(make_pair("ValidW", "W Selection"));
    }else{
      evtTypes.push_back(Channel("3e", 0));
      evtTypes.push_back(Channel("2e1$\\mu$", 1));
      evtTypes.push_back(Channel("1e2$\\mu$", 2));
      evtTypes.push_back(Channel("3$\\mu$", 3));
      levels.push_back(make_pair("MET", "\\MET"));
    }
  }else if(infile.find("HadVZ") != string::npos){
    if(mode == -1){
      evtTypes.push_back(Channel("All Channels", 0,2));
    }else{
      evtTypes.push_back(Channel("ee", 0));
      evtTypes.push_back(Channel("$\\mu$$\\mu$", 2));
    }
    levels.push_back(make_pair("Zpt", "Z Selection"));
    levels.push_back(make_pair("VMass", "V Selection"));
  }  

  //%Cory: Updated DATE
  cout<<"\\begin{table}[!h] \\centering \\begin{tabular}{|c||*{"<<levels.size()*evtTypes.size()<<"}{c|}} \\hline"<<endl;

  if(mode == -1){
    cout<<"Sample ";
    for(unsigned level=0; level<levels.size(); ++level){
      cout<<" & "<<levels[level].second;//Print cut names
    }
  }else{//don't do if doing all channels
    cout<<"\\multirow{2}{*}{Sample} ";
    for(unsigned level=0; level<levels.size(); ++level){
      cout<<" & \\multicolumn{"<<evtTypes.size()<<"}{c|}{$"<<levels[level].second<<"$}";//Print Cut Names
    }
    cout<<" \\\\ "<<endl;
    cout<<" \\cline{2-"<<levels.size()*evtTypes.size()+1<<"}"<<endl;
    for(unsigned level=0; level<levels.size(); ++level){
      for(unsigned channel=0; channel<evtTypes.size(); ++channel){//multi col
        cout<<" & "<<evtTypes[channel].name;//Print cut names
      }
    }
  }
  cout<<" \\\\ \\hline"<<endl;

  for(unsigned i=0; i<Samples.size(); ++i){
    //if(Samples[i].second.find("Data") == string::npos) continue;//Cory: test only data

    if(Samples[i].second == "BREAK"){
      cout<<" \\hline"<<endl;
      continue;
    }

    cout<<Samples[i].second;//Print Sample names
    
    for(unsigned level=0; level<levels.size(); ++level){
      for(unsigned channel=0; channel<evtTypes.size(); ++channel){
        float tot(0), sigma(0);
        /*
        //Use trees
        string tName = "tEvts_" + levels[level].first;
        string cuts = "(weight)*(ZMass > 71.188)*(ZMass < 111.188)";
        int nbins = evtTypes[channel].sub.size();
        for(int ibin=0; ibin<nbins; ++ibin){
          cuts += Form("*(EvtType == %i)", evtTypes[channel].sub[ibin]);
        }
        //cout<<"Cuts are "<<cuts<<endl;
        TTree* tree = getTree(f, Samples[i].first, tName); assert (tree || !(std::cerr << "Failed getting : " << tName << endl));
        Value val = GetNEvtsAndError(tree, cuts);
        float tot = val.val;
        float sigma = val.err;
        */

        if(mode == -1 && 0){
          //Use EvtType Histos
          string hName = "hEvtType_" + levels[level].first;//Cory: Should be modified to also look at + and - Ws only
          //hName = "hEvtTypeP_" + levels[level].first;//Cory: Should be modified to also look at + and - Ws only
          //hName = "hEvtTypeM_" + levels[level].first;//Cory: Should be modified to also look at + and - Ws only

          TH1F* hist = get_sum_of_hists(f, Samples[i].first, hName);
        
          int nbins = evtTypes[channel].sub.size();
          for(int ibin=0; ibin<nbins; ++ibin){
            int bin = hist->FindBin(evtTypes[channel].sub[ibin]);
            if(bin == hist->GetNbinsX()+1) cout<<"\n\nUsing an overflow bin\n\n";
            tot += hist->GetBinContent(bin);///////Cory: is this right?, += or =?
            sigma = AddInQuad(hist->GetBinError(bin), sigma);
          }//binOfffset
        }else{
          //Use Zmass by channel histos
          int nbins = evtTypes[channel].sub.size();
          for(int ibin=0; ibin<nbins; ++ibin){
            int & a = evtTypes[channel].sub[ibin];
            string hName;
            if(levels[level].second == "Preselection"){
              if(mode == -1 && ibin > 0) break;//Don't double count
              hName = "hMET_" + levels[level].first;//Cory: Should be modified to also look at + and - Ws only
            }else if(levels[level].second == "Z Selection"){
              if(mode == -1 && ibin > 0) break;//Don't double count
              if(mode == -2 && ibin%2==1 ) break;//Don't double count
              if(mode == -1){
                hName = "hZMass_" + levels[level].first;
              }else{
                if(a < 1.5) hName = "hZeeMass_" + levels[level].first;
                else        hName = "hZmmMass_" + levels[level].first;
              }
            }else{
              hName = Form("hZ%ie%imMass_",3-a,a) + levels[level].first;
            }
            TH1F* hist = get_sum_of_hists(f, Samples[i].first, hName);
            //cout<<"Hname is "<<hName<<endl;
            int bin1 = 0;//hist->FindBin(0);
            int bin2 = hist->GetNbinsX()+1;
            double t(0),s(0);
            t = hist->IntegralAndError(bin1, bin2, s);

            tot += t;
            sigma = AddInQuad(s, sigma);
          }//binOfffset
        }
        
        cout<<" & "<<Value(tot,sigma);
        if(Samples[i].second.find("Data") == string::npos) cout<<" $\\pm$ "<<Value(sigma);//Only print error for MC
      }//loop over evtTypes
    }//loop over cuts
    cout<<" \\\\ \\hline"<<endl;
  }//loop over samples
  cout<<"\\end{tabular} \\end{table}"<<endl;
}
