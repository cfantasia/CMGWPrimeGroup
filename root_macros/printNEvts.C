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
  Samples.push_back(make_pair(vector<string>(1,"DYJetsToLL"),"$Z+jets$"));
  Samples.push_back(make_pair(vector<string>(1,"TTJets"), "$t\\bar{t}$"));
  Samples.push_back(make_pair(vector<string>(1,"ZZ"), "$ZZ$"));
  Samples.push_back(make_pair(vector<string>(1,"GVJets"), "$V\\gamma$"));
  Samples.push_back(make_pair(vector<string>(1, "WJetsToLNu"), "$W+jets$"));
  Samples.push_back(make_pair(vector<string>(1,"WWTo2L2Nu"), "$WW$"));
  Samples.push_back(make_pair(vector<string>(1,"WZJetsTo3LNu"), "$WZ$"));

  cout << std::fixed << std::setprecision(1);
    
  //total bkg
  vector<string> bkg;
  for(unsigned i=0; i<Samples.size(); ++i){
    for(unsigned j=0; j<Samples[i].first.size(); ++j){
      bkg.push_back(Samples[i].first[j]);
    }
  }
  Samples.push_back(make_pair(bkg, "Total Background"));

  Samples.push_back(make_pair(vector<string>(1, "data"), "Data"));

  //signal
  Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-300"), "300"));
  Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-400"), "400"));
  //Samples.push_back(make_pair(vector<string>(1, "WprimeToWZTo3LNu_M-500"), "500"));
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
  
  //////////////////////////////

  vector<pair<string, string> > levels;
  if(evtType == -1){
    levels.push_back(make_pair("MinNLeptons", "Preselection"));
    levels.push_back(make_pair("ValidZ", "Z Selection"));
  }
  levels.push_back(make_pair("ValidW", "W Selection"));
  levels.push_back(make_pair("MET", "\\MET"));
  levels.push_back(make_pair("ValidWZCand", "WZ"));

  cout<<" Sample ";
  for(unsigned level=0; level<levels.size(); ++level) cout<<" & "<<levels[level].second;
  cout<<" \\\\ \\hline"<<endl;

  for(unsigned i=0; i<Samples.size(); ++i){
    cout<<Samples[i].second;

    for(unsigned level=0; level<levels.size(); ++level){
      float tot = 0;
      for(unsigned subsam=0; subsam<Samples[i].first.size(); ++subsam){
        if(evtType == -1){
          string hist_name = Samples[i].first[subsam] + "/hNumEvts";
          TH1F* hist = (TH1F*) f->Get(hist_name.c_str());
          int bin = hist->GetXaxis()->FindBin(levels[level].first.c_str());
          tot += hist->GetBinContent(bin);
          //tot += hist->Integral();
        }else{
          string hist_name = Samples[i].first[subsam] + "/hEvtType_" + levels[level].first;
          TH1F* hist = (TH1F*) f->Get(hist_name.c_str()); 
          if(!hist){
            cout<<"\n\nDidn't find histo "<<hist_name<<endl;
            abort();
          }
          tot += hist->GetBinContent(evtType+1);//underflow
        }

      }//loop over subsamples
      cout<<" & "<<tot;
    }//loop over cuts
    cout<<" \\\\ \\hline"<<endl;
  }//loop over samples
}
