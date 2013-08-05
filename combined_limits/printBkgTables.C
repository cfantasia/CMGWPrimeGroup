/* 
   Author: Cory Fantasia
   Script to print important numbers in LaTex form
   To Use: 
   root -b -l -q 'printBkgTable.C+("../../../WprimeWZ.root")'
*/

//#include "consts.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "../root_macros/common.h"
#include "makeLimitCard.C"
#include "TGraphAsymmErrors.h"
#include "../root_macros/CMSStyle.C"
#include "TMultiGraph.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLegend.h"

void
printBkgTables(string inName="../../../WprimeWZ.root"){
  TFile *fIn = TFile::Open(inName.c_str(), "read"); assert(fIn);
  //Make Sig Eff Graph
  TGraphErrors* geff   [nch];
  TGraphErrors* gefferr[nch];
  for(int ch=0; ch<nch; ch++){
    geff[ch]    = MakeSigEffGraph   (fIn, Form("EvtType == %i", ch));
    gefferr[ch] = MakeSigEffErrGraph(fIn, Form("EvtType == %i", ch));
  }

  TGraph* gxsec = new TGraph("../Limits/xSec_WZ.dat", "%lg %lg");
  TGraph* gxsecCh = new TGraph(*gxsec);
  Scale(gxsecCh, 1./4.);

  vector<string> bkgSamples = BkgSamples();

  //map<string, Value> sysErr[4];
  for(int ch=0; ch<nch; ++ch){
    //Open limits file
    TTree* tLimits = new TTree("tLimits", "Limits");
    tLimits->ReadFile(Form("nLimit_WprimeWZCh%i_MarkovChainMC.txt",ch));

    cout<<"\\begin{table}[!htbp] \\centering \\begin{tabular}{|c||*{8}{c|}} \\hline"<<endl;
    cout<<"Mass point & WZ & ZZ & Z\\gamma & DY & TT & $Sum_{MC}$ \\\\"<<endl;
  
    map<string, Value> yields;
    float nLims = tLimits->Draw("Mass:ExpLimit:ObsLimit", "1 || Mass%100 == 0", "para goff");
    for(int iLims=0; iLims<nLims; ++iLims){
      const int      mass   = tLimits->GetVal(0)[iLims];
      const Double_t expLim = tLimits->GetVal(1)[iLims];
      const Double_t obsLim = tLimits->GetVal(2)[iLims];

      string analysisCuts = AnalysisCuts(mass);
      double minLt = LtCut(mass);
      float winWidth = WindowWidth(mass);
      float minWind = mass - winWidth/2.;
      float maxWind = mass + winWidth/2.;
      
      string chCuts = Form("EvtType == %i", ch);
      string cuts = "(" + analysisCuts + " && " + chCuts + ") * weight";

      yields = getYields(fIn, mass, cuts, bkgSamples, gxsecCh, geff[ch], gefferr[ch]);
      //sysErr[ch] = getSysErrs(mass, bkgSamples, yields, ch);

      cout<<setiosflags(ios::fixed) << setprecision(0)
          <<"\\PWpr "<<mass;
      const Value & vBkg = yields["bkg"];
      const Value & vWZ  = yields["WZJetsTo3LNu"];
      const Value & vZZ  = yields["ZZ"];
      const Value & vZG  = yields["GVJets"];
      const Value & vTT  = yields["TTJets"];
      const Value & vDY  = yields["DYJetsToLL"];
      
      cout<<std::fixed << std::setprecision(1)
          <<" & "<<vWZ<<" $\\pm$ "<<Value(vWZ.err)
          <<" & "<<vZZ<<" $\\pm$ "<<Value(vZZ.err)
          <<" & "<<vZG<<" $\\pm$ "<<Value(vZG.err)
          <<" & "<<vDY<<" $\\pm$ "<<Value(vDY.err)
          <<" & "<<vTT<<" $\\pm$ "<<Value(vTT.err)
          <<" & "<<vBkg<<" $\\pm$ "<<Value(vBkg.err);//right one 
      cout<<" \\\\" <<endl;
 
   
    }//Mass Loop
    cout<<"\\hline"<<endl<<"\\end{tabular} "<<endl;
    cout<<"\\caption{Search windows for each \\PWpr mass point giving the number of expected background events from Monte Carlo, number of events observed, number of expected signal events from MC, signal efficiency, and the expected and"<<endl
        <<"observed exclusion limits on $\\sigma \\times BR(\\PWpr \\rightarrow"<<endl
        <<"3 \\ell \\nu)$ at 95\% C.L. for the $"<<3-ch<<"e"<<ch<<"\\mu$ channel. Errors indicated are statistical only.}"<<endl;
    cout<<"\\label{tab:BackgroundWindow-"<<3-ch<<"e"<<ch<<"mu"<<endl;
    cout<<"\\end{table}"<<endl;
    delete tLimits;

  }//ch loop

  
  //Now do sum table
  map<string, Value> yieldsSum;//, sysErrSum;

  TTree* tLimits = new TTree("tLimits", "Limits");
  tLimits->ReadFile("nLimit_WprimeWZ_MarkovChainMC.txt");
  
  cout<<"\\begin{table}[!htbp] \\centering \\begin{tabular}{|c||*{8}{c|}} \\hline"<<endl;
  cout<<"Mass point & WZ & ZZ & Z\\gamma & DY & TT & $Sum_{MC}$ \\\\"<<endl;
  
  float nLims = tLimits->Draw("Mass:ExpLimit:ObsLimit", "1 || Mass%100 == 0", "para goff");
  for(int iLims=0; iLims<nLims; ++iLims){
    const int      mass   = tLimits->GetVal(0)[iLims];
    const Double_t expLim = tLimits->GetVal(1)[iLims];
    const Double_t obsLim = tLimits->GetVal(2)[iLims];
    string analysisCuts = AnalysisCuts(mass);
    double minLt = LtCut(mass);
    float winWidth = WindowWidth(mass);
    float minWind = mass - winWidth/2.;
    float maxWind = mass + winWidth/2.;
    
    string chCuts = Form("EvtType > -1");
    string cuts = "(" + analysisCuts + " && " + chCuts + ") * weight";
    
    TGraphErrors* geffSum      = MakeSigEffGraph   (fIn, "1", 0, 0, 4.);
    TGraphErrors* gefferrSum   = MakeSigEffErrGraph(fIn, "1", 0, 0, 4.);


    yieldsSum = getYields(fIn, mass, cuts, bkgSamples, gxsec, geffSum, gefferrSum);
    //sysErrSum = getSysErrSum(sysErr, bkgSamples);
    
    cout<<setiosflags(ios::fixed) << setprecision(0)
        <<"\\PWpr "<<mass;
    const Value & vBkg = yieldsSum["bkg"];
    const Value & vWZ  = yieldsSum["WZJetsTo3LNu"];
    const Value & vZZ  = yieldsSum["ZZ"];
    const Value & vZG  = yieldsSum["GVJets"];
    const Value & vTT  = yieldsSum["TTJets"];
    const Value & vDY  = yieldsSum["DYJetsToLL"];
    
    cout<<std::fixed << std::setprecision(1)
        <<" & "<<vWZ<<" $\\pm$ "<<Value(vWZ.err)
        <<" & "<<vZZ<<" $\\pm$ "<<Value(vZZ.err)
        <<" & "<<vZG<<" $\\pm$ "<<Value(vZG.err)
        <<" & "<<vDY<<" $\\pm$ "<<Value(vDY.err)
        <<" & "<<vTT<<" $\\pm$ "<<Value(vTT.err)
        <<" & "<<vBkg<<" $\\pm$ "<<Value(vBkg.err);//right one 
    cout<<" \\\\" <<endl;
  }//Mass Loop
  cout<<"\\hline"<<endl<<"\\end{tabular}"<<endl;
  cout<<"\\caption{Search windows for each \\PWpr mass point giving the number of expected background events from Monte Carlo, number of events observed, number of expected signal events from MC, signal efficiency, and the expected and"<<endl
      <<"observed exclusion limits on $\\sigma \\times BR(\\PWpr \\rightarrow"<<endl
      <<"3 \\ell \\nu)$ at 95\% C.L. for the combined channels. Errors indicated are statistical only.}"<<endl;
  cout<<"\\label{tab:BackgroundWindow-combined}"<<endl;
  cout<<"\\end{table}"<<endl;
  delete tLimits;

    ////////////////


}
