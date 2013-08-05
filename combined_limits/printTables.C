/* 
   Author: Cory Fantasia
   Script to print important numbers in LaTex form
   To Use: 
   root -b -l -q 'printTable.C+("../../../WprimeWZ.root")'
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
printTables(string inName="../../../WprimeWZ.root"){
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

  TGraphAsymmErrors *gData[4], *gBkg[4];
  for(int ch=0; ch<4; ++ch){
    gData[ch] = new TGraphAsymmErrors();
    gBkg [ch] = new TGraphAsymmErrors();
  }

  map<string, Value> sysErr[4];
  for(int ch=0; ch<nch; ++ch){
    //Open limits file
    TTree* tLimits = new TTree("tLimits", "Limits");
    tLimits->ReadFile(Form("nLimit_WprimeWZCh%i_MarkovChainMC.txt",ch));

    cout<<"\\begin{table}[!htbp] \\centering \\begin{tabular}{|c||*{8}{c|}} \\hline"<<endl;
    cout<<"Mass point & $L_{T}$ & Window & $N_{Bkg}^{MC}$ & $N_{Data}$ & $N_{Sig}$ & $\\varepsilon_{Sig}$ & $\\sigma^{Exp.}_{Limit}$ & $\\sigma^{Obs.}_{Limit}$ \\\\"<<endl;
  
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
      sysErr[ch] = getSysErrs(mass, bkgSamples, yields, ch);

      cout<<setiosflags(ios::fixed) << setprecision(0)
          <<"\\PWpr "<<mass
          <<" & "<<minLt
          <<" & "<<minWind<<"-"<<maxWind;
      const Value & vBkg = yields["bkg"];
      const Value & vData = yields["data"];
      const Value & vSig = yields["sig"];
      const Value & vSigEff = yields["sigeff"];
      
      const Value & vBkgSys = sysErr[ch]["bkg"];
      const Value & vSigEffSys = sysErr[ch]["sigeff"];
      cout<<std::fixed << std::setprecision(1)
          <<" & "<<vBkg<<" $\\pm$ "<<Value(vBkg.err)//right one 
        //<<" & "<<vBkg<<" $\\pm$ "<<Value(0.)//delete me
        //<<" $\\pm$ "<<Value(vBkgSys.err)//Comment this line to hide Sys Error
        //<<" $=$ "<<Value(AddInQuad(vBkg.err,vBkgSys.err))
          <<" & "<<(int)vData.val
//          <<std::fixed << std::setprecision(1)
          <<" & "<<vSig<<" $\\pm$ "<<Value(vSig.err)
//          <<std::fixed << std::setprecision(1)
          <<" & "<<vSigEff*100<<" $\\pm$ "<<Value(vSigEff.err*100)//right one
        //<<" & "<<vSigEff*100<<" $\\pm$ "<<Value(0.)//delete me 
        //<<" $\\pm$ "<<Value(vSigEffSys.err*100)//Comment this line to hide Sys Errors
        //<<" $=$ "<<Value(AddInQuad(vSigEff.err*100,vSigEffSys.err*100))
          <<" & "<<Value(expLim, -2)
          <<" & "<<Value(obsLim, -2);
      cout<<" \\\\" <<endl;
 
    
      ///////Make Plots
      if(vData.val) gData[ch]->SetPoint     (iLims, mass, vData.val);
      
      if(vBkg.val){
        gBkg[ch]->SetPoint     (iLims, mass, vBkg.val);
        gBkg[ch]->SetPointError(iLims, 0., 0., vBkg.err, vBkg.err);
      }

    }//Mass Loop
    cout<<"\\hline"<<endl<<"\\end{tabular} "<<endl;
    cout<<"\\caption{Search windows for each \\PWpr mass point giving the number of expected background events from Monte Carlo, number of events observed, number of expected signal events from MC, signal efficiency, and the expected and"<<endl
        <<"observed exclusion limits on $\\sigma \\times BR(\\PWpr \\rightarrow"<<endl
        <<"3 \\ell \\nu)$ at 95\% C.L. for the $"<<3-ch<<"e"<<ch<<"\\mu$ channel. Errors indicated are statistical only.}"<<endl;
    cout<<"\\label{tab:MassWindows-combined}"<<endl;
    cout<<"\\end{table}"<<endl;
    delete tLimits;

  }//ch loop

  
  //Now do sum table
  map<string, Value> yieldsSum, sysErrSum;

  TTree* tLimits = new TTree("tLimits", "Limits");
  tLimits->ReadFile("nLimit_WprimeWZ_MarkovChainMC.txt");
  
  cout<<"\\begin{table}[!htbp] \\centering \\begin{tabular}{|c||*{8}{c|}} \\hline"<<endl;
  cout<<"Mass point & $L_{T}$ & Window & $N_{Bkg}^{MC}$ & $N_{Data}$ & $N_{Sig}$ & $\\varepsilon_{Sig}$ & $\\sigma^{Exp.}_{Limit}$ & $\\sigma^{Obs.}_{Limit}$ \\\\"<<endl;
  
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
    sysErrSum = getSysErrSum(sysErr, bkgSamples);
    
    cout<<setiosflags(ios::fixed) << setprecision(0)
        <<"\\PWpr "<<mass
        <<" & "<<minLt
        <<" & "<<minWind<<"-"<<maxWind;
    const Value & vBkg    = yieldsSum["bkg"];
    const Value & vData   = yieldsSum["data"];
    const Value & vSig    = yieldsSum["sig"];
    const Value & vSigEff = yieldsSum["sigeff"];

    const Value & vBkgSys = sysErrSum["bkg"];
    const Value & vSigEffSys = sysErrSum["sigeff"];
    cout<<std::fixed << std::setprecision(1)
        <<" & "<<vBkg<<" $\\pm$ "<<Value(vBkg.err) //right one
      //<<" & "<<vBkg<<" $\\pm$ "<<Value(0.)//delrte me 
      //<<" $\\pm$ "<<Value(vBkgSys.err)//Comment this line to hide Sys Error
      //<<" $=$ "<<Value(AddInQuad(vBkg.err,vBkgSys.err))
        <<" & "<<(int)vData.val
//          <<std::fixed << std::setprecision(1)
        <<" & "<<vSig<<" $\\pm$ "<<Value(vSig.err)
//          <<std::fixed << std::setprecision(1)
        <<" & "<<vSigEff*100<<" $\\pm$ "<<Value(vSigEff.err*100)//right one
      //<<" & "<<vSigEff*100<<" $\\pm$ "<<Value(0.)//delete me
      //<<" $\\pm$ "<<Value(vSigEffSys.err*100/4.)//Comment this line to hide Sys Error (Factor of 4 for naive combo of channels)
      //<<" $=$ "<<Value(AddInQuad(vSigEff.err*100,vSigEffSys.err*100/4.))
        <<" & "<<Value(expLim, -2)
        <<" & "<<Value(obsLim, -2);
    cout<<" \\\\" <<endl;
  }//Mass Loop
  cout<<"\\hline"<<endl<<"\\end{tabular}"<<endl;
  cout<<"\\caption{Search windows for each \\PWpr mass point giving the number of expected background events from Monte Carlo, number of events observed, number of expected signal events from MC, signal efficiency, and the expected and"<<endl
      <<"observed exclusion limits on $\\sigma \\times BR(\\PWpr \\rightarrow"<<endl
      <<"3 \\ell \\nu)$ at 95\% C.L. for the combined channels. Errors indicated are statistical only.}"<<endl;
  cout<<"\\label{tab:MassWindows-combined}"<<endl;
  cout<<"\\end{table}"<<endl;
  delete tLimits;

    ////////////////
    //Finalize Plots
    ////////////////

  gErrorIgnoreLevel = kWarning;
  CMSstyle();
  gStyle->SetOptStat(0);

  TCanvas* c1 = new TCanvas();
  c1->Divide(2,2);
  TMultiGraph* mg[4];
  TLegend *legend = new TLegend(0.4,0.2,0.51,0.4,"");
  prepLegend(legend);
  for(int ch=0; ch<nch; ++ch){
    c1->cd(ch+1);
    c1->SetLogy(1);
    mg[ch] = new TMultiGraph("mg", ";M_{W'} (GeV); N_{Evts}");

    gBkg[ch]->SetMarkerStyle(kFullTriangleUp);
    gBkg[ch]->SetMarkerColor(ch+2);
    gBkg[ch]->SetLineColor(ch+2);
    legend->AddEntry(gBkg [ch], Form("Bkg  %ie%i#mu", 3-ch, ch), "PE");
    mg[ch]->Add(gBkg [ch]);

    gData[ch]->SetMarkerStyle(kFullCircle);
    gData[ch]->SetMarkerColor(ch+2);
    gData[ch]->SetLineColor(ch+2);
    legend->AddEntry(gData[ch], Form("Data %ie%i#mu", 3-ch, ch), "PE");
    mg[ch]->Add(gData[ch]);

    mg[ch]->SetMinimum(0. );
    mg[ch]->Draw("AP");
    legend->Draw();
  }
  //mg->SetMaximum(1.1);

  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextSize(0.05);
  latexLabel.SetTextFont(42);
  latexLabel.DrawLatex(0.33, 0.96, "CMS Preliminary 2012");
  latexLabel.DrawLatex(0.21, 0.85, "#sqrt{s} = 8 TeV");

  //string outName = applyAnalysisCuts ? (applyAnalysisCuts==1 ? "EffVsMass_LtCut.png" : "EffVsMass_AnalysisCuts.png") : "EffVsMass.png";
  string outName = "NEvtsVsMass.png";
  c1->SaveAs(outName.c_str());


}
