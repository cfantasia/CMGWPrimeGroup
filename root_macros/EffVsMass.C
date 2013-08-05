//Author: Cory Fantasia 2012
//Purpose: Plot efficiency versus W' mass
//Usage: root -b -l -q 'EffVsMass.C+(string inName, int applyAnalysisCuts)'
//e.g. : root -b -l -q 'EffVsMass.C+("../../../WprimeWZ.root", 1)'

#include <iostream>
#include <iomanip>

#include "TLatex.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TError.h"
#include "TStyle.h"
#include "TLine.h"

#include "common.h"
#include "CMSStyle.C"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"

using namespace std;

void EffVsMass(string inName, int applyAnalysisCuts){

  gErrorIgnoreLevel = kWarning;
  CMSstyle();
  gStyle->SetOptStat(0);

  TFile *fin = TFile::Open(inName.c_str(), "read"); assert(fin);
  TGraphAsymmErrors* g[4];
  for(int ch=0; ch<4; ++ch){
    g[ch] = new TGraphAsymmErrors();
  }

  int ipoint = 0;
  for(int mass=200; mass <= 2000; mass += 100){
    string name = Form("WprimeToWZTo3LNu_M-%i",mass);
    
    float sample_tot(0), sample_pass(0);
    
    TH1F* hInfo = (TH1F*) fin->Get(Form("%s/hFileInfo", name.c_str()));
    if(!hInfo) continue;  
    float sample_weight = GetSampleInfo(hInfo, "Sample Weight");
    sample_tot = GetSampleInfo(hInfo, "Number of Events Produced");
    sample_tot *= (4./9.); //remove tau contribution from denom
    
    //get tree from file
    TTree* t = getTree(fin, name, "tEvts_MET"); assert(t);
    string analysisCuts = applyAnalysisCuts ? (applyAnalysisCuts == 1 ? Form("Lt > %.0f", LtCut(mass)) : AnalysisCuts(mass)) : "1";

    for(int ch=0; ch<4; ++ch){
      float tot = sample_tot / 4;
      float pass = 0;

      float pass_alt = t->Draw("WZMass", Form("weight*(EvtType == %i)*(%s)",ch, analysisCuts.c_str()), "goff");
      int n = t->GetSelectedRows();
      for(int ientry=0; ientry<n; ++ientry){
        float weight = t->GetW()[ientry];
        pass += weight;
      }
      pass /= sample_weight;
      sample_pass += pass;
      
      float mean    = tot>0 ? pass / tot : 0.;
      float errUp   = TEfficiency::ClopperPearson(tot, pass, 0.68, true) - mean;
      float errDown = mean - TEfficiency::ClopperPearson(tot, pass, 0.68, false);
      
      errDown = std::max(errDown, (float) 0.);
      float err = (errUp + errDown)/2;
      printf(" %s %ie%i\\mu &  %4.2f & %4.2f or %4.2f & %.2f \\pm %.2f\\%% \\\\\n", name.c_str(), 3-ch, ch, tot, pass, pass_alt, mean*100, err*100);

      g[ch]->SetPoint(ipoint, mass, mean);
      g[ch]->SetPointError(ipoint, 0., 0., errDown, errUp);
    }
    ipoint++;
  }

  TCanvas* c1 = new TCanvas();
  TMultiGraph* mg = new TMultiGraph("mg", ";M_{W'} (GeV); #varepsilon");
  TLegend *legend = new TLegend(0.4,0.2,0.51,0.4,"");
  prepLegend(legend);
  for(int ch=0; ch<4; ++ch){
    g[ch]->SetMarkerColor(ch+2);
    g[ch]->SetLineColor(ch+2);
    legend->AddEntry(g[ch], Form("%ie%i#mu", 3-ch, ch), "PE");
    mg->Add(g[ch]);
  }
  mg->SetMaximum(1.1);
  mg->SetMinimum(0. );
  mg->Draw("AP");
  legend->Draw();

  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextSize(0.05);
  latexLabel.SetTextFont(42);
  latexLabel.DrawLatex(0.33, 0.96, "CMS Simulation 2012");
  latexLabel.DrawLatex(0.21, 0.85, "#sqrt{s} = 8 TeV");

  string outName = applyAnalysisCuts ? (applyAnalysisCuts==1 ? "EffVsMass_LtCut.png" : "EffVsMass_AnalysisCuts.png") : "EffVsMass.png";
  c1->SaveAs(outName.c_str());
}


