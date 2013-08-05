//Author: Cory Fantasia 2012
//Purpose: Print final background yield plot
//Usage: root -b -l -q 'printYieldPlot.C+("../../../WprimeWZ.root")'

#include <fstream>
#include "../root_macros/common.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"

const int nch = 4;
void
printYieldPlot(const string origName="../../../WprimeWZ.root"){
  TFile *fOrig = TFile::Open(origName.c_str(), "read"); assert(fOrig);

  map<string, TH1F*> hists[nch];

  vector<string> bkgSamples = BkgSamples();
  for(int iBkg=0; iBkg<(int)bkgSamples.size(); ++iBkg){
    for(int ch=0; ch<nch; ++ch){
      hists[ch][bkgSamples[iBkg]] = new TH1F(Form("%s_ch%i", bkgSamples[iBkg].c_str(), ch), "", 40, 0, 2000);
      for(int mass=0; mass<=2000; mass+=50){
        //cout<<bkgSamples[iBkg]<<" "<<ch<<" "<<mass<<endl;

        string cuts = Form("weight*(%s)*(EvtType == %i)", AnalysisCuts(mass).c_str(), ch);
        
        Value origYield = GetNEvtsAndError(fOrig, bkgSamples[iBkg], "tEvts_MET", cuts);
        
        hists[ch][bkgSamples[iBkg]]->Fill(mass, origYield.val);
      }
    }
  }
  //Data
  for(int ch=0; ch<nch; ++ch){
    hists[ch]["data"] = new TH1F(Form("%s_ch%i", "data", ch), "", 40, 0, 2000);
    for(int mass=0; mass<=2000; mass+=50){
      //cout<<"data"<<" "<<ch<<" "<<mass<<endl;
      
      string cuts = Form("weight*(%s)*(EvtType == %i)", AnalysisCuts(mass).c_str(), ch);
      
      Value origYield = GetNEvtsAndError(fOrig, "data", "tEvts_MET", cuts);
      
      hists[ch]["data"]->Fill(mass, origYield.val);
    }
  }
      

  TCanvas* c = new TCanvas("c", "");
  c->Divide(2,2);

  TLegend *leg = new TLegend(0.5, 0.55,0.9, 0.89,"");
  prepLegend(leg);
  leg->SetNColumns(2);
  leg->SetColumnSeparation(0.05);

  THStack* stacks[nch];
  for(int ch=0; ch<nch; ++ch){
    c->cd(1+ch)->SetLogy();
    stacks[ch] = new THStack(Form("BkgCh%i", ch), Form("Background Yield By Mass Point (%s)", binLatex(ch).c_str()));
    for(int iBkg=bkgSamples.size()-1; iBkg>=0; --iBkg){
      int fillColor;
      if(iBkg==0) fillColor = kOrange-2;
      if(iBkg==1) fillColor = kOrange+7;
      if(iBkg==2) fillColor = kViolet+2;
      if(iBkg==3) fillColor = kGray;
      if(iBkg==4) fillColor = kOrange+3;

      string legName;
      if(iBkg==0) legName = "WZ";
      if(iBkg==1) legName = "Z+Jets";
      if(iBkg==2) legName = "t\\bar{t}";
      if(iBkg==3) legName = "Z#gamma";
      if(iBkg==4) legName = "ZZ";

      hists[ch][bkgSamples[iBkg]]->SetFillColor(fillColor);
      //hists[ch][bkgSamples[iBkg]]->SetLineColor(fillColor);
      stacks[ch]->Add(hists[ch][bkgSamples[iBkg]]);
      if(ch==0) leg->AddEntry(hists[ch][bkgSamples[iBkg]], legName.c_str(), "F");
    }
    //Data
    hists[ch]["data"]->SetMarkerStyle(20);
    //hists[ch]["data"]->SetLineColor(fillColor);
    if(ch==0) leg->AddEntry(hists[ch]["data"], "Data", "PE");
    

    stacks[ch]->Draw("hist");
    stacks[ch]->SetMaximum(50);
    stacks[ch]->SetMinimum(1e-2);
    hists[ch]["data"]->Draw("E SAME");
    leg->Draw();
  }
  
  c->SaveAs("BkgYieldByMassPoint.png");

}
