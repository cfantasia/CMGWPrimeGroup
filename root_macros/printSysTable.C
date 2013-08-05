//Author: Cory Fantasia 2012
//Purpose: Print latex tables of systematic errors
//Usage: root -b -l -q 'printSysTable.C+("SysSigPDF.dat", "SysBkgPDF.dat", "SysPDF.pdf")'

#include <fstream>
#include "common.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMultiGraph.h"

typedef pair<TGraphErrors*, TGraphErrors*> gPair;

const int nch = 4;
void
printSysTable(string sigName, string bkgName, string outName){
  TCanvas* c = new TCanvas("c", outName.c_str());
  TMultiGraph *mg = new TMultiGraph("mg", ";M_{W', #rho_{TC}} (GeV);Relative Systematic Uncertainty (%)");
  TLegend *leg = new TLegend(0.6, 0.79,0.9, 0.89,"");
  prepLegend(leg);
  leg->SetNColumns(2);
  leg->SetColumnSeparation(0.05);

  TGraphErrors *gSig[nch], *gBkg[nch];
  TGraphErrors *gSigErr[nch], *gBkgErr[nch];

  for(int ch=0; ch<nch; ch++){
    string format = "%lg";//mass
    for(int i=0; i<ch; i++) format += " %*lg %*lg";//value and err
    format += " %lg %lg";
    gSig[ch] = new TGraphErrors(sigName.c_str(), format.c_str()); 
    gBkg[ch] = new TGraphErrors(bkgName.c_str(), format.c_str()); 
    gSigErr[ch] = new TGraphErrors(sigName.c_str(), format.c_str()); 
    gBkgErr[ch] = new TGraphErrors(bkgName.c_str(), format.c_str()); 

    gSig[ch]->SetLineColor(ch+2);
    gBkg[ch]->SetLineColor(ch+2);

    gSig[ch]->SetMarkerColor(ch+2);
    gBkg[ch]->SetMarkerColor(ch+2);

    gSig[ch]->SetMarkerStyle(kOpenCircle);
    gBkg[ch]->SetMarkerStyle(kFullTriangleUp);
    
    
    //Scale the graphs
    float scale = 100.;//100%
    for (int i=0;i<gSig[ch]->GetN();i++){
      gSig[ch]->GetY()[i] *= scale;
      gSig[ch]->GetEY()[i] *= scale;
      gSigErr[ch]->GetY()[i] = gSig[ch]->GetEY()[i];
    }
    
    for (int i=0;i<gBkg[ch]->GetN();i++){
      gBkg[ch]->GetY()[i] *= scale;
      gBkg[ch]->GetEY()[i] *= scale;
      gBkgErr[ch]->GetY()[i] = gBkg[ch]->GetEY()[i];
    }

    mg->Add(gSig[ch]);
    mg->Add(gBkg[ch]);
  
    leg->AddEntry(gSig[ch], Form("Sig %s", bin(ch).c_str()), "EP");
    leg->AddEntry(gBkg[ch], Form("Bkg %s", bin(ch).c_str()), "EP");

  }
  mg->SetMinimum(0.);
  mg->Draw("alp");
  leg->Draw();
  
  //c->SetLogy();
  c->SaveAs(outName.c_str());
  

  //Now Print Latex Table 
  //cout<<"//%Cory: Updated DATE"<<endl;
  cout<<"\\begin{table}[!h] \\centering \\begin{tabular}{|c||*{8}{c|}} \\hline"<<endl;
  cout<<"\\multirow{2}{*}{Mass} & \\multicolumn{4}{c|}{Signal Relative Systematic (\\%)} & \\multicolumn{4}{c|}{Background Relative Systematic (\\%)} \\\\"<<endl;
  for(int ch=0; ch<nch; ch++) cout<<" & "<<bin(ch);//Signal cols
  for(int ch=0; ch<nch; ch++) cout<<" & "<<bin(ch);//Bkg cols
  cout<<"\\\\ \\hline"<<endl;
  for(int mass=200; mass<=1000; mass+=100){
    cout<<Form("%i ",mass);
    for(int ch=0; ch<nch; ch++) cout<<Form(" & %.1f $\\pm$ %.1f", gSig[ch]->Eval(mass), gSigErr[ch]->Eval(mass));
    for(int ch=0; ch<nch; ch++) cout<<Form(" & %.1f $\\pm$ %.1f", gBkg[ch]->Eval(mass), gBkgErr[ch]->Eval(mass));
    cout<<"\\\\"<<endl;
  }
  cout<<"\\hline"<<endl<<"\\end{tabular} "<<endl;
  cout<<"\\caption{"<<outName<<" Systematic uncertainties for signal and background broken down by channel.}"<<endl;
  cout<<"\\label{tab:sys-"<<outName<<"}"<<endl;
  cout<<"\\end{table}"<<endl;

}
