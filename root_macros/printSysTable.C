//Usage: root -b -l -q 'printSysTable.C+("SysSigPDF.dat", "SysBkgPDF.dat", "SysPDF.pdf")'

#include <fstream>
#include "common.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMultiGraph.h"

void
printSysTable(string sigName, string bkgName, string outName){
  
  TGraphErrors* gSig = new TGraphErrors(sigName.c_str(), "%lg %lg %lg"); 
  TGraphErrors* gBkg = new TGraphErrors(bkgName.c_str(), "%lg %lg %lg"); gBkg->SetLineColor(kRed);

  //Scale the graphs
  float scale = 100.;//100%
  for (int i=0;i<gSig->GetN();i++){
    gSig->GetY()[i] *= scale;
    gSig->GetEY()[i] *= scale;
  }

  for (int i=0;i<gBkg->GetN();i++){
    gBkg->GetY()[i] *= scale;
    gBkg->GetEY()[i] *= scale;
  }

  TCanvas* c = new TCanvas("c", outName.c_str());
  TMultiGraph *mg = new TMultiGraph("mg", ";M_{W', #rho_{TC}} (GeV);Relative Systematic Uncertainty");
  mg->Add(gSig,"LE*");
  mg->Add(gBkg,"LE*");
  
  TLegend *leg = new TLegend(0.8, 0.79,0.9, 0.89,"");
  prepLegend(leg);
  leg->AddEntry(gSig, "Signal", "LE");
  leg->AddEntry(gBkg, "Background", "LE");

  //mg->SetMinimum(0.00001);
  mg->Draw("a*");
  leg->Draw();

  c->SaveAs(outName.c_str());


  //Now Print Latex Table 
  //cout<<"//%Cory: Updated DATE"<<endl;
  cout<<"\\begin{table}[!h] \\centering \\begin{tabular}{|c|c|c|} \\hline"<<endl;
  cout<<" Mass & Signal Relative Systematic (%) & Background Relative Systematic (%)\\\\ \\hline"<<endl;
  for(int mass=200; mass<=2000; mass+=100){
    string outString = Form(" %i & %.4f & %.4f \\\\", mass, gSig->Eval(mass), gBkg->Eval(mass));
    cout<<outString<<endl;
  }
  cout<<"\\end{tabular} \\end{table}"<<endl;


}
