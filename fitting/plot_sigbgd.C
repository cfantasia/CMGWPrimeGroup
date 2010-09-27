#include <TFile.h>
#include <TH1F.h>
#include <TPad.h>
#include <TF1.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveStats.h>

#include <string>
#include <iostream>

using std::cout; using std::endl; using std::string;

unsigned mass_option = 1; // 1, 2 -> 1.0, 1.5
unsigned f_option = 2; // 1->fix at 1.0, 2->fix at 2.0, 3->float

void plot_sigbgd()
{
  assert(mass_option >= 1 && mass_option <= 2);
  assert(f_option >= 1 && f_option <= 3);
  
  gStyle->SetOptStat(0);

  TFile * file = TFile::Open("Wprime_analysis.root");
  string histo_lumi = "lumi_ipb";
  TH1F * lumi_ipb = (TH1F* )file->Get(histo_lumi.c_str());
  float Lumi_ipb = lumi_ipb->GetBinContent(1);
  file->Close();

  char filein[256]; char fileout[256];
  string energy = ""; string fudge = "";
  string en = "";
  if(mass_option == 1)
    { energy = "oneTeV"; en = "1.0";}
  else if(mass_option == 2)
    { energy = "onefiveTeV"; en = "1.5";}

  if(f_option == 1) 
    fudge = "f_1.0";
  else if(f_option == 2) 
    fudge = "f_2.0";
  else if(f_option == 3) 
    fudge = "f_float";

  //  sprintf(filein, "fit_sigbgd_%s_%s.root", energy.c_str(), fudge.c_str());
  sprintf(filein, "fit_sigbgd.root");
  sprintf(fileout, "fit_sigbgd_%s_%s_log.eps",energy.c_str(), fudge.c_str());
  //sprintf(fileout, "fit_sigbgd_%s_%s_log.pdf",energy.c_str(), fudge.c_str());

  TFile * f = new TFile(filein);
  if(f->IsZombie())
    {
      cout << " Oops! Cannot open file " << filein << endl;
      return;
    }

  TH1F * tot = (TH1F *) f->Get("tot");
  if(!tot)
    {
      cout << " Oops! Cannot open histogram tot " << endl;
      return;
    }
  tot->SetTitle("");

  tot->SetMarkerStyle(4); tot->SetMarkerSize(1.2);


  tot->Draw(); gPad->Update();

  char desc[512]; sprintf(desc, " %4.2f pb^{-1}", Lumi_ipb);
  char desc1[512]; 
  if(f_option == 1)
    sprintf(desc1, "f = 1.0");
  else if(f_option == 2)
    sprintf(desc1, "f = 2.0");

  float xo = -9; float yo = -9;
  float x1 = -9; float y1 = -9;
  if(mass_option == 1)
    {
      xo = 195; yo = 1;
      x1 = 200; y1 = 0.6;
    }
  else if(mass_option == 2)
    {
      xo = 230; yo = 90;
      x1 = 235; y1 = 50;
    }

  TLatex * lo = new TLatex(xo, yo, desc);
  lo->SetTextSize(0.04);

  TCanvas * c1 = new TCanvas();
  c1->SetLogy(1);

  tot->Draw("e");
  lo->Draw();

  c1->SaveAs(fileout);



}
