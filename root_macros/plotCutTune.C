#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TStyle.h>
#include <TLine.h>

#include <string>
#include <iostream>

#include <set>

using std::string; using std::cout; using std::endl;

enum MASS_POINT {NoSig = 0, OneTeV, OneFiveTeV, TwoTeV};
enum CUT_TYPE {JETVETO = 0, ISO_SUMPT, ISO_DR, CHI2};

const unsigned N_points = 3; // one point below and one above the reference
const unsigned N_masses = 4;  // no signal, 1.0 TeV, 1.5 TeV, 2.0 TeV
string desc_mass[N_masses] = {"No signal", "1.0 TeV", "1.5 TeV", "2.0 TeV"};
string desc_jet[N_points] = {"Jet-Veto 80 GeV", "Jet-Veto 100 GeV", "Jet-Veto 120 GeV"};
string desc_sumpt[N_points] = {"Isolation 3 GeV", "Isolation 10 GeV", "Isolation 25 GeV"};
string desc_dr[N_points] = {"Delta-R = 0.2", "Delta-R = 0.3", "Delta-R = 0.4"};
string desc_chi2[N_points] = {"chi2/Ndof = 3","chi2/Ndof = 10","chi2/Ndof = 100"};

const unsigned N_cuts = 4; 
string desc_cuts[N_cuts] = {"Jet-activity veto", "SumPt isolation", 
			    "Delta-R isolation", "Chi2/Ndof"}; 

TH1F * tot[N_masses][N_points][N_cuts]; // total (sig + bgd) distributions

TFile * file_def = TFile::Open("Wprime_analysis_loose_default.root");
TFile * file_j80 = TFile::Open("Wprime_analysis_loose_JVeto80.root");
TFile * file_j120 = TFile::Open("Wprime_analysis_loose_JVeto120.root");
TFile * file_i3 = TFile::Open("Wprime_analysis_loose_SumPt3.root");  
TFile * file_i25 = TFile::Open("Wprime_analysis_loose_SumPt25.root");  
TFile * file_dr2 = TFile::Open("Wprime_analysis_loose_dR_0.2.root");
TFile * file_dr4 = TFile::Open("Wprime_analysis_loose_dR_0.4.root");
TFile * file_chi2_3 = TFile::Open("Wprime_analysis_loose_chi2_3.root");
TFile * file_chi2_100 = TFile::Open("Wprime_analysis_loose_chi2_100.root");

TFile * file_jet[N_points] = {file_j80, file_def, file_j120};
TFile * file_sumpt[N_points] = {file_i3, file_def, file_i25};
TFile * file_dr[N_points] = {file_dr2, file_def, file_dr4};
TFile * file_chi2[N_points] = {file_chi2_3, file_def, file_chi2_100};

// histograms with signal/background distributions per input file
TH1F * w = 0;
TH1F * qcd = 0;
TH1F * top = 0;
TH1F * z = 0;
TH1F * wp10 = 0;
TH1F * wp15 = 0;
TH1F * wp20 = 0;

void doPlots(unsigned point_j, CUT_TYPE cut_k, MASS_POINT mass_i,TFile * _file0);
void getHistos(TFile * _file0);
bool badHisto(TH1F * h, string s);

// function doPlots gets called several times; need to call Sumw2 on the 
// histograms only once per input file; use files_read_in to keep track of whether
// a file has been read before
std::set<string> files_read;

void plotCutTune()
{
  gStyle->SetOptStat(00000);

  for(unsigned i = 0; i != N_masses; ++i){ // loop over mass points
    
    for (unsigned k = 0; k != N_cuts; ++k){ // loop over cuts

      TCanvas * c1 = new TCanvas();

      for(unsigned j = 0; j != N_points; ++j)
	{
	  if(CUT_TYPE(k) == JETVETO)
	    doPlots(j, JETVETO, MASS_POINT(i), file_jet[j]);
	  else if(CUT_TYPE(k) == ISO_SUMPT)
	    doPlots(j, ISO_SUMPT, MASS_POINT(i), file_sumpt[j]);
	  else if(CUT_TYPE(k) == ISO_DR)
	    doPlots(j, ISO_DR, MASS_POINT(i), file_dr[j]);
	  else if(CUT_TYPE(k) == CHI2)
	    doPlots(j, CHI2, MASS_POINT(i), file_chi2[j]);
	}
      
      char hname0[128]; char hname2[128];
      if(CUT_TYPE(k) == JETVETO)
	{
	  sprintf(hname0, "jet80_%d", i);
	  sprintf(hname2, "jet120_%d", i);
	}
      else if(CUT_TYPE(k) == ISO_SUMPT)
	{
	  sprintf(hname0, "sumpt3_%d", i);
	  sprintf(hname2, "sumpt25_%d", i);
	}
      else if(CUT_TYPE(k) == ISO_DR)
	{
	  sprintf(hname0, "dr6_%d", i);
	  sprintf(hname2, "dr14_%d", i);
	}
      else if(CUT_TYPE(k) == CHI2)
	{
	  sprintf(hname0, "chi2_3_%d", i);
	  sprintf(hname2, "chi2_100_%d", i);
	}

      string htitle = "Relative distribution change for " + 
	desc_cuts[k] + " variations (" + desc_mass[i] + ")";
      TH1F * eff0 = new TH1F(hname0, htitle.c_str(), 
			     tot[i][0][k]->GetNbinsX(),
			     tot[i][0][k]->fXaxis.GetXmin(),
			     tot[i][0][k]->fXaxis.GetXmax());
      TH1F * eff2 = new TH1F(hname2, htitle.c_str(), 
			     tot[i][0][k]->GetNbinsX(),
			     tot[i][0][k]->fXaxis.GetXmin(),
			     tot[i][0][k]->fXaxis.GetXmax());
      c1->Divide(1,2);
      c1->cd(1);
      gPad->SetLogy();
      tot[i][1][k]->Draw();
      c1->cd(2);
      gPad->SetLogy(0);
      eff0->Divide(tot[i][0][k], tot[i][1][k], 1, 1, "B");
      eff2->Divide(tot[i][2][k], tot[i][1][k], 1, 1, "B");
      tot[i][0][k]->SetLineColor(kBlue);
      eff0->SetLineColor(kBlue);
      tot[i][2][k]->SetLineColor(kRed);
      eff2->SetLineColor(kRed);
      eff0->GetXaxis()->SetTitle("Muon p_{T} (GeV/c)");
      eff0->SetMinimum(0.8);
      eff0->SetMaximum(1.2);
      eff0->Draw();
      eff2->Draw("same");
      TLine *line1 = new TLine(eff0->fXaxis.GetXmin(), 1, 
			       eff0->fXaxis.GetXmax(), 1); 
      line1->Draw();

      TLegend * lg = new TLegend(0.59, 0.67, 0.89, 0.89);
      lg->SetTextSize(0.03);
      lg->SetBorderSize(0);
      lg->SetFillColor(0);
      for(int j = N_points - 1; j >= 0; --j)
	{
	  if(CUT_TYPE(k) == JETVETO)
	    lg->AddEntry(tot[i][j][k], desc_jet[j].c_str(), "L");
	  else if(CUT_TYPE(k) == ISO_SUMPT)
	    lg->AddEntry(tot[i][j][k], desc_sumpt[j].c_str(), "L");
	  else if(CUT_TYPE(k) == ISO_DR)
	    lg->AddEntry(tot[i][j][k], desc_dr[j].c_str(), "L");
	  else if(CUT_TYPE(k) == CHI2)
	    lg->AddEntry(tot[i][j][k], desc_chi2[j].c_str(), "L");
	}

      lg->Draw();

      string file(hname0); file += ".gif";
      c1->SaveAs(file.c_str());
      //  delete c1;
      

    } // loop over cuts

  } // loop over mass points

}


void doPlots(unsigned point_j, CUT_TYPE cut_k, MASS_POINT mass_i, TFile * _file0)
{
  unsigned i = unsigned(mass_i);
  unsigned j = point_j;
  unsigned k = unsigned(cut_k);

  getHistos(_file0);
  
  string filename = _file0->GetName();
  // see comment in declartion of files_read
  if(files_read.find(filename) == files_read.end())
    {
      // do this only once per input file 
      wp10->Sumw2(); wp15->Sumw2(); wp20->Sumw2();
      top->Sumw2(); z->Sumw2(); qcd->Sumw2(); w->Sumw2();

      files_read.insert(filename);
    }
      

  string htitle = "Total distribution (" + desc_mass[i] + ")";
  char hname[128];

  if(mass_i == NoSig)
    sprintf(hname, "NoSig_%d_%d", j, k);
  else if(mass_i == OneTeV)
    sprintf(hname, "OneTeV_%d_%d", j, k);
  else if(mass_i == OneFiveTeV)
    sprintf(hname, "OneFiveTeV_%d_%d", j, k);
  else if(mass_i == TwoTeV)
    sprintf(hname, "TwoTeV_%d_%d", j, k);
    
  tot[i][j][k] = new TH1F(hname, htitle.c_str(), w->GetNbinsX(), 
			  w->fXaxis.GetXmin(), w->fXaxis.GetXmax());
  
  if(mass_i == OneTeV)
    tot[i][j][k]->Add(wp10);
  else if(mass_i == OneFiveTeV)
    tot[i][j][k]->Add(wp15);
  else if(mass_i == TwoTeV)
    tot[i][j][k]->Add(wp20);
    

  tot[i][j][k]->Add(top);
  tot[i][j][k]->Add(z);
  tot[i][j][k]->Add(qcd);
  tot[i][j][k]->Add(w);
  tot[i][j][k]->GetXaxis()->SetTitle("Muon p_{T} (GeV/c)");

  // set to zero all entries with (# of ev) < 0.5 to avoid extra large error bars
  const float epsilon = 0.5;
  for(int y = 1; y <= w->GetNbinsX(); ++y)
    if(tot[i][j][k]->GetBinContent(y) < epsilon)
      tot[i][j][k]->SetBinContent(y, 0);

 
}

bool badHisto(TH1F * h, string s)
{
  if(!h)
    {
      cout << " *** Oops! Didn't find " << s << " histogram... " << endl;
      return true;
    }

  return false;
}

void getHistos(TFile * _file0)
{

  string algo = "glb"; // "tev";
  string histo = "hPT" + algo + "_" + "qual";
  string histoW = "W/" + histo;
  w = (TH1F* ) _file0->Get(histoW.c_str());
  if(badHisto(w, "W"))
    return;
  
  string histoQ = "QCD/" + histo;
  qcd = (TH1F* )_file0->Get(histoQ.c_str());
  if(badHisto(qcd, "QCD"))
    return;
  
  string histoZ = "Z/" + histo;
  z = (TH1F* )_file0->Get(histoZ.c_str());
  if(badHisto(z, "Z"))
    return;
  
  string histoT = "Top/" + histo;
  top = (TH1F* )_file0->Get(histoT.c_str());
  if(badHisto(top, "Top"))
    return;
  
  string histo10 = "wprime10/" + histo;
  wp10 = (TH1F* )_file0->Get(histo10.c_str());
  if(badHisto(wp10, "wprime10"))
    return;
  
  string histo15 = "wprime15/" + histo;
  wp15 = (TH1F* )_file0->Get(histo15.c_str());
  if(badHisto(wp15, "wprime15"))
    return;
  
  string histo20 = "wprime20/" + histo;
  wp20 = (TH1F* )_file0->Get(histo20.c_str());
  if(badHisto(wp20, "wprime20"))
    return;
}
