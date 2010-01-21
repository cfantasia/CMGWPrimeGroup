#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TStyle.h>
#include <TLine.h>

#include <string>
#include <iostream>

using std::string; using std::cout; using std::endl;

enum MASS_POINT {NoSig = 0, OneTeV, OneFiveTeV, TwoTeV};
enum ALGO_TYPE {GLOBAL = 0, TRACKER, TRK_1ST};

const unsigned N_algos = 3;
string algos[N_algos] = {"glb", "trk", "tev"};
string desc_algo[N_algos] = {" Global Muons", "Tracker-only Muons",
			"Tracker + 1st station Muons"};

const unsigned N_masses = 4;  // no signal, 1.0 TeV, 1.5 TeV, 2.0 TeV
string desc_mass[N_masses] = {"No signal", "1.0 TeV", "1.5 TeV", "2.0 TeV"};

// histograms with signal/background distributions
TH1F * w = 0;
TH1F * qcd = 0;
TH1F * top = 0;
TH1F * z = 0;
TH1F * wp10 = 0;
TH1F * wp15 = 0;
TH1F * wp20 = 0;

void getHistos(TFile * _file0, string algo, bool plus_flag);
void doPlots(TFile * _file0, MASS_POINT mass_i, unsigned algo_j, bool plus_flag);

// total (sig + bgd) distributions
TH1F * tot_plus[N_algos] = {0};
TH1F * tot_minus[N_algos] = {0};
TH1F * asym[N_algos] = {0};
TH1F * tot[N_algos] = {0};

void plotChargeAsym()
{
  string input_file = "Wprime_analysis.root";
  TFile *_file0 = TFile::Open(input_file.c_str());
  if(!_file0 || _file0->IsZombie())
    {
      cout << " *** Ooops! Cannot find input file " << input_file << endl;
      return;
    }

  gStyle->SetOptStat(00000);

  for(unsigned i = 0; i != N_masses; ++i) // loop over mass points
    {
      for(unsigned j = 0; j != N_algos; ++j)
	{  // loop over algorithms
	  doPlots(_file0, MASS_POINT(i), j, true); // positive charge
	  doPlots(_file0, MASS_POINT(i), j, false); // negative charge

	  char hname[128]; string htitle;

	  sprintf(hname, "tot_%s_%d", algos[j].c_str(), i);
	  htitle = "Total distribution for " + desc_algo[j] + " (" + 
	    desc_mass[i] + ")";
	  tot[j] = new TH1F(hname, htitle.c_str(), 
			    tot_plus[j]->GetNbinsX(),
			    tot_plus[j]->fXaxis.GetXmin(),
			    tot_plus[j]->fXaxis.GetXmax());
 	  tot[j]->Add(tot_plus[j], tot_minus[j]);
	  sprintf(hname, "asym_%s_%d", algos[j].c_str(), i);
	  htitle = "Charged asymmetry for " + desc_algo[j] + " (" + 
	    desc_mass[i] + ")";
	  asym[j] = (TH1F *) tot_plus[j]->GetAsymmetry(tot_minus[j]);
	  asym[j]->SetName(hname);
	  asym[j]->SetTitle(htitle.c_str());

	  TCanvas * c1 = new TCanvas();
	  c1->Divide(1,2);
	  c1->cd(1);
	  gPad->SetLogy();
	  tot[j]->Draw();
	  c1->cd(2);
	  asym[j]->Draw();

	  TLine *line1 = new TLine(asym[j]->fXaxis.GetXmin(), 0, 
				   asym[j]->fXaxis.GetXmax(), 0); 
	  line1->SetLineColor(kRed);
	  line1->Draw();

	  string file(hname); file += ".gif";
	  c1->SaveAs(file.c_str());

	  // delete c1;

	} // loop over algorithms


    } // loop over mass points

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

void doPlots(TFile * _file0, MASS_POINT mass_i, unsigned algo_j, bool plus_flag)
{

  getHistos(_file0, algos[algo_j], plus_flag);

  // do this only once per file (therefore: for one mass point)
  if(mass_i == NoSig)
    {
      wp10->Sumw2(); wp15->Sumw2(); wp20->Sumw2();
      top->Sumw2(); z->Sumw2(); qcd->Sumw2(); w->Sumw2();
    }
  
  char hname[128]; string htitle;
  if(plus_flag)
    {
      sprintf(hname, "plus_%s_%d", algos[algo_j].c_str(), int(mass_i));
      htitle = "Positively";
    }
  else
    {
      sprintf(hname, "minus_%s_%d", algos[algo_j].c_str(), int(mass_i));
      htitle = "Negatively";
    }
  htitle += " charged distribution for " + desc_algo[algo_j] + " (" + 
    desc_mass[int(mass_i)] + ")";
  TH1F * htemp = new TH1F(hname, htitle.c_str(), w->GetNbinsX(), 
			  w->fXaxis.GetXmin(), w->fXaxis.GetXmax());

  htemp->Add(top);
  htemp->Add(z);
  htemp->Add(qcd);
  htemp->Add(w);
 
  if(mass_i == OneTeV)
    htemp->Add(wp10);
  else if(mass_i == OneFiveTeV)
    htemp->Add(wp15);
  else if(mass_i == TwoTeV)
    htemp->Add(wp20);

  htemp->GetXaxis()->SetTitle("Muon p_{T} (GeV/c)");

  // set to zero all entries with (# of ev) < 0.5 to avoid extra large error bars
  const float epsilon = 0.5;
  for(int y = 1; y <= w->GetNbinsX(); ++y)
    if(htemp->GetBinContent(y) < epsilon)
      htemp->SetBinContent(y, 0);

  if(plus_flag)
    tot_plus[algo_j] = htemp;
  else
    tot_minus[algo_j] = htemp;

}

void getHistos(TFile * _file0, string algo, bool plus_flag)
{

  string postfix;
  if(plus_flag)
    postfix = "plus";
  else
    postfix = "minus";

  string histo = "hPT" + algo + "_" + postfix;
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
