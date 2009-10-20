// ROOT stuff
#include <TROOT.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1F.h>

// wprime stuff
#include "UserCode/CMGWPrimeGroup/interface/wprimeEvent.h"
#include "UserCode/CMGWPrimeGroup/root_macros/input_file.h"
#include "UserCode/CMGWPrimeGroup/root_macros/loadCuts.h"
#include "UserCode/CMGWPrimeGroup/root_macros/util.h"

// std stuff
#include <iostream>

using std::cout; using std::endl; using std::vector;

TH1F * hPTplus[Num_trkAlgos] = {0};
TH1F * hPTminus[Num_trkAlgos] = {0};
//TH1F * hPTasym[Num_trkAlgos] = {0};

static float Nexp_evt = 0;
static float Nexp_evt_cut[Num_trkAlgos] = {0};

static void initHistos();
static void saveHistos(TFile * fout, string dir);
static void printSummary();
static void fillHistos(const wprime::Muon * mu, float weight, 
		       bool * survived_cut, bool glb_flag, bool trk_flag, 
		       bool tev_flag);

void GetChargePtDistribution(const vector<wprime::InputFile>& files, 
			   TFile *fout, string dir)
{
  initHistos();
  int Nfiles = files.size();

  // counter initialization
  Nexp_evt = 0;
  for(int j = 0; j != Num_trkAlgos; ++j)
    Nexp_evt_cut[j] = 0;

  for(int tr = 0; tr != Nfiles; ++tr){
    if(!files[tr].tree)
      continue;
    
    cout << " Processing sample " << files[tr].description << endl;
    wprime::Event * ev = new wprime::Event();
    files[tr].tree->SetBranchAddress("wp", &ev);
    
    int nevents = files[tr].tree->GetEntries();
    float weight = files[tr].weight;

    int Num_surv_cut[Num_trkAlgos] = {0};
    for(int i = 0; i != nevents; ++i)
      { // event loop
	files[tr].tree->GetEntry(i);
	
	bool survived_cut[Num_trkAlgos] = {false};
	
       	int nmuons = ev->mu->GetLast() + 1;

	bool HLT = ev->HLT_Mu9;
	bool OneMuon = (NmuAboveThresh(OneMuPtTrackCut, ev) == 1);

	//	bool JetActivity = (NjetAboveThresh(EtJetCut, ev) > 0);

       	for(int j = 0; j != nmuons; ++j)
	  { // loop over muons
	    wprime::Muon * mu = (wprime::Muon *) ev->mu->At(j);

	    // trigger decision
 	    if(!HLT)continue;

	    bool ptrange_glb = mu->global.p.Pt() >=minPtMu 
	      && mu->global.p.Pt() <=maxPtMu;
	    bool ptrange_trk = mu->tracker.p.Pt() >=minPtMu 
	      && mu->tracker.p.Pt() <=maxPtMu;
	    bool ptrange_tev = mu->tev_1st.p.Pt() >=minPtMu 
	      && mu->tev_1st.p.Pt() <=maxPtMu;

	    // only one muon with tracker pt above OneMuPtTrackCut
	    if(!OneMuon)continue;

	    // bin-2 corresponds to delta-R = 0.3
	    // Sum-Pt isolation: SumPtCut threshold
	    if(mu->SumPtIso[2] > SumPtCut )continue;
	    
	    bool JetActivity = (NjetAboveThresh(EtJetCut, Delta_Phi, 
						mu, ev) > 0);

	    // jet veto: no jet above EtJetCut
	    if(JetActivity)continue;

	    // quality track: tracker pt > 60 GeV + chi^2/ndof < 10 and |eta| < 1.8
	    if(mu->tracker.p.Pt() < PtTrackCut)continue;

	    bool check_glb = ((mu->global.chi2 / mu->global.ndof)<Chi2Cut)
	      && TMath::Abs(mu->global.p.Eta()) < Muon_Eta_Cut;

	    bool check_trk = ((mu->tracker.chi2 / mu->tracker.ndof) 
			      < Chi2Cut)
	      && TMath::Abs(mu->tracker.p.Eta()) < Muon_Eta_Cut;

	    bool	check_tev = ((mu->tev_1st.chi2 / mu->tev_1st.ndof) 
				     < Chi2Cut)
	      && TMath::Abs(mu->tev_1st.p.Eta()) < Muon_Eta_Cut;

	    fillHistos(mu, weight, survived_cut, 
			ptrange_glb && check_glb,ptrange_trk && check_trk,
			ptrange_tev && check_tev);

	  } // loop over muons

	for(int j = 0; j != Num_trkAlgos; ++j)
	  if(survived_cut[j])
	    ++Num_surv_cut[j];
	
      } // event loop

    Nexp_evt += nevents * weight; // total # of events (before any cuts)

    for(int j = 0; j != Num_trkAlgos; ++j)
      Nexp_evt_cut[j] += Num_surv_cut[j] * weight;
    
    delete ev;
    
  } //end loop over files

  printSummary();
  saveHistos(fout, dir);
}

void printSummary()
{
  cout << " Total # of expected events = " << Nexp_evt << endl;

  for(int j = 0; j != Num_trkAlgos; ++j)
    {
      cout << "\n Algorithm:" << algo_desc[j] << endl;
      cout << " Expected evts = " << Nexp_evt_cut[j];
      float eff, deff;
      getEff(eff, deff, Nexp_evt_cut[j], Nexp_evt);
      cout << ", Relative eff = "<<eff*100 << " +- " << deff*100 << "%";
      cout << endl;
    } // loop over tracking algorithms 
}

void saveHistos(TFile * fout, string dir)
{
  fout->cd(); fout->mkdir(dir.c_str()); fout->cd(dir.c_str());

  for(int j = 0; j != Num_trkAlgos; ++j)
    {
      hPTplus[j]->Write();
      hPTminus[j]->Write();
    }
  return;
}


void initHistos()
{
  hPTplus[0]= new TH1F("hPTglb_plus","(+) Global Muon Pt qual",
			  nBinPtMu,minPtMu,maxPtMu);
  hPTplus[1]= new TH1F("hPTtrk_plus","(+) Tracker Muon Pt qual",
			  nBinPtMu,minPtMu,maxPtMu);
  hPTplus[2] = new TH1F("hPTtev_plus","(+) TeV-1st Muon Pt qual",
			   nBinPtMu,minPtMu,maxPtMu);

  hPTminus[0]= new TH1F("hPTglb_minus","(-) Global Muon Pt qual",
			  nBinPtMu,minPtMu,maxPtMu);
  hPTminus[1]= new TH1F("hPTtrk_minus","(-) Tracker Muon Pt qual",
			  nBinPtMu,minPtMu,maxPtMu);
  hPTminus[2] = new TH1F("hPTtev_minus","(-) TeV-1st Muon Pt qual",
			   nBinPtMu,minPtMu,maxPtMu);
}

void fillHistos(const wprime::Muon * mu, float weight, 
		       bool * survived_cut, bool glb_flag, bool trk_flag, 
		       bool tev_flag)
{
  if(glb_flag)
    {
      if(mu->global.q > 0)
	hPTplus[0]->Fill(mu->global.p.Pt(), weight);
      else
	hPTminus[0]->Fill(mu->global.p.Pt(), weight);
      survived_cut[0] = true;
    }
  if(trk_flag)
    {
      if(mu->tracker.q > 0)
	hPTplus[1]->Fill(mu->tracker.p.Pt(), weight);
      else
	hPTminus[1]->Fill(mu->tracker.p.Pt(), weight);
      survived_cut[1] = true;
    }
  if(tev_flag)
    {
      if(mu->tev_1st.q > 0)
	hPTplus[2]->Fill(mu->tev_1st.p.Pt(), weight);
      else
	hPTminus[2]->Fill(mu->tev_1st.p.Pt(), weight);
      survived_cut[2] = true;
    }
}
