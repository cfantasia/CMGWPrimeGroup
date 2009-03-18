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

// std stuff
#include <iostream>


using std::cout; using std::endl; using std::vector;


static const unsigned  nBinPtMu = 380;
static const float minPtMu = 100;
static const float  maxPtMu = 2000;

// cone size and SumPt for isolation
static const unsigned nBinCone = wprime::N_CONESIZE;
static const float minCone = wprime::MIN_CONE;
static const float maxCone = wprime::MAX_CONE;

const float EtJetCut = 80;
const float SumPtCut = 3;
const float PtTrackCut = 60 ;
const float OneMuPtTrackCut = 20;

const int Num_histo_sets = 5; // one new set after each cut
const int Num_trkAlgos = 3; // global, tracker, tev_1st

TH1F * hPT[Num_histo_sets][Num_trkAlgos] = {0};

string cut_desc[Num_histo_sets]= {" Pt within range:", " 1 muon only:",
				  " Isolation:", " Jet Veto:", " Quality:"};

string algo_desc[Num_trkAlgos] = {" global", " tracker", " TeV-1st"};

float Nexp_evt = 0;
float Nexp_evt_cut[Num_histo_sets][Num_trkAlgos] = {0};

void initHistos();
void saveHistos(TFile * fout, string dir);
void printSummary();
void getEff(float & eff, float & deff, float Num, float Denom);

void GetMuonPtDistribution(const vector<wprime::InputFile>& files, 
			   TFile *fout, string dir)
{

  initHistos();
  int Nfiles = files.size();

  // counter initialization
  Nexp_evt = 0;
  for(int i = 0; i != Num_histo_sets; ++i)
    for(int j = 0; j != Num_trkAlgos; ++j)
      Nexp_evt_cut[i][j] = 0;

  for(int tr = 0; tr != Nfiles; ++tr){
    if(!files[tr].tree)
      continue;
    
    cout << " Processing sample " << files[tr].description << endl;

    wprime::Event * ev = new wprime::Event();
    files[tr].tree->SetBranchAddress("wp", &ev);
    
    int nevents = files[tr].tree->GetEntries();
    float weight = files[tr].weight;

    int Num_surv_cut[Num_histo_sets][Num_trkAlgos] = {0};
    for(int i = 0; i != nevents; ++i)
      { // event loop
	files[tr].tree->GetEntry(i);
	
	bool survived_cut[Num_histo_sets][Num_trkAlgos] = {false};
	
       	int nmuons = ev->mu->GetLast() + 1;

	bool HLT = ev->HLT_Mu5;

	bool OneMuon = (NmuAboveThresh(OneMuPtTrackCut, ev) == 1);

	bool JetActivity = (NjetAboveThresh(EtJetCut, ev) > 0);

       	for(int j = 0; j != nmuons; ++j)
	  { // loop over muons
	    wprime::Muon * mu = (wprime::Muon *) ev->mu->At(j);

	    int index = 0;
	
	    bool eta_LT_21 = TMath::Abs(mu->global.p.Eta()) < 2.1;

	    hPT[index][0]->Fill(mu->global.p.Pt(), weight);
	    hPT[index][1]->Fill(mu->tracker.p.Pt(), weight);
	    hPT[index][2]->Fill(mu->tev_1st.p.Pt(), weight);
	    
	    bool ptrange_glb = mu->global.p.Pt() >=minPtMu 
	      && mu->global.p.Pt() <=maxPtMu && eta_LT_21;
	    if(ptrange_glb)
	      survived_cut[index][0] = true;

	    bool ptrange_trk = mu->tracker.p.Pt() >=minPtMu 
	      && mu->tracker.p.Pt() <=maxPtMu && eta_LT_21;
	    if(ptrange_trk)
	      survived_cut[index][1] = true;

	    bool ptrange_tev = mu->tev_1st.p.Pt() >=minPtMu 
	      && mu->tev_1st.p.Pt() <=maxPtMu && eta_LT_21;

	    if(ptrange_tev)
	      survived_cut[index][2] = true;

	    ++index;

	    /*
	    // trigger decision
	    if(!HLT)continue;
	    
	    hPT[index][0]->Fill(mu->global.p.Pt(), weight);
	    hPT[index][1]->Fill(mu->tracker.p.Pt(), weight);
	    hPT[index][2]->Fill(mu->tev_1st.p.Pt(), weight);

	    survived_cut[index][0] = survived_cut[index][1] = survived_cut[index][2] = true;
	    ++index;
	    */

	    // only one muon with tracker pt above OneMuPtTrackCut
	    if(!OneMuon)continue;

	    hPT[index][0]->Fill(mu->global.p.Pt(), weight);
	    hPT[index][1]->Fill(mu->tracker.p.Pt(), weight);
	    hPT[index][2]->Fill(mu->tev_1st.p.Pt(), weight);
	    // only keep muons with pt within range
	    if(ptrange_glb)
	      survived_cut[index][0] = true;
	    if(ptrange_trk)
	      survived_cut[index][1] = true;
	    if(ptrange_tev)
	      survived_cut[index][2] = true;
	    ++index;

	    // bin-2 corresponds to delta-R = 0.3
	    // Sum-Pt isolation: SumPtCut threshold
	    if(mu->SumPtIso[2] > SumPtCut )continue;

	    hPT[index][0]->Fill(mu->global.p.Pt(), weight);
	    hPT[index][1]->Fill(mu->tracker.p.Pt(), weight);
	    hPT[index][2]->Fill(mu->tev_1st.p.Pt(), weight);

	    // only keep muons with pt within range
	    if(ptrange_glb)
	      survived_cut[index][0] = true;
	    if(ptrange_trk)
	      survived_cut[index][1] = true;
	    if(ptrange_tev)
	      survived_cut[index][2] = true;
	    ++index;

	    // jet veto: no jet above EtJetCut
	    if(JetActivity)continue;

	    hPT[index][0]->Fill(mu->global.p.Pt(), weight);
	    hPT[index][1]->Fill(mu->tracker.p.Pt(), weight);
	    hPT[index][2]->Fill(mu->tev_1st.p.Pt(), weight);

	    // only keep muons with pt within range
	    if(ptrange_glb)
	      survived_cut[index][0] = true;
	    if(ptrange_trk)
	      survived_cut[index][1] = true;
	    if(ptrange_tev)
	      survived_cut[index][2] = true;

	    ++index;

	    // quality track: tracker pt > 60 GeV + chi^2/ndof < 10
	    if(mu->tracker.p.Pt() < PtTrackCut)continue;

	    bool check = ((mu->global.chi2 / mu->global.ndof) < 10);
	    if(check && ptrange_glb)
	      {
		hPT[index][0]->Fill(mu->global.p.Pt(), weight);
		survived_cut[index][0] = true;
	      }
	    check = ((mu->tracker.chi2 / mu->tracker.ndof) < 10);
	    if(check && ptrange_trk)
	      {
		hPT[index][1]->Fill(mu->tracker.p.Pt(), weight);
		survived_cut[index][1] = true;
	      }
	    check = ((mu->tev_1st.chi2 / mu->tev_1st.ndof) < 10);
	    if(check && ptrange_tev)
	      {
		hPT[index][2]->Fill(mu->tev_1st.p.Pt(), weight);
		survived_cut[index][2] = true;
	      }
	    

	  } // loop over muons

	for(int i = 0; i != Num_histo_sets; ++i)
	  for(int j = 0; j != Num_trkAlgos; ++j)
	    if(survived_cut[i][j])
	      ++Num_surv_cut[i][j];
	
      } // event loop

    Nexp_evt += nevents * weight; // total # of events (before any cuts)

    for(int i = 0; i != Num_histo_sets; ++i)
      for(int j = 0; j != Num_trkAlgos; ++j)
	Nexp_evt_cut[i][j] += Num_surv_cut[i][j] * weight;
    
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
      for(int i = 0; i != Num_histo_sets; ++i)
	{
	  cout << cut_desc[i] << " expected evts = " << Nexp_evt_cut[i][j];
	  float eff, deff;
	  if(i == 0)
	    getEff(eff, deff, Nexp_evt_cut[i][j], Nexp_evt);
	  else
	    getEff(eff, deff, Nexp_evt_cut[i][j], Nexp_evt_cut[i-1][j]);
	  cout << ", Relative eff = "<<eff*100 << " +- " << deff*100 << "%";
	  getEff(eff, deff, Nexp_evt_cut[i][j], Nexp_evt);
	  cout << ", Absolute eff = "<< eff*100 << " +- " << deff*100 << "%"
	       << endl;
	} // loop over different cuts
    } // loop over tracking algorithms 

}

void getEff(float & eff, float & deff, float Num, float Denom)
{
  eff = Num/Denom;
  deff = TMath::Sqrt(eff * (1-eff)/Denom);
}

void saveHistos(TFile * fout, string dir)
{
  fout->cd(); fout->mkdir(dir.c_str()); fout->cd(dir.c_str());

  for(int i = 0; i != Num_histo_sets; ++i)
    for(int j = 0; j != Num_trkAlgos; ++j)
      hPT[i][j]->Write();
  
  return;
}


void initHistos()
{
  int index = 0;
  hPT[index][0] = new TH1F("hPTglb_all","Global Muon Pt",nBinPtMu,minPtMu,maxPtMu);
  hPT[index][1] = new TH1F("hPTtrk_all","Tracker Muon Pt",nBinPtMu,minPtMu,maxPtMu);
  hPT[index][2] = new TH1F("hPTtev_all","TeV-1st Muon Pt",nBinPtMu,minPtMu,maxPtMu);
  ++index;
  
  /*
    hPT[index][0]= new TH1F("hPTglb_trig","Global Muon Pt with HLT",nBinPtMu,minPtMu,maxPtMu);
    hPT[index][0]= new TH1F("hPTtrk_trig","Tracker Muon Pt with HLT",nBinPtMu,minPtMu,maxPtMu);
    hPT[index][0]= new TH1F("hPTtev_trig","TeV-1st Muon Pt with HLT",nBinPtMu,minPtMu,maxPtMu);
    ++index;
  */

  hPT[index][0]= new TH1F("hPTglb_1mu","Global Muon Pt 1 muon",nBinPtMu,minPtMu,maxPtMu);
  hPT[index][1]= new TH1F("hPTtrk_1mu","Tracker Muon Pt 1 muon",nBinPtMu,minPtMu,maxPtMu);
  hPT[index][2]= new TH1F("hPTtev_1mu","TeV-1st Muon Pt 1 muon",nBinPtMu,minPtMu,maxPtMu);
  ++index;

  hPT[index][0]= new TH1F("hPTglb_iso","Global Muon Pt iso",nBinPtMu,minPtMu,maxPtMu);
  hPT[index][1]= new TH1F("hPTtrk_iso","Tracker Muon Pt iso",nBinPtMu,minPtMu,maxPtMu);
  hPT[index][2]= new TH1F("hPTtev_iso","TeV-1st Muon Pt iso",nBinPtMu,minPtMu,maxPtMu);
  ++index;
  
  hPT[index][0]= new TH1F("hPTglb_jveto","Global Muon Pt jet veto",nBinPtMu,minPtMu,maxPtMu);
  hPT[index][1]= new TH1F("hPTtrk_jveto","Tracker Muon Pt jet veto",nBinPtMu,minPtMu,maxPtMu);
  hPT[index][2]= new TH1F("hPTtev_jveto","TeV-1st Muon Pt jet veto",nBinPtMu,minPtMu,maxPtMu);
  ++index;
  
  hPT[index][0]= new TH1F("hPTglb_qual","Global Muon Pt qual",nBinPtMu,minPtMu,maxPtMu);
  hPT[index][1]= new TH1F("hPTtrk_qual","Tracker Muon Pt qual",nBinPtMu,minPtMu,maxPtMu);
  hPT[index][2] = new TH1F("hPTtev_qual","TeV-1st Muon Pt qual",nBinPtMu,minPtMu,maxPtMu);
  ++index;
}
