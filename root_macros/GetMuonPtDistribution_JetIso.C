// ROOT stuff
#include <TROOT.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

// wprime stuff
#include "UserCode/CMGWPrimeGroup/interface/wprimeEvent.h"
#include "UserCode/CMGWPrimeGroup/root_macros/input_file.h"
#include "UserCode/CMGWPrimeGroup/root_macros/util.h"

// std stuff
#include <iostream>


using std::cout; using std::endl; using std::vector;

// jet distributions (# of jets and Et)
static const unsigned nBinNJets = 50;
static const float minNJets = -0.5;
static const float maxNJets = 49.5;
static const unsigned nBinEtJets = 150;
static const float minEtJets = 0;
static const float maxEtJets = 300;
// jet-activity veto (# of jets and Et)
static const unsigned nBinEtJets_veto = 10;
static const float minEtJets_veto = 0;
static const float maxEtJets_veto = 200;
static const unsigned nBinNJets_veto = 10;
static const float minNJets_veto = -0.5;
static const float maxNJets_veto = 9.5;

static const unsigned nBinSumPt = 120;
static const float minSumPt = 0;
static const float maxSumPt = 600;


void GetMuonPtDistribution_JetIso(const vector<wprime::InputFile>& files, 
				  TFile *fout, string dir)
{
  
  TH1F* hNMu = new TH1F("NMu","Nb. muons in the event",10,0.,10.);
  TH1F* hPtMaxMu= new TH1F("ptMaxMu","Pt mu",nBinPtMu,minPtMu,maxPtMu);
  
  TH1F* hPtMaxMuTrackVeto= new TH1F("ptMaxMuTrackVeto","Pt mu track Veto",nBinPtMu,minPtMu,maxPtMu); 
  
  
  TH1F* hPtMaxMuJetVeto= new TH1F("ptMaxMuJetVeto","Pt mu JetVeto",nBinPtMu,minPtMu,maxPtMu); 
  
  TH1F* hPtMaxMuTrackVetoJetVeto= new TH1F("ptMaxMuTrackVetoJetVeto","Pt mu TrackVeto and JetVeto ",nBinPtMu,minPtMu,maxPtMu); 
  
  hPtMaxMu->SetLineColor(2); hPtMaxMuTrackVeto->SetLineColor(4);
  hPtMaxMuJetVeto->SetLineColor(6); 
  //hPtMaxMuTrackVetoJetVeto->SetLineColor(1);

  int Nfiles = files.size();
  
  for(int tr = 0; tr != Nfiles; ++tr){
    if(!files[tr].tree)
      continue;
    
    wprime::Event * ev = new wprime::Event();
    files[tr].tree->SetBranchAddress("wp", &ev);
    
    int nevents = files[tr].tree->GetEntries();
    for(int i = 0; i != nevents; ++i)
      { // event loop
	files[tr].tree->GetEntry(i);
	
	int nmuons = ev->mu->GetLast() + 1;
	hNMu->Fill(nmuons,files[tr].weight);
	float ptMaxMu = 0.0;
	float ptMaxMuTrackVeto  = 0.0;
	for(int j = 0; j != nmuons; ++j)
	  { // loop over muons
	    wprime::Muon * mu = (wprime::Muon *) ev->mu->At(j);
	    float ptmu =  mu->global.p.Pt();
	    if(ptmu > ptMaxMu) {ptMaxMu = ptmu;}
	    
	    if( mu->tracker.p.Pt() > PtTrackCut )
	      {
		if(ptmu > ptMaxMuTrackVeto) {ptMaxMuTrackVeto = ptmu;}
	      }
	  } // loop over muons

	if(ptMaxMu > 0 ){hPtMaxMu->Fill(ptMaxMu,files[tr].weight);}
	if(ptMaxMuTrackVeto > 0 )
	  {hPtMaxMuTrackVeto->Fill(ptMaxMuTrackVeto,files[tr].weight);}
	
	int njet = ev->jet->GetLast() + 1;
	int nJetEt = 0;
	for(int i =0 ;i != njet; ++i){
	  TLorentzVector*	 ajet = ( TLorentzVector*) ev->jet->At(i);
	  if ( ajet->Pt() > EtJetCut ) {nJetEt++;}
	}
	if(nJetEt == 0){
	  if(ptMaxMu > 0 ){hPtMaxMuJetVeto->Fill(ptMaxMu,files[tr].weight);}
	  if(ptMaxMuTrackVeto > 0) {hPtMaxMuTrackVetoJetVeto->Fill(ptMaxMuTrackVeto, files[tr].weight);}
	}

      } // event loop
    
    delete ev;
    
  } //end loop over files

  fout->cd(); fout->mkdir(dir.c_str()); fout->cd(dir.c_str());
  hPtMaxMu->Write(); hPtMaxMuTrackVetoJetVeto->Write(); 
  hPtMaxMuJetVeto->Write();  hPtMaxMuTrackVeto->Write();

  return;

}
