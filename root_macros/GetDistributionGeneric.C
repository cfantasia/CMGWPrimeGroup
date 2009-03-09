// ROOT stuff
#include <TROOT.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>

#include <TH1F.h>
// wprime stuff
#include "UserCode/CMGWPrimeGroup/interface/wprimeEvent.h"

// std stuff
#include <iostream>


using std::cout; using std::endl;


  static const unsigned  nBinPtMu = 1500;
  static const float minPtMu = 0;
  static const float  maxPtMu = 3000;

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

  // cone size and SumPt for isolation
  static const unsigned nBinCone = wprime::N_CONESIZE;
  static const float minCone = wprime::MIN_CONE;
  static const float maxCone = wprime::MAX_CONE;
  static const unsigned nBinSumPt = 120;
  static const float minSumPt = 0;
  static const float maxSumPt = 600;





float PtJetCut = 80;
float PtTrackCut = 60 ;

int GetDistributionGeneric(const int Nfiles, const string* fileName, const float* weight , TFile *fout, string dir)
{

TH1F* hNMu = new TH1F("NMu","Nb. muons in the event",10,0.,10.);
TH1F* hPtMaxMu= new TH1F("ptMaxMu","Pt mu",nBinPtMu,minPtMu,maxPtMu);

TH1F* hPtMaxMuTrackVeto= new TH1F("ptMaxMuTrackVeto","Pt mu track Veto",nBinPtMu,minPtMu,maxPtMu); 


TH1F* hPtMaxMuJetVeto= new TH1F("ptMaxMuJetVeto","Pt mu JetVeto",nBinPtMu,minPtMu,maxPtMu); 

TH1F* hPtMaxMuTrackVetoJetVeto= new TH1F("ptMaxMuTrackVetoJetVeto","Pt mu TrackVeto and JetVeto ",nBinPtMu,minPtMu,maxPtMu); 

hPtMaxMu->SetLineColor(2); hPtMaxMuTrackVeto->SetLineColor(4);
 hPtMaxMuJetVeto->SetLineColor(6); //hPtMaxMuTrackVetoJetVeto->SetLineColor(1);


 
 TFile* file[Nfiles];
 for(int i =0;i<Nfiles;i++){file[i]= 0;}
 for(int i =0;i<Nfiles;i++){file[i]= new TFile(fileName[i].c_str());}
 TTree* tree[Nfiles];
 for(int i =0;i<Nfiles;i++){
   if(!(file[i]->IsOpen())){cout<<"Missing file: "<<fileName[i]<<" !!!!!! "<<endl; continue;}
   tree[i]=(TTree*)file[i]->Get("StdMu/wprime");
  if(!tree[i]->GetBranch("wp"))
    {
      cerr << " *** Can't find wp branch! in file:"<< fileName[i]<< endl;
      return -1;
    }

 }


 for(int tr =0;tr<Nfiles;tr++){
   if(!(file[tr]->IsOpen())){continue;}
  wprime::Event * ev = new wprime::Event();
  tree[tr]->SetBranchAddress("wp", &ev);


  int nevents = tree[tr]->GetEntries();
  cout << " Opened file with " << nevents << " entries" << endl;
  for(int i = 0; i != nevents; ++i)
    { // event loop
      tree[tr]->GetEntry(i);

      int nmuons = ev->mu->GetLast() + 1;
      hNMu->Fill(nmuons,weight[tr]);
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
//          cout << " Muon = " << j << ", # of std-mu hits = " << mu->Nmu_hits
//	       << endl;
//	  cout << " tracker pt = " << mu->tracker.p.Pt()
//	       << ", global pt = " << mu->global.p.Pt()
//	       << ", tev 1st pt = " << mu->tev_1st.p.Pt() << endl;

	} // loop over muons
      if(ptMaxMu > 0 ){hPtMaxMu->Fill(ptMaxMu,weight[tr]);}
      if(ptMaxMuTrackVeto > 0 ){hPtMaxMuTrackVeto->Fill(ptMaxMuTrackVeto,weight[tr]);}

      int njet = ev->jet->GetLast() + 1;
      int nJetEt = 0;
      for(int i =0 ;i != njet; ++i){
	TLorentzVector*	 ajet = ( TLorentzVector*) ev->jet->At(i);
	if ( ajet->Pt() > PtJetCut ) {nJetEt++;}
      }
      if(nJetEt == 0){
	if(ptMaxMu > 0 ){hPtMaxMuJetVeto->Fill(ptMaxMu,weight[tr]);}
	if(ptMaxMuTrackVeto > 0) {hPtMaxMuTrackVetoJetVeto->Fill(ptMaxMuTrackVeto,weight[tr]);}
      }
    } // event loop

 }//end loop over files

fout->cd(); fout->mkdir(dir.c_str()); fout->cd(dir.c_str());
hPtMaxMu->Write(); hPtMaxMuTrackVetoJetVeto->Write(); 
hPtMaxMuJetVeto->Write();  hPtMaxMuTrackVeto->Write();

  return 0;
}
