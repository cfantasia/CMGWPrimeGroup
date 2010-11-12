// ROOT stuff
#include <TROOT.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>

// wprime stuff
#include "UserCode/CMGWPrimeGroup/interface/wprimeEvent.h"

// std stuff
#include <iostream>

using std::cout; using std::endl;

int read_wprime()
{
  TFile *file0 = new TFile("wprime_1TeV.root");
  if(file0->IsZombie())
    {
      cerr << " *** Input file does not exist! " << endl;
      return -1;
    }
  TTree * wp = (TTree *) file0->Get("StdMu/wprime");
  
  if(!wp->GetBranch("wp"))
    {
      cerr << " *** Can't find wp branch! " << endl;
      return -2;
    }

  wprime::Event * ev = new wprime::Event();
  wp->SetBranchAddress("wp", &ev);


  int nevents = wp->GetEntries();
  cout << " Opened file with " << nevents << " entries" << endl;
  //  for(int i = 0; i != nevents; ++i)
  for(int i = 0; i != 10; ++i)
    { // event loop
      wp->GetEntry(i);
      cout << " ==============================================================\n";
      cout << " Event = " << i << " Evt # " << ev->evt_no << " run # " 
	   << ev->run_no << endl;
      //      cout << " HLT_Mu5 = " << bool(ev->HLT_Mu5) << endl;
      int nmuons = ev->mu->GetLast() + 1;
	cout << " # of muons = " << nmuons << endl;

      for(int j = 0; j != nmuons; ++j)
	{ // loop over muons
	  wprime::Muon * mu = (wprime::Muon *) ev->mu->At(j);
	  cout << " Muon = " << j << ", # of std-mu hits = " << mu->Nmu_hits
	       << endl;
	  cout << " tracker pt = " << mu->tracker.p.Pt()
	       << ", global pt = " << mu->global.p.Pt()
	       << ", tev 1st pt = " << mu->tpfms.p.Pt() << endl;

	} // loop over muons

      

    } // event loop

  return 0;
}
