// ROOT49 stuff
#include <TROOT.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1F.h>

#include <map>

// wprime stuff
#include "UserCode/CMGWPrimeGroup/interface/wprimeEvent.h"
#include "UserCode/CMGWPrimeGroup/root_macros/input_file.h"
#include "UserCode/CMGWPrimeGroup/root_macros/loadCutsAndThresholds.h"

using std::cout; using std::endl; using std::string; using std::map;

//methods in this macro:
void printHighPtMuon(const wprime::Event * ev, const wprime::Muon* mu,
		     const bool fill_entry[]);
void printSummary_MuonPt(ofstream & out, const string& dir, 
                         float Nexp_evt, 
			 const float Nexp_evt_cut[][Num_trkAlgos]);
void defineHistos_MuonPt();
void defineHistos_MuonEta();
void defineHistos_MuonPhi();
void defineHistos_MuonJetDPhi();
void defineHistos_MuonIso();
void defineHistos_TMass();
void defineHistos_MuonChargePt();
void defineHistos(int option);
void tabulateMe(int Num_surv_cut[][Num_trkAlgos], int& cut_index, 
		float weight,const wprime::Event* ev, const wprime::Muon* mu,
                const bool fill_entry[], int option, 
		bool accountme[][Num_trkAlgos]);
void fillHistos_MuonPt(int index, float weight, const TLorentzVector **p,
                       const bool fill_entry[]);
void fillHistos_MuonEta(int index, float weight, const TLorentzVector **p,
                       const bool fill_entry[]);
void fillHistos_MuonPhi(int index, float weight, const TLorentzVector **p,
                        const bool fill_entry[]);

void fillHistos_MuonJetDPhi(int index, float weight, const wprime::Event * ev,
                        const TLorentzVector **p, const bool fill_entry[]);
void fillHistos_TMass(int index, float weight, const wprime::Event * ev,
		      const TLorentzVector **p, const bool fill_entry[]);
void fillHistos_MuonIso(int index, float weight, const wprime::Muon* mu,
			const bool fill_entry[]);
void fillHistos_MuonChargePt(int index, float weight, const Int_t * Q, 
			     const TLorentzVector **p,const bool fill_entry[]);
void saveHistos_MainAnalysis();
void saveHistos_MuonChargePt();
void saveHistos(TFile * fout, string dir, int option);
void gatherFileBasicInfo(const wprime::InputFile& file,
                         wprime::Event*& ev, int & nevents, float & weight);

// Studies:
void GetMuonPtDistribution(const wprime::InputFile& file, 
			   float Nexp_evt_cut[][Num_trkAlgos],
                           float & Nexp_evt, int option,
			   bool highestPtMuonOnly);


//--------------------------------------------------------------------------
void printSummary_MuonPt(ofstream & out, const string& sample, 
                         float Nexp_evt, 
			 const float Nexp_evt_cut[][Num_trkAlgos])
//------------------------------------------------------------------------
{
    cout << "\n Sample: " << sample << endl;
    cout << " Total # of expected events = " << Nexp_evt << endl;

    for (int mual = 0; mual != Num_trkAlgos; ++mual){
      cout << " Algorithm: " << algo_desc_long[mual] << endl;
      for(int i = 0; i < Num_selection_cuts; ++i){

	//print in a txt file the final total number of events
	if(i == Num_selection_cuts - 1)
	  out << algo_desc_long[mual] << " " << sample << " " 
	      << Nexp_evt_cut[i][mual] << endl;
	
	cout << " Cut # " << i << ": " << cuts_desc_long[i] 
	     <<", expected # of evts = " << Nexp_evt_cut[i][mual];
	
	//calculate efficiencies
	float eff, deff;
	if(i == 0)
	  getEff(eff, deff, Nexp_evt_cut[i][mual], Nexp_evt);
	else
	  getEff(eff, deff, Nexp_evt_cut[i][mual], 
		 Nexp_evt_cut[i-1][mual]);
	cout << ", Relative eff = "<<eff*100 << " +- " << deff*100 << "%";
	getEff(eff, deff, Nexp_evt_cut[i][mual], Nexp_evt);
	cout << ", Absolute eff = "<< eff*100 << " +- " << deff*100 << "%"
	     << endl;
            
      } // loop over different cuts
    }//loop over different muon algos

}//------------- printSummary_MuonPt()


//--------------------------------------------------------------
void defineHistos_MuonPt()
{
//--------------------------------------------------------------
#if debugme
  cout << " define Muon pT histos" << endl;
#endif

  for(int cut = 0; cut != Num_selection_cuts; ++cut)
    for(int algo = 0; algo != Num_trkAlgos; ++algo)
      {
	string name = "hPT" + algo_desc_short[algo] + "_" + cuts_desc_short[cut];
	string title = algo_desc_long[algo] + " muon p_{T} with " + cuts_desc_long[cut];
	hPT[cut][algo] = new TH1F(name.c_str(), title.c_str(), nBinPtMu,minPtMu,maxPtMu);
      }
  
}//---------defineMuonPtHistos()



//--------------------------------------------------------------
void defineHistos_TMass()
{
//--------------------------------------------------------------
#if debugme
  cout << " define TMass histos" << endl;
#endif

  for(int cut = 0; cut != Num_selection_cuts; ++cut)
    for(int algo = 0; algo != Num_trkAlgos; ++algo)
      {
	string name = "hTM" + algo_desc_short[algo] + "_" + cuts_desc_short[cut];
	string title = algo_desc_long[algo] + " Transv. Mass with " + cuts_desc_long[cut];
	hTM[cut][algo] = new TH1F(name.c_str(), title.c_str(), nBinTmMu,minTmMu,maxTmMu);
      }
  
}//---------defineMuonPtHistos()



//--------------------------------------------------------------
void defineHistos_MuonEta()
{
//--------------------------------------------------------------
#if debugme
  cout << " define Muon ETA histos" << endl;
#endif

  for(int cut = 0; cut != Num_selection_cuts; ++cut)
    for(int algo = 0; algo != Num_trkAlgos; ++algo)
      {
	string name = "hETA" + algo_desc_short[algo] + "_" + cuts_desc_short[cut];
	string title = algo_desc_long[algo] + " muon #eta with " + cuts_desc_long[cut];
	hETA[cut][algo] = new TH1F(name.c_str(), title.c_str(), nBinEtaMu,minEtaMu,maxEtaMu);
      }
  
}//---------defineMuonEtaHistos()



//--------------------------------------------------------------
void defineHistos_MuonPhi()
{
//--------------------------------------------------------------
#if debugme
  cout << " define Muon PHI histos" << endl;
#endif

  for(int cut = 0; cut != Num_selection_cuts; ++cut)
    for(int algo = 0; algo != Num_trkAlgos; ++algo)
      {
	string name = "hPHI" + algo_desc_short[algo] + "_" + cuts_desc_short[cut];
	string title = algo_desc_long[algo] + " muon #phi with " + cuts_desc_long[cut];
	hPHI[cut][algo] = new TH1F(name.c_str(), title.c_str(), nBinPhiMu,minPhiMu,maxPhiMu);
      }
  
}//---------defineMuonPhiHistos()



//--------------------------------------------------------------
void defineHistos_MuonJetDPhi()
{
//--------------------------------------------------------------
#if debugme
  cout << " define Muon-Jet DPhi histos" << endl;
#endif

  for(int cut = 0; cut != Num_selection_cuts; ++cut)
    for(int algo = 0; algo != Num_trkAlgos; ++algo)
      {
	string name = "hMJDPHI" + algo_desc_short[algo] + "_" + cuts_desc_short[cut];
	string title = algo_desc_long[algo] + " muon-jet #Delta#phi with " + cuts_desc_long[cut];
	hMJDPHI[cut][algo] = new TH1F(name.c_str(), title.c_str(), nBinDPhiMu,minDPhiMu,maxDPhiMu);
      }
  
}//---------defineMuonJetDPhiHistos()



//--------------------------------------------------------------
void defineHistos_MuonIso()
{
//--------------------------------------------------------------
#if debugme
  cout << " define Muon ISO histos" << endl;
#endif

  for(int cut = 0; cut != Num_selection_cuts; ++cut)
    for(int algo = 0; algo != Num_trkAlgos; ++algo)
      {
	string name = "hISO" + algo_desc_short[algo] + "_" + cuts_desc_short[cut];
	string title = algo_desc_long[algo] + " muon isol with " + cuts_desc_long[cut];
	hISO[cut][algo] = new TH1F(name.c_str(), title.c_str(), nBinIsoMu,minIsoMu,maxIsoMu);
      }
  
}//---------defineMuonPhiHistos()


//--------------------------------------------------------------
void defineHistos_MuonChargePt()
{
//--------------------------------------------------------------
#if debugme
  cout << " define histos MuonChargePt " << endl;
#endif

  for(int algo = 0; algo != Num_trkAlgos; ++algo)
    {
      string name = "hPT" + algo_desc_short[algo] + "_" + "plus";
      string title = " (+) " + algo_desc_long[algo] + " muon p_{T} qual ";
      hPTplus[algo] = new TH1F(name.c_str(), title.c_str(), nBinPtMu,minPtMu,maxPtMu);
      name = "hPT" + algo_desc_short[algo] + "_" + "minus";
      title = " (-) " + algo_desc_long[algo] + " muon p_{T} qual ";
      hPTminus[algo] = new TH1F(name.c_str(), title.c_str(), nBinPtMu,minPtMu,maxPtMu);
    }
    
  return;

}//------defineHistos_MuonChargePt



// Define histograms depending on the type of study
//--------------------------------------------------------------
void defineHistos(int option)
{
//--------------------------------------------------------------

#if debugme
    cout << " define histos with option = " << option<<endl;
#endif

    if(option == 1)
      {
        defineHistos_MuonPt();
        defineHistos_MuonEta();
        defineHistos_MuonPhi();
        defineHistos_MuonJetDPhi();
        defineHistos_MuonIso();
        defineHistos_TMass();
      }
    else if(option == 2) 
      defineHistos_MuonChargePt();
    else 
    {
      cout << " WARNING!!! No histos were defined for option = " << option  
	   << endl;
      return;
    }

}//----defineHistos()


//tabulate results after the cut has been passed
//TODO: Improve flexibility
//-----------------------------------------------------------
void tabulateMe(int Num_surv_cut[][Num_trkAlgos], int& cut_index, 
                float weight, const wprime::Event * ev,
                const wprime::Muon* mu,
                const bool fill_entry[], int option,
                bool accountme[][Num_trkAlgos])
{
//-----------------------------------------------------------
#if debugme
  cout << " Tabulating results for cut_index = " << cut_index << endl;
#endif

  //if the accountme switch is on,
  //increase the number of events passing the cuts
  //and turn the switch off so we don't count more than once per event
  for(int j = 0; j != Num_trkAlgos; ++j)
    if (accountme[cut_index][j] && fill_entry[j]) 
      {
	++(Num_surv_cut[cut_index][j]);
	accountme[cut_index][j] = false;
      }
  
  const TLorentzVector * P[Num_trkAlgos] = {
    &(mu->global.p), &(mu->tracker.p), &(mu->tpfms.p), &(mu->cocktail.p),
    &(mu->picky.p),  &(mu->tmr.p)};

  //fill the histograms
  if(option == 1) {
    fillHistos_MuonPt(cut_index,weight,P, fill_entry);
    fillHistos_MuonEta(cut_index,weight,P, fill_entry);
    fillHistos_MuonPhi(cut_index,weight,P, fill_entry);
    fillHistos_MuonJetDPhi(cut_index,weight, ev, P, fill_entry);
    fillHistos_MuonIso(cut_index,weight, mu, fill_entry);
    fillHistos_TMass(cut_index,weight, ev, P, fill_entry);
  }
  else if(option == 2) 
    {
      const Int_t charge[Num_trkAlgos] = {
	mu->global.q, mu->tracker.q, mu->tpfms.q, mu->cocktail.q,
	mu->picky.q, mu->tmr.q};
      
      fillHistos_MuonChargePt(cut_index, weight, charge, P, fill_entry);
    }
  else 
    cout << " WARNING!! I can't tabulate anything for option = " << option 
	 << endl;
  
#if debugme
  cout << " cut_index is now = " << cut_index << endl;
#endif

}//-----tabulateMe


//Fill Histograms
//-----------------------------------------------------------
void fillHistos_MuonPt(int index, float weight, const TLorentzVector **p,
                       const bool fill_entry[])
{
//-----------------------------------------------------------
  for(int algo = 0; algo != Num_trkAlgos; ++algo)
    {
      if(fill_entry[algo])
	hPT[index][algo]->Fill( (p[algo])->Pt(), weight);
    }
   
}//fillHistos_MuonPt


//-----------------------------------------------------------
void fillHistos_MuonEta(int index, float weight, const TLorentzVector **p,
                       const bool fill_entry[])
{
//-----------------------------------------------------------
  for(int algo = 0; algo != Num_trkAlgos; ++algo)
    {
      if(fill_entry[algo])
	hETA[index][algo]->Fill( (p[algo])->Eta(), weight);
    }
}//fillHistos_MuonEta

//-----------------------------------------------------------
void fillHistos_MuonPhi(int index, float weight, const TLorentzVector **p,
                       const bool fill_entry[])
{
//-----------------------------------------------------------
  for(int algo = 0; algo != Num_trkAlgos; ++algo)
    {
      if(fill_entry[algo])
	hPHI[index][algo]->Fill( (p[algo])->Phi(), weight);
    }

}//fillHistos_MuonPhi

//-----------------------------------------------------------
void fillHistos_MuonJetDPhi(int index, float weight, const wprime::Event * ev,
			    const TLorentzVector **p, const bool fill_entry[])
{
//-----------------------------------------------------------
  for(int algo = 0; algo != Num_trkAlgos; ++algo)
    {
      if(fill_entry[algo])
	hMJDPHI[index][algo]->Fill( XJetDPhi(*(p[algo]), ev), weight);
    }
  
}//fillHistos_MuonJetDPhi

//-----------------------------------------------------------
void fillHistos_TMass(int index, float weight, const wprime::Event * ev,
		      const TLorentzVector **p, const bool fill_entry[])
{
//-----------------------------------------------------------

// NB: this is most likely correct only for global muons; 
// corrections to (pf)MET are needed for the rest of the tracking algorithms
  for(int algo = 0; algo != Num_trkAlgos; ++algo)
    {
      if(fill_entry[algo])
	hTM[index][algo]->Fill(TMass( *(p[algo]), ev->pfmet), weight);
    }

}//fillHistos_TMass


//-----------------------------------------------------------
void fillHistos_MuonIso(int index, float weight, const wprime::Muon* mu,
                       const bool fill_entry[])
{
//-----------------------------------------------------------

  if(!mu) {
#if debugme
    cout << " No muon found in the event " << endl;
#endif
    return;
  }

  for(int algo = 0; algo != Num_trkAlgos; ++algo){
    if(fill_entry[algo])
      hISO[index][algo]->Fill(CombRelIsolation(mu,deltaRIsoIndex),weight);
  }
  
}//fillHistos_MuonIso

//Fill Histograms
//-----------------------------------------------------------
void fillHistos_MuonChargePt(int index, float weight, const Int_t * Q, 
			     const TLorentzVector **p, const bool fill_entry[])
{
//-----------------------------------------------------------

  //only fill the required histograms
  if (index < (Num_selection_cuts - Num_histo_sets_chargePt)) return;

  for(int algo = 0; algo != Num_trkAlgos; ++algo)
    {
      if(fill_entry[algo])
	{
	  if(Q[algo] > 0)
	    hPTplus[algo]->Fill( (p[algo])->Pt(), weight);
	  else
	    hPTminus[algo]->Fill( (p[algo])->Pt(), weight);
	}
    }
    
}//---------------fillHistos_MuonChargePt()



//------------------------------------------------------------------------
void saveHistos_MainAnalysis()
{
//------------------------------------------------------------------------
#if debugme
  cout << " Saving MuonPt histos " << endl;
#endif
  for(int i = 0; i != Num_selection_cuts; ++i)
      for(int j = 0; j != Num_trkAlgos; ++j){
        hPT[i][j]->Write();
        hETA[i][j]->Write(); 
        hPHI[i][j]->Write();  
        hMJDPHI[i][j]->Write();
        hISO[i][j]->Write();   
        hTM[i][j]->Write();
      }
}//-----saveHistos_MainAnalysis()
  
        
//------------------------------------------------------------------------
void saveHistos_MuonChargePt()
{
//------------------------------------------------------------------------
#if debugme
  cout << " Saving MuonChargePt histos " << endl;
#endif
    
  for(int j = 0; j != Num_trkAlgos; ++j)
    {
      hPTplus[j]->Write();
      hPTminus[j]->Write();
    }
  return;
}//-----saveHistos_MuonPtJetIso


//------------------------------------------------------------------------
void saveHistos(TFile * fout, string sample, int option)
{
//------------------------------------------------------------------------

  fout->cd(); 
  fout->mkdir(sample.c_str()); 
  fout->cd(sample.c_str());

  if(option == 1) saveHistos_MainAnalysis();

  else if(option == 2) saveHistos_MuonChargePt();
  else 
    {
      cout << " Something went wrong when saving with option " << option 
	   << " ...Quiting..." << endl;
      abort();
    }
  return;

}//saveHistos




//Routine to grab the information from a given file
//---------------------------------------------------------------------------
void gatherFileBasicInfo(const wprime::InputFile& file,
                         wprime::Event*& ev, int& nevents, float& weight)
{
//---------------------------------------------------------------------------

  ev = new wprime::Event();
  file.tree->SetBranchAddress("wp", &ev);
  nevents = file.tree->GetEntries();
  weight = file.weight;
#if debugme
  cout << " Number of events in the file = " << nevents << endl;
#endif
}//---------gatherFileBasicInfo()


// function signature: expects the event, a muon and an array with length equal to 
// Num_trkAlgos, to be updated with true/false, depending on selection cut;
//
// the function returs false if rest of selection cuts should be skipped (e.g. when
// the trigger has falied the event, or there are more than one muons in the event,
// or the jet-activity is vetoing the event; should return true if the rest of the 
// cuts should be executed (e.g. when quality cuts fail only for one particular
// tracking algorithm)
typedef bool (*funcPtr)(const wprime::Event *, const wprime::Muon*, bool * );

// key: cuts_desc_short[i], value: function pointer corresponding to selection cut
typedef map<string, funcPtr> selection_map;


// determine before event-loop the order in which cuts are to be executed
void setupCutOrder(selection_map & cuts)
{
  cout << "\n Cuts will be applied in this order: " << endl;
  for(int cut_i = 0; cut_i != Num_selection_cuts; ++cut_i)
    { // loop over selection cuts
      string arg = cuts_desc_short[cut_i];
      cout << " Cut #" << (cut_i+1) << ": " << cuts_desc_long[cut_i]
	   << " (" << arg << ") " << endl;
      
      if(arg == "hlt")cuts[arg] = &PassedHLT;
      else if(arg == "ptrange")cuts[arg] = &MuonPtWithinRange;
      else if(arg == "qual")cuts[arg] = &GoodQualityMuon;
      else if(arg == "1mu")cuts[arg] = &OnlyOneHighTrackPtMuon;
      else if(arg == "iso")cuts[arg] = &IsolatedMuon;
      else if(arg == "jet")cuts[arg] = &NoJetActivity;
      else
	{
	  cout << " Oops! Don't understand how to prepare for cut nicknamed as "
	       << arg << endl;
	  abort();
	}
    } // loop over selection cuts

  cout << endl;

}

//Main Analysis study
//---------------------------------------------------------------------------
void GetMuonPtDistribution(const wprime::InputFile& file, 
			   float Nexp_evt_cut[][Num_trkAlgos],
                           float & Nexp_evt, int option,
			   bool highestPtMuonOnly)
{
//---------------------------------------------------------------------------

  //gather file basic info
  wprime::Event * ev = 0;
  int nevents = 0; float weight = 0.0;
  gatherFileBasicInfo(file, ev, nevents, weight);
  
  // counter for (unweighted) events after cuts
  int Num_surv_cut[Num_selection_cuts][Num_trkAlgos] = {{0}};

  // the idea is to be able to swap selection cuts without having to modify both
  // the algo_desc_short array AND the code looping over the cuts 
  // in the very same order
  selection_map cuts;
  setupCutOrder(cuts);

  // Loop over events 
  for(int i = 0; i != nevents; ++i){ // event loop
#if debugme
    cout << " ########## Processing event #: " << i+1 << endl;
#endif
    
    file.tree->GetEntry(i);
    
    // switch to help us keep track of whether a muon has already
    // been found in current event surviving the ith-cut;
    // this will ensure that we increase Num_surv_cut maximum once per evet
    // whereas we nevertheless fill the histograms 
    // for every muon surviving the i-th cut
    bool accountme[Num_selection_cuts][Num_trkAlgos];
    for(int mm = 0; mm != Num_selection_cuts; ++mm)
      for(int nn = 0; nn != Num_trkAlgos; ++nn)
	accountme[mm][nn] = true;
    
    int iMuMin = 0; int iMuMax = ev->mu->GetLast() + 1;
    
    int iMu = -1;
    if(highestPtMuonOnly)
      {
	GetTheHardestMuon(ev, iMu);
	if(iMu >=0) // there is at least one muon in the event
	  {
	    iMuMin = iMu;
	    iMuMax = iMu + 1;
	  }
      }

    wprime::Muon* theMu = 0;
    //loop over muons
    for (int mi = iMuMin; mi != iMuMax; ++mi){//loop over muons
      
#if debugme
      cout << " ########## Processing muon #: " << mi+1 << endl;
#endif
      //get the muon
      theMu = (wprime::Muon *) ev->mu->At(mi);

      bool fill_entry[Num_trkAlgos] = {true, true, true, true, true, true};

#if dumpHighPtMuons
      bool HighPtMuon = false;
#endif

      for(int cut_index = 0; cut_index != Num_selection_cuts; ++cut_index)
	{ // loop over selection cuts
	  string arg = cuts_desc_short[cut_index];

	  // call to function [as implemented in setupCutOder]
	  bool survived_cut = (cuts[arg])(ev, theMu, fill_entry);
	  if(!survived_cut)break; // skip rest of selection cuts

	  tabulateMe(Num_surv_cut, cut_index, weight, ev, theMu, 
		     fill_entry, option, accountme);

#if dumpHighPtMuons
	  if(file.samplename == "data" && cut_index == Num_selection_cuts-1)
	    {
	      for(int i = 0; i != Num_trkAlgos; ++i)
		if(fill_entry[i])
		  {
		    HighPtMuon = true;
		    break;
		  }

	      //	      if(HighPtMuon || ev->run_no == 143657)
	      if(HighPtMuon)
		printHighPtMuon(ev, theMu, fill_entry);
	    }
#endif      

	} // loop over selection cuts

    } //muon loop
    
  } //event loop
  
  //Number of expected events for each cut (weighted)
  for(int ii = 0; ii != Num_selection_cuts; ++ii)
    for(int mual = 0; mual != Num_trkAlgos; ++mual)
      Nexp_evt_cut[ii][mual] += Num_surv_cut[ii][mual] * weight;
  
  //total # of events (before any cuts)
  Nexp_evt += nevents * weight; 
  
  delete ev;

}//-----------GetMuonPtDistribution()      

void printHighPtMuon(const wprime::Event * ev, const wprime::Muon* theMu,
		     const bool fill_entry[])
{
  cout << " Run # = " << ev->run_no << " Event # = " 
       << ev->evt_no << " LS = " << ev->LS_no << endl;
  cout << " glb pt = " << theMu->global.p.Pt()
       << " trk pt = " << theMu->tracker.p.Pt()
       << " tpfms pt = " << theMu->tpfms.p.Pt()
       << " ctail pt = " << theMu->cocktail.p.Pt()
       << " picky pt = " << theMu->picky.p.Pt()
       << " tmr pt = " << theMu->tmr.p.Pt() << endl;
  cout << " glb dpt = " << theMu->global.dpt
       << " trk dpt = " << theMu->tracker.dpt
       << " tpfms dpt = " << theMu->tpfms.dpt
       << " ctail dpt = " << theMu->cocktail.dpt
       << " picky dpt = " << theMu->picky.dpt
       << " tmr dpt = " << theMu->tmr.dpt << endl;
  cout << " eta =  " << theMu->tracker.p.Eta()
       << " phi = " << theMu->tracker.p.Phi() << endl;
  cout << " # of strip layers = " << theMu->tracker.Nstrip_layer 
       << ", # of pixel layers = " << theMu->tracker.Npixel_layer << endl;
  cout << " # of strip hits = " << theMu->tracker.NsiStrip_hits<<endl;
  cout << " # of strip layers w/o measurement = " 
       << theMu->tracker.Nstrip_layerNoMeas 
       << ", # of pixel layers w/o measurement= " 
       << theMu->tracker.Npixel_layerNoMeas << endl;
  cout << " # of strip hits = " << theMu->tracker.NsiStrip_hits
       << " # of pixel hits = " << theMu->tracker.Npixel_hits
       << " # of muon hits = " << theMu->tracker.Nmuon_hits
       << " # of standalone muon hits = " << theMu->Nmu_hits << endl;
  cout << " Global: " << theMu->AllGlobalMuons
       << " Tracker: " << theMu->AllTrackerMuons
       << " Standalone: " << theMu->AllStandAloneMuons
       << " Global prompt tight: " << theMu->GlobalMuonPromptTight << endl;

  cout << " Survives all cuts? ";
  for(int i = 0; i != Num_trkAlgos; ++i)
    cout << algo_desc_short[i] << ": " << fill_entry[i]
	 << " ";
  cout << endl << endl;
  
}


//---------------------------------------------------------------------------
void GetDistributionGeneric(const wprime::InputFile & file, 
                            TFile *fout, ofstream & out, 
                            const int option = 1, 
			    const bool highestPtMuonOnly = true)
{
//---------------------------------------------------------------------------
#if debugme
  cout << " $$$$$$$$$$$$$$$$$$$$GetDistributionGeneric " << endl;
#endif

  // Define histograms according to the type of study
  defineHistos(option);
  
  //initialize counters to be used in studies
  float Nexp_evt = 0;
  float Nexp_evt_cut[Num_selection_cuts][Num_trkAlgos] = {{0}};

  if(!file.tree)
    return;

  cout << "\n Processing sample " << file.description << endl;
  
  if (option == 1 || option == 2) 
    GetMuonPtDistribution(file, Nexp_evt_cut, Nexp_evt,option,
			  highestPtMuonOnly);
  else 
    {
      cout << " Option " << option << " is invalid, quiting..." << endl;
      abort();
    }
  
  // Print the results if needed according to study case
  if (option == 1 || option == 2) 
    printSummary_MuonPt(out, file.samplename, Nexp_evt, Nexp_evt_cut);
  else  
    cout << " Nothing to print for this study " << endl;
  
  //save histograms according to study case
  saveHistos(fout, file.samplename, option);  
}





