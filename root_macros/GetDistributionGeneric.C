// ROOT49 stuff
#include <TROOT.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1F.h>

// wprime stuff
#include "UserCode/CMGWPrimeGroup/interface/wprimeEvent.h"
#include "UserCode/CMGWPrimeGroup/root_macros/input_file.h"
#include "UserCode/CMGWPrimeGroup/root_macros/loadCutsAndThresholds.h"

using std::cout; using std::endl; using std::vector; using std::string;

//methods in this macro:
void printSummary_MuonPt(ofstream & out, const string& dir, 
                         float Nexp_evt, 
			 const float Nexp_evt_cut[][Num_trkAlgos]);
void defineHistos_MuonPt();
void defineHistos_MuonChargePt();
void defineHistos(int option);
void tabulateMe(int Num_surv_cut[][Num_trkAlgos], int& cut_index, 
		float weight, const wprime::Muon* mu,
                const bool fill_entry[], int option, 
		bool accountme[][Num_trkAlgos]);
void fillHistos_MuonPt(int index, float weight,const wprime::Muon* mu,
                       const bool fill_entry[]);
void fillHistos_MuonChargePt(int index, float weight,const wprime::Muon* mu,
			     const bool fill_entry[]);
void saveHistos_MuonPt();
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
      for(int i = 0; i < Num_histo_sets; ++i){

	//print in a txt file the final total number of events
	if(i == Num_histo_sets - 1)
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

  for(int cut = 0; cut != Num_histo_sets; ++cut)
    for(int algo = 0; algo != Num_trkAlgos; ++algo)
      {
	string name = "hPT" + algo_desc_short[algo] + "_" + cuts_desc_short[cut];
	string title = algo_desc_long[algo] + " muon p_{T} with " + cuts_desc_long[cut];
	hPT[cut][algo] = new TH1F(name.c_str(), title.c_str(), nBinPtMu,minPtMu,maxPtMu);
      }
  
}//---------defineMuonPtHistos()


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
    defineHistos_MuonPt();
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
		float weight, const wprime::Muon* mu,
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
  
  //fill the histograms
  if(option == 1) 
    fillHistos_MuonPt(cut_index,weight,mu, fill_entry);
  else if(option == 2) 
    fillHistos_MuonChargePt(cut_index, weight, mu, fill_entry);
  else 
    cout << " WARNING!! I can't tabulate anything for option = " << option 
	 << endl;
  
  //since the event has passed the cut,
  //increase the cut_index for the next cut
  ++cut_index;

#if debugme
  cout << " cut_index is now = " << cut_index << endl;
#endif

}//-----tabulateMe



//Fill Histograms
//-----------------------------------------------------------
void fillHistos_MuonPt(int index, float weight,const wprime::Muon* mu,
                       const bool fill_entry[])
{
//-----------------------------------------------------------

  if(!mu) {
#if debugme
    cout << " No muon found in the event " << endl;
#endif
    return;
  }
  
  if(fill_entry[0])
    hPT[index][0]->Fill(mu->global.p.Pt(), weight);
  
  if(fill_entry[1])
    hPT[index][1]->Fill(mu->tracker.p.Pt(), weight);
  
  if(fill_entry[2])
    hPT[index][2]->Fill(mu->tpfms.p.Pt(), weight);

  if(fill_entry[3])
    hPT[index][3]->Fill(mu->cocktail.p.Pt(), weight);

  if(fill_entry[4])
    hPT[index][4]->Fill(mu->picky.p.Pt(), weight);

  if(fill_entry[5])
    hPT[index][5]->Fill(mu->tmr.p.Pt(), weight);
  
}//fillHistos_MuonPt


//Fill Histograms
//-----------------------------------------------------------
void fillHistos_MuonChargePt(int index, float weight, const wprime::Muon* mu,
			     const bool fill_entry[])
{
//-----------------------------------------------------------

  //do nothing if there is no muon in the event
  if(!mu) {
#if debugme
    cout << " No muon found in the event " << endl;
#endif
    return;
  }
  
  //only fill the required histograms
  if (index < (Num_histo_sets - Num_histo_sets_chargePt)) return;
  
  if (fill_entry[0]){
    if(mu->global.q > 0)
      hPTplus[0]->Fill(mu->global.p.Pt(), weight);
    else
      hPTminus[0]->Fill(mu->global.p.Pt(), weight);
  }
  if (fill_entry[1]){
    if(mu->tracker.q > 0)
      hPTplus[1]->Fill(mu->tracker.p.Pt(), weight);
    else
      hPTminus[1]->Fill(mu->tracker.p.Pt(), weight);
  }
  if (fill_entry[2]){
    if(mu->tpfms.q > 0)
      hPTplus[2]->Fill(mu->tpfms.p.Pt(), weight);
    else
      hPTminus[2]->Fill(mu->tpfms.p.Pt(), weight);
  }
  if (fill_entry[3]){
    if(mu->cocktail.q > 0)
      hPTplus[3]->Fill(mu->cocktail.p.Pt(), weight);
    else
      hPTminus[3]->Fill(mu->cocktail.p.Pt(), weight);
  }
  if (fill_entry[4]){
    if(mu->picky.q > 0)
      hPTplus[4]->Fill(mu->picky.p.Pt(), weight);
    else
      hPTminus[4]->Fill(mu->picky.p.Pt(), weight);
  }
  if (fill_entry[5]){
    if(mu->tmr.q > 0)
      hPTplus[5]->Fill(mu->tmr.p.Pt(), weight);
    else
      hPTminus[5]->Fill(mu->tmr.p.Pt(), weight);
  }
    
}//---------------fillHistos_MuonChargePt()



//------------------------------------------------------------------------
void saveHistos_MuonPt()
{
//------------------------------------------------------------------------
#if debugme
  cout << " Saving MuonPt histos " << endl;
#endif
    
  for(int i = 0; i != Num_histo_sets; ++i)
    for(int j = 0; j != Num_trkAlgos; ++j)
      hPT[i][j]->Write();
  
}//-----saveHistos_MuonPt


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

  if(option == 1) saveHistos_MuonPt();
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



//Main Analysis study
// TODO: Maybe rename the method to a more general one because
// in principle we can plot any variable in the event with this function, 
// not only muon pT. The tabulateMe function can be re-written to
// be a templated function.
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
  int Num_surv_cut[Num_histo_sets][Num_trkAlgos] = {{0}};
  
  //Loop over events:
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
    bool accountme[Num_histo_sets][Num_trkAlgos];
    for(int mm = 0; mm != Num_histo_sets; ++mm)
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
      
      //cut index to keep track of the order of cuts
      int cut_index = 0;
      //get the muon
      theMu = (wprime::Muon *) ev->mu->At(mi);
      
      bool fill_entry[Num_trkAlgos] = {true, true, true, true, true, true};
      
      // apply cuts
      //>>>>>>>>>>CUT 1
      if (!PassedHLT(ev)) continue;
      tabulateMe(Num_surv_cut, cut_index, weight, theMu, fill_entry,
		 option, accountme);
      
      //>>>>>>>>>>CUT 2
      CheckMuonPtInRange(theMu, fill_entry, minPtMu, maxPtMu);
      tabulateMe(Num_surv_cut, cut_index, weight, theMu, fill_entry,
		 option, accountme);
       
      //>>>>>>>>>>CUT 3
      if (!OnlyOneHighTrackPtMuon(ev, OneMuPtTrackCut)) continue;
      tabulateMe(Num_surv_cut, cut_index, weight, theMu, fill_entry,
		 option, accountme);
      
      //>>>>>>>>>>CUT 4
      //if (!SumPtIsolation(theMu, deltaRIsoIndex, SumPtCut)) continue;
      if (!CombRelIsolation(theMu, deltaRIsoIndex, CombRelCut)) continue;
      tabulateMe(Num_surv_cut, cut_index, weight, theMu, fill_entry,
		 option, accountme);

      //>>>>>>>>>>CUT 5
      if (ExceedMaxNumJetsOpposedToMu(MaxNjetsAboveThresh, EtJetCut, 
				      Delta_Phi, theMu,ev)) continue;
      tabulateMe(Num_surv_cut, cut_index, weight, theMu, fill_entry,
                 option, accountme);

      
      //>>>>>>>>>>CUT 6
      CheckQuality(theMu, fill_entry, PtTrackCut, Chi2Cut,Muon_Eta_Cut);
      tabulateMe(Num_surv_cut, cut_index, weight, theMu, fill_entry,
		 option, accountme);

#if dumpHighPtMuons
      if(file.samplename == "data")
	{
	  if(theMu->tracker.p.Pt() > minPtMu)
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
	      cout << " # of strip layers = " << theMu->tracker.Nstrip_layer << ", # of pixel layers = " << theMu->tracker.Npixel_layer << endl;
	      cout << " # of strip hits = " << theMu->tracker.NsiStrip_hits
	      cout << " # of strip layers w/o measurement = " << theMu->tracker.Nstrip_layerNoMeas << ", # of pixel layers w/o measurement= " << theMu->tracker.Npixel_layerNoMeas << endl;
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
	}
#endif
      
    }//muon loop
    
  }//event loop
    
  //Number of expected events for each cut (weighted)
  for(int ii = 0; ii != Num_histo_sets; ++ii)
    for(int mual = 0; mual != Num_trkAlgos; ++mual)
      Nexp_evt_cut[ii][mual] += Num_surv_cut[ii][mual] * weight;
  
  //total # of events (before any cuts)
  Nexp_evt += nevents * weight; 

  delete ev;
  
}//-----------GetMuonPtDistribution()


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
  float Nexp_evt_cut[Num_histo_sets][Num_trkAlgos] = {{0}};

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





