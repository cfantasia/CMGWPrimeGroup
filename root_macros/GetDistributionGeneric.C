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
void fillHistos_MuonPt(int index, float weight,const wprime::Muon* mu,
                       const bool fill_entry[]);
void fillHistos_MuonEta(int index, float weight,const wprime::Muon* mu,
                       const bool fill_entry[]);
void fillHistos_MuonPhi(int index, float weight,const wprime::Muon* mu,
                        const bool fill_entry[]);

void fillHistos_MuonJetDPhi(int index, float weight,
                        const wprime::Event * ev,
                        const wprime::Muon* mu,
                            const bool fill_entry[]);
void fillHistos_TMass(int index, float weight,
                        const wprime::Event * ev,
                        const wprime::Muon* mu,
                            const bool fill_entry[]);
void fillHistos_MuonIso(int index, float weight,
                        const wprime::Muon* mu,
                            const bool fill_entry[]);
void fillHistos_MuonChargePt(int index, float weight,const wprime::Muon* mu,
			     const bool fill_entry[]);
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
void defineHistos_TMass()
{
//--------------------------------------------------------------
#if debugme
  cout << " define TMass histos" << endl;
#endif

  for(int cut = 0; cut != Num_histo_sets; ++cut)
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

  for(int cut = 0; cut != Num_histo_sets; ++cut)
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

  for(int cut = 0; cut != Num_histo_sets; ++cut)
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

  for(int cut = 0; cut != Num_histo_sets; ++cut)
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

  for(int cut = 0; cut != Num_histo_sets; ++cut)
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

    if(option == 1){
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
  
  //fill the histograms
  if(option == 1) {
    fillHistos_MuonPt(cut_index,weight,mu, fill_entry);
    fillHistos_MuonEta(cut_index,weight,mu, fill_entry);
    fillHistos_MuonPhi(cut_index,weight,mu, fill_entry);
    fillHistos_MuonJetDPhi(cut_index,weight, ev, mu, fill_entry);
    fillHistos_MuonIso(cut_index,weight, mu, fill_entry);
    fillHistos_TMass(cut_index,weight, ev, mu, fill_entry);
  }
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


//-----------------------------------------------------------
void fillHistos_MuonEta(int index, float weight,const wprime::Muon* mu,
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
    hETA[index][0]->Fill(mu->global.p.Eta(), weight);
  
  if(fill_entry[1])
    hETA[index][1]->Fill(mu->tracker.p.Eta(), weight);
  
  if(fill_entry[2])
    hETA[index][2]->Fill(mu->tpfms.p.Eta(), weight);

  if(fill_entry[3])
    hETA[index][3]->Fill(mu->cocktail.p.Eta(), weight);

  if(fill_entry[4])
    hETA[index][4]->Fill(mu->picky.p.Eta(), weight);

  if(fill_entry[5])
    hETA[index][5]->Fill(mu->tmr.p.Eta(), weight);
  
}//fillHistos_MuonEta




//-----------------------------------------------------------
void fillHistos_MuonPhi(int index, float weight,const wprime::Muon* mu,
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
    hPHI[index][0]->Fill(mu->global.p.Phi(), weight);
  
  if(fill_entry[1])
    hPHI[index][1]->Fill(mu->tracker.p.Phi(), weight);
  
  if(fill_entry[2])
    hPHI[index][2]->Fill(mu->tpfms.p.Phi(), weight);

  if(fill_entry[3])
    hPHI[index][3]->Fill(mu->cocktail.p.Phi(), weight);

  if(fill_entry[4])
    hPHI[index][4]->Fill(mu->picky.p.Phi(), weight);

  if(fill_entry[5])
    hPHI[index][5]->Fill(mu->tmr.p.Phi(), weight);
  
}//fillHistos_MuonPhi




//-----------------------------------------------------------
void fillHistos_MuonJetDPhi(int index, float weight,
                        const wprime::Event * ev,
                        const wprime::Muon* mu,
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
    hMJDPHI[index][0]->Fill(XJetDPhi(mu->global.p, ev), weight);
  
  if(fill_entry[1])
    hMJDPHI[index][1]->Fill(XJetDPhi(mu->tracker.p, ev), weight);
  
  if(fill_entry[2])
    hMJDPHI[index][2]->Fill(XJetDPhi(mu->tpfms.p, ev), weight);

  if(fill_entry[3])
    hMJDPHI[index][3]->Fill(XJetDPhi(mu->cocktail.p, ev), weight);

  if(fill_entry[4])
    hMJDPHI[index][4]->Fill(XJetDPhi(mu->picky.p, ev), weight);

  if(fill_entry[5])
    hMJDPHI[index][5]->Fill(XJetDPhi(mu->tmr.p, ev), weight);
  
}//fillHistos_MuonJetDPhi



//-----------------------------------------------------------
void fillHistos_TMass(int index, float weight,
                        const wprime::Event * ev,
                        const wprime::Muon* mu,
                       const bool fill_entry[])
{
//-----------------------------------------------------------

  if(!mu) {
#if debugme
    cout << " No muon found in the event " << endl;
#endif
    return;
  }
  
  //for now just fill for global as we are not sure
  //how the met corrections are done
  if(fill_entry[0])
    hTM[index][0]->Fill(TMass(mu->global.p, ev->pfmet), weight);
  
  if(fill_entry[1])
    hTM[index][1]->Fill(TMass(mu->tracker.p, ev->pfmet), weight);
  
  if(fill_entry[2])
    hTM[index][2]->Fill(TMass(mu->tpfms.p, ev->pfmet), weight);

  if(fill_entry[3])
    hTM[index][3]->Fill(TMass(mu->cocktail.p, ev->pfmet), weight);

  if(fill_entry[4])
    hTM[index][4]->Fill(TMass(mu->picky.p, ev->pfmet), weight);

  if(fill_entry[5])
    hTM[index][5]->Fill(TMass(mu->tmr.p, ev->pfmet), weight);
  
}//fillHistos_TMass


//-----------------------------------------------------------
void fillHistos_MuonIso(int index, float weight,
                        const wprime::Muon* mu,
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
      hISO[index][0]->Fill(CombRelIsolation(mu, deltaRIsoIndex),weight);
  
  if(fill_entry[1])
    hISO[index][1]->Fill(CombRelIsolation(mu, deltaRIsoIndex), weight);
  
  if(fill_entry[2])
    hISO[index][2]->Fill(CombRelIsolation(mu, deltaRIsoIndex), weight);

  if(fill_entry[3])
    hISO[index][3]->Fill(CombRelIsolation(mu, deltaRIsoIndex), weight);

  if(fill_entry[4])
    hISO[index][4]->Fill(CombRelIsolation(mu, deltaRIsoIndex), weight);

  if(fill_entry[5])
    hISO[index][5]->Fill(CombRelIsolation(mu, deltaRIsoIndex), weight);
  
}//fillHistos_MuonIso







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
void saveHistos_MainAnalysis()
{
//------------------------------------------------------------------------
#if debugme
  cout << " Saving MuonPt histos " << endl;
#endif
  for(int i = 0; i != Num_histo_sets; ++i)
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
      tabulateMe(Num_surv_cut, cut_index, weight, ev, theMu, fill_entry,
		 option, accountme);
      
      //>>>>>>>>>>CUT 2
      CheckMuonPtInRange(theMu, fill_entry, minPtMu, maxPtMu);
      tabulateMe(Num_surv_cut, cut_index, weight, ev, theMu, fill_entry,
		 option, accountme);

             //>>>>>>>>>>CUT 3
      CheckQuality(theMu, fill_entry, PtTrackCut, Chi2Cut,Muon_Eta_Cut);
      tabulateMe(Num_surv_cut, cut_index, weight, ev, theMu, fill_entry,
		 option, accountme);

      //>>>>>>>>>>CUT 4
      if (!OnlyOneHighTrackPtMuon(ev, OneMuPtTrackCut)) continue;
      tabulateMe(Num_surv_cut, cut_index, weight, ev, theMu, fill_entry,
		 option, accountme);
      
      //>>>>>>>>>>CUT 5
      if (!Isolation(theMu, deltaRIsoIndex, CombRelCut)) continue;
      tabulateMe(Num_surv_cut, cut_index, weight, ev, theMu, fill_entry,
		 option, accountme);

      //>>>>>>>>>>CUT 6
      if (ExceedMaxNumJetsOpposedToMu(MaxNjetsAboveThresh, EtJetCut, 
				      Delta_Phi, theMu,ev)) continue;
      tabulateMe(Num_surv_cut, cut_index, weight, ev, theMu, fill_entry,
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
	      cout << " # of strip hits = " << theMu->tracker.NsiStrip_hits<<endl;
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





