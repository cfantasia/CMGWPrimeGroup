#include "UserCode/CMGWPrimeGroup/interface/HadronicVZAnalyzer.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember)) 

using namespace std;
HadronicVZAnalyzer::HadronicVZAnalyzer(){}
HadronicVZAnalyzer::HadronicVZAnalyzer(const edm::ParameterSet & cfg, WPrimeUtil * wprimeUtil){

  wprimeUtil_ = wprimeUtil;
  assert(wprimeUtil_);

  SetLogFile(cfg.getParameter<string>("logFile"));
  SetCandEvtFile(cfg.getParameter<string>("candEvtFile"));

  intOptions_["report"] = cfg.getParameter<uint>("reportAfter");
  intOptions_["verbose"] = cfg.getParameter<bool>("debugme");
  doPreselect_ = cfg.getParameter<bool>("preselect");
  intOptions_["events"] = 0;

  debugme = cfg.getParameter<bool>("debugme");

  muonsLabel_ = cfg.getParameter<string>("muons");
  jetsLabel_ = cfg.getParameter<string>("jets");

  muonAlgo_ = cfg.getParameter<uint>("muonAlgo");
  if(debugme) cout<<"Using muon algo "<<muonAlgo_<<endl;

  
  hltEventLabel_ = cfg.getParameter<string>("hltEventTag");
  pileupLabel_ = cfg.getParameter<string>("pileupTag");

  triggersToUse_          = cfg.getParameter<vstring>("triggersToUse");

  PDGMUON = 13;
  PDGELEC = 11;
  PDGW = 24;
  PDGZ = 23;
  PDGZPRIME = 33;
  PDGWPRIME = 34;

  PI    = TMath::Pi();
  TWOPI = TMath::TwoPi();
  NOCUT = 9e9;
 
// +++++++++++++++++++Event characteristics
  

// +++++++++++++++++++General Cut values
  maxNumZs = cfg.getParameter<uint>("maxNumZs");
  minNLeptons = cfg.getParameter<uint>("minNLeptons");
  minNJets = cfg.getParameter<uint>("minNJets");
  maxNJets = cfg.getParameter<uint>("maxNJets");
  minLeadPt = cfg.getParameter<double>("minLeadPt");
  maxAngleBetweenJets = cfg.getParameter<double>("maxAngleBetweenJets");

// +++++++++++++++++++Z Cuts
  minZpt = cfg.getParameter<double>("minZpt");
  minZmass = cfg.getParameter<double>("minZmass");
  maxZmass = cfg.getParameter<double>("maxZmass");


// +++++++++++++++++++Hadronic Boson Cuts
  minHadVpt = cfg.getParameter<double>("minHadVPt");
  minHadVmass = cfg.getParameter<double>("minHadVmass");
  maxHadVmass = cfg.getParameter<double>("maxHadVmass");

// +++++++++++++++++++Muon General Cuts
  maxMuonEta = cfg.getParameter<double>("maxMuonEta");
  minMuonLoosePt = cfg.getParameter<double>("minMuonLoosePt");
  minMuonTightPt = cfg.getParameter<double>("minMuonTightPt");
//VBTF Recommended Cuts
  maxMuonDxy = cfg.getParameter<double>("maxMuonDxy");
  maxMuonNormChi2 = cfg.getParameter<double>("maxMuonNormChi2");
  minMuonNPixHit = cfg.getParameter<uint>("minMuonNPixHit");
  minMuonNTrkHit = cfg.getParameter<uint>("minMuonNTrkHit");
  minMuonStations = cfg.getParameter<uint>("minMuonStations");
  minMuonHitsUsed = cfg.getParameter<uint>("minMuonHitsUsed");

// +++++++++++++++++++Jet General Cuts
  minJetPt = cfg.getParameter<double>("minJetPt");
  maxJetEta = cfg.getParameter<double>("maxJetEta");
  maxJetNHF = cfg.getParameter<double>("maxJetNHF");
  maxJetNEF = cfg.getParameter<double>("maxJetNEF");
  minJetCHF = cfg.getParameter<double>("minJetCHF");
  maxJetCEF = cfg.getParameter<double>("maxJetCEF");
  minJetnumConst = cfg.getParameter<uint>("minJetnumConst");
  minJetcMult = cfg.getParameter<uint>("minJetcMult");
}

HadronicVZAnalyzer::~HadronicVZAnalyzer(){
  outCandEvt.close(); 
  outLogFile.close(); 
}


/// Declare Histograms
//--------------------------------------------------------------
void HadronicVZAnalyzer::Declare_Histos(TFileDirectory & dir)
{
  // Extend later.
  printf("Declare histos\n");
  h_HadVZMass = dir.make<TH1F>("h_HadVZMass","h_HadVZMass",250,0.0,2500.0);

  //New histos
  
  h_Zmuon1_pt = dir.make<TH1F>("h_Zmuon1_pt", "h_Zmuon1_pt", 100, 0.0, 1000.0);
  h_Zmuon1_eta = dir.make<TH1F>("h_Zmuon1_eta", "h_Zmuon1_eta", 40, -5.0, 5.0);
  h_Zmuon1_phi = dir.make<TH1F>("h_Zmuon1_phi", "h_Zmuon1_phi", 20, -4.0, 4.0);
  h_Zmuon2_pt = dir.make<TH1F>("h_Zmuon2_pt", "h_Zmuon2_pt", 100, 0.0, 1000.0);
  h_Zmuon2_eta = dir.make<TH1F>("h_Zmuon2_eta", "h_Zmuon2_eta", 40, -5.0, 5.0);
  h_Zmuon2_phi = dir.make<TH1F>("h_Zmuon2_phi", "h_Zmuon2_phi", 20, -4.0, 4.0);


  h_Zmuon1_VZCut_pt = dir.make<TH1F>("h_Zmuon1_VZCut_pt", "h_Zmuon1_VZCut_pt", 100, 0.0, 1000.0);
  h_Zmuon1_VZCut_eta = dir.make<TH1F>("h_Zmuon1_VZCut_eta", "h_Zmuon1_VZCut_eta", 40, -5.0, 5.0);
  h_Zmuon1_VZCut_phi = dir.make<TH1F>("h_Zmuon1_VZCut_phi", "h_Zmuon1_VZCut_phi", 20, -4.0, 4.0);
  h_Zmuon2_VZCut_pt = dir.make<TH1F>("h_Zmuon2_VZCut_pt", "h_Zmuon2_VZCut_pt", 100, 0.0, 1000.0);
  h_Zmuon2_VZCut_eta = dir.make<TH1F>("h_Zmuon2_VZCut_eta", "h_Zmuon2_VZCut_eta", 40, -5.0, 5.0);
  h_Zmuon2_VZCut_phi = dir.make<TH1F>("h_Zmuon2_VZCut_phi", "h_Zmuon2_VZCut_phi", 20, -4.0, 4.0);

  h_deltaR_muon1muon2 = dir.make<TH1F>("h_deltaR_muon1muon2", "h_deltaR_muon1muon2", 50, 0., 5.);

  h_deltaR_HadVmuon1 = dir.make<TH1F>("h_deltaR_HadVmuon1", "h_deltaR_HadVmuon1", 50, 0., 5.);
  h_deltaR_HadVmuon2 = dir.make<TH1F>("h_deltaR_HadVmuon2", "h_deltaR_HadVmuon2", 50, 0., 5.);



  h_jet1_pt = dir.make<TH1F>("h_jet1_pt", "h_jet1_pt", 100, 0.0, 1000.0);
  h_jet1_eta = dir.make<TH1F>("h_jet1_eta", "h_jet1_eta", 40, -5.0, 5.0);
  h_jet1_phi = dir.make<TH1F>("h_jet1_phi", "h_jet1_phi", 20, -4.0, 4.0);
  h_jet2_pt = dir.make<TH1F>("h_jet2_pt", "h_jet2_pt", 100, 0.0, 1000.0);
  h_jet2_eta = dir.make<TH1F>("h_jet2_eta", "h_jet2_eta", 40, -5.0, 5.0);
  h_jet2_phi = dir.make<TH1F>("h_jet2_phi", "h_jet2_phi", 20, -4.0, 4.0);
  h_jet1_mass = dir.make<TH1F>("h_jet1_mass", "h_jet1_mass", 100, 0.0, 500.);
  h_jet2_mass = dir.make<TH1F>("h_jet2_mass", "h_jet2_mass", 100, 0.0, 500.);


  h_deltaR_jet1muon1 = dir.make<TH1F>("h_deltaR_jet1muon1", "h_deltaR_jet1muon1", 50, 0., 5.);
  h_deltaR_jet1muon2 = dir.make<TH1F>("h_deltaR_jet1muon2", "h_deltaR_jet1muon2", 50, 0., 5.);
  h_deltaR_jet2muon1 = dir.make<TH1F>("h_deltaR_jet2muon1", "h_deltaR_jet2muon1", 50, 0., 5.);
  h_deltaR_jet2muon2 = dir.make<TH1F>("h_deltaR_jet2muon2", "h_deltaR_jet2muon2", 50, 0., 5.);

  h_HadVZpt = dir.make<TH1F>("h_HadVZpt", "h_HadVZpt", 30, 0.0, 300.0);
  h_HadVZeta = dir.make<TH1F>("h_HadVZeta", "h_HadVZeta", 40, -5.0, 5.0);
  h_HadVZphi = dir.make<TH1F>("h_HadVZphi", "h_HadVZphi", 20, -4.0, 4.0);
  
  h_muons_pt = dir.make<TH1F>("h_muons_pt", "h_muons_pt", 100, 0.0, 1000.);
  h_muons_eta = dir.make<TH1F>("h_muons_eta", "h_muons_eta", 40, -5.0, 5.);  
  h_muons_phi = dir.make<TH1F>("h_muons_phi", "h_muons_phi", 20, -4.0, 4.);


  h_jets_pt = dir.make<TH1F>("h_jets_pt", "h_jets_pt", 100, 0.0, 1000.);
  h_jets_eta = dir.make<TH1F>("h_jets_eta", "h_jets_eta", 40, -5.0, 5.);  
  h_jets_phi = dir.make<TH1F>("h_jets_phi", "h_jets_phi", 20, -4.0, 4.);
  

  h_jet_HadV_pt = dir.make<TH1F>("h_jet_HadV_pt", "h_jet_HadV_pt", 100, 0.0, 1000.);
  h_jet_HadV_eta = dir.make<TH1F>("h_jet_HadV_eta", "h_jet_HadV_eta", 40, -5.0, 5.);  
  h_jet_HadV_phi = dir.make<TH1F>("h_jet_HadV_phi", "h_jet_HadV_phi", 20, -4.0, 4.);


  h_jet_VZCut_pt = dir.make<TH1F>("h_jet_VZCut_pt", "h_jet_VZCut_pt", 100, 0.0, 1000.);
  h_jet_VZCut_eta = dir.make<TH1F>("h_jet_VZCut_eta", "h_jet_VZCut_eta", 40, -5.0, 5.);  
  h_jet_VZCut_phi = dir.make<TH1F>("h_jet_VZCut_phi", "h_jet_VZCut_phi", 20, -4.0, 4.);
  

  h_jet_mult = dir.make<TH1F>("h_jet_mult", "h_jet_mult", 10, -0.5, 9.5); 
  h_jet_mult_inc = dir.make<TH1F>("h_jet_mult_incl", "h_jet_mult_incl", 10, -0.5, 9.5); 


  // 2D Piotr histos

  h_dptpt_vs_pt = dir.make<TH2F>("h_dptpt_vs_pt", "h_dptpt_vs_pt", 50, 0., 0.5, 50, 0.0, 500.);
  h_dptpt2_vs_pt = dir.make<TH2F>("h_dptpt2_vs_pt", "h_dptpt2_vs_pt", 150, 0., 0.1, 150, 0.0, 500.);

  h_dptpt_vs_invpt = dir.make<TH2F>("h_dptpt_vs_invpt","h_dptpt_vs_invpt", 50, 0., 0.5, 50, 0., 0.1);
  h_dptpt2_vs_invpt = dir.make<TH2F>("h_dptpt2_vs_invpt","h_dptpt2_vs_invpt", 150, 0., 0.5, 150, 0., 0.1);
  h_dptpt_vs_eta = dir.make<TH2F>("h_dptpt_vs_eta","h_dptpt_vs_eta", 50, 0., 0.5, 50, -5., 5.);
  h_dptpt2_vs_eta = dir.make<TH2F>("h_dptpt2_vs_eta","h_dptpt2_vs_eta", 150, 0., 0.5, 150, -5., 5.);

  //



}//Declare_Histos

//Fill Histograms
//-----------------------------------------------------------
void HadronicVZAnalyzer::Fill_Histos(int index, float weight)
{
  // For now we have only one histo, so no need for this function.
  printf("Filling Histos\n");
}//Fill_Histos

void 
HadronicVZAnalyzer::eventLoop(edm::EventBase const & event){
  ClearEvtVariables();

  // Preselection - skip events that don't look promising
  if (doPreselect_){
    if(debugme) cout<<"Testing Preselection...\n";
    // We could setup some preselection here. To be implemented.
  }

  //  rhoFastJet_ = getProduct<double>(event,"kt6PFJets:rho");


  // Get leptons
  const vector<pat::Muon    > patMuons     = getProduct<vector<pat::Muon    > >(event, muonsLabel_);
  // Get jets
  const vector<pat::Jet     > patJets      = getProduct<vector<pat::Jet     > >(event, jetsLabel_);

  if(debugme){
    printf("    Contains: %i muon(s)\n",
	   (int)patMuons.size());
    printf("    Contains: %i jet(s)\n",
	   (int)patJets.size());
  }



  // Make vectors of leptons passing various criteria
  // Loop over muons, and see if they pass the TeVMuon criteria  
  for (size_t i = 0; i < patMuons.size(); i++) {
    muons_.push_back(TeVMuon(patMuons[i],muonAlgo_));   


    //New simple implementation of muon cuts
    //Chi2 cut dropped

    bool TightMuonCuts = PassMuonTightPtCut(&muons_[i]) && PassMuonGlobalCut(&muons_[i]) && PassMuonNpixhitCut(&muons_[i]) && PassMuonNtrkhitCut(&muons_[i]) && PassMuonHitsUsedCut(&muons_[i]) && PassMuonStationsCut(&muons_[i]) && PassMuonEtaCut(&muons_[i]) && PassMuonDxyCut(&muons_[i]);

    //  bool LooseMuonCuts = PassMuonLoosePtCut(&muons_[i]) && PassMuonGlobalCut(&muons_[i]) && PassMuonNpixhitCut(&muons_[i]) && PassMuonNtrkhitCut(&muons_[i]) && PassMuonHitsUsedCut(&muons_[i]) && PassMuonStationsCut(&muons_[i]) && PassMuonEtaCut(&muons_[i]) && PassMuonDxyCut(&muons_[i]);

    //    bool LooseMuonCuts = PassMuonLoosePtCut(&muons_[i]) && PassMuonGlobalCut(&muons_[i]) && PassMuonNpixhitCut(&muons_[i]) && PassMuonNtrkhitCut(&muons_[i]) && PassMuonHitsUsedCut(&muons_[i]) && PassMuonStationsCut(&muons_[i]) && PassMuonEtaCut(&muons_[i]) && PassMuonDxyCut(&muons_[i]);
    
bool LooseMuonCuts = PassMuonLoosePtCut(&muons_[i]) && PassMuonGlobalCut(&muons_[i]) && PassMuonNpixhitCut(&muons_[i]) && PassMuonHitsUsedCut(&muons_[i]) && PassMuonStationsCut(&muons_[i]) && PassMuonEtaCut(&muons_[i]) && PassMuonDxyCut(&muons_[i]);


    if (LooseMuonCuts)
      looseMuons_.push_back(muons_[i]);
    
    if (TightMuonCuts)
      tightMuons_.push_back(muons_[i]);




    //    if (PassMuonCut(&muons_[i]))
    //  looseMuons_.push_back(muons_[i]);

    //New stuff

    /*    const float pu = MuonPU(muons_[i]);
    if (looseMuon_(muons_[i], muonResult_, pu))
      looseMuons_.push_back(muons_[i]);

    if (tightMuon_(muons_[i], muonResult_, pu))
      tightMuons_.push_back(muons_[i]);
    //
    */
  }

  if (debugme){
    printf("    Contains: %i looseMuons_(s)\n",
	   (int)looseMuons_.size());

    printf("    Contains: %i tightMuons_(s)\n",
	   (int)tightMuons_.size());

  }


  if (looseMuons_.size() > 0)
    {
      sort(looseMuons_.begin(), looseMuons_.end(), highestMuonPt());

      if(debugme)
	{
	  for (uint k=0; k<looseMuons_.size(); k++)
	    {
	      cout << "looseMuons_ " << k << " has pt of " << looseMuons_.at(k).pt() <<  endl;
	    }
	}

    }

  if (tightMuons_.size() > 0)
    {
      sort(tightMuons_.begin(), tightMuons_.end(), highestMuonPt());
    }
  
  //Fill histos for all muons in the event
  for (size_t nM=0; nM<looseMuons_.size(); nM++)
    {
      h_muons_pt->Fill(looseMuons_.at(nM).pt());
      h_muons_eta->Fill(looseMuons_.at(nM).eta());
      h_muons_phi->Fill(looseMuons_.at(nM).phi());

      double trackpT = looseMuons_.at(nM).track()->pt();
      double trackpT2 = trackpT*trackpT;
      double invpt = 1/trackpT;
      double trackpTError = looseMuons_.at(nM).track()->ptError();

      double dptpt = trackpTError/trackpT;
      double dptpt2 =trackpTError/trackpT2;

      //  cout << "Muon number " << nM << " in the event had ptError = " << trackpTError << " and track pt " << trackpT << " and dpt/pt " << dptpt << " and dpt/pt2 " << dptpt2 << endl;

      h_dptpt_vs_pt->Fill(dptpt, trackpT);
      h_dptpt2_vs_pt->Fill(dptpt2, trackpT);

      h_dptpt_vs_invpt->Fill(dptpt, invpt);
      h_dptpt2_vs_invpt->Fill(dptpt2, invpt);

      h_dptpt_vs_eta->Fill(dptpt, looseMuons_.at(nM).track()->eta());
      h_dptpt2_vs_eta->Fill(dptpt2, looseMuons_.at(nM).track()->eta());


      /*
      h_dptpt_vs_pt->Fill((looseMuons_.at(nM).track()->ptError())/(looseMuons_.at(nM).track()->pt()),looseMuons_.at(nM).pt());
      h_dptpt2_vs_pt->Fill((looseMuons_.at(nM).track()->ptError())/((looseMuons_.at(nM).track()->pt())*(looseMuons_.at(nM).track()->pt())),looseMuons_.at(nM).pt());
      */
      
    }


  //Fill muon histos for loose muons
  int oneMu = 0;
  int twoMu = 0;
  if (looseMuons_.size() > 0)
    {
      /*      h_Zmuon1_pt->Fill(looseMuons_.at(0).pt());
      h_Zmuon1_eta->Fill(looseMuons_.at(0).eta());
      h_Zmuon1_phi->Fill(looseMuons_.at(0).phi());*/
      oneMu=1;
      if (looseMuons_.size() > 1)
	{
	  /*
	  h_Zmuon2_pt->Fill(looseMuons_.at(1).pt());
	  h_Zmuon2_eta->Fill(looseMuons_.at(1).eta());
	  h_Zmuon2_phi->Fill(looseMuons_.at(1).phi());	  */
	  twoMu=1;
	  //  h_deltaR_muon1muon2->Fill(deltaR(looseMuons_.at(0).eta(), looseMuons_.at(0).phi(), looseMuons_.at(1).eta(), looseMuons_.at(1).phi()));

	}

    }
  

  // Loop over jets, and see if they pass the jet criteria
  for (size_t i = 0; i < patJets.size(); ++i) {
    if (PassJetCut(&patJets[i]))
      jets_.push_back(patJets[i]);
  }
  
  h_jet_mult->Fill(jets_.size());


  h_jet_mult_inc->Fill(0);
    

  if (jets_.size()>=1)
    {
      h_jet_mult_inc->Fill(1);
    }


  if (jets_.size()>=2)
    {
      h_jet_mult_inc->Fill(2);
    }

  if (jets_.size()>=3)
    {
      h_jet_mult_inc->Fill(3);
    }

  if (jets_.size()>=4)
    {
      h_jet_mult_inc->Fill(4);
    }

  if (jets_.size()>=5)
    {
      h_jet_mult_inc->Fill(5);
    }

  if (jets_.size()>=6)
    {
      h_jet_mult_inc->Fill(6);
    }


  if (jets_.size()>=7)
    {
      h_jet_mult_inc->Fill(7);
    }

  if (jets_.size()>=8)
    {
      h_jet_mult_inc->Fill(8);
    }

  if (jets_.size()>=9)
    {
      h_jet_mult_inc->Fill(9);
    }

  if (jets_.size()>=10)
    {
      h_jet_mult_inc->Fill(10);
    }




  for (size_t nJ=0; nJ<jets_.size(); nJ++)
    {
      h_jets_pt->Fill(jets_.at(nJ).pt());
      h_jets_eta->Fill(jets_.at(nJ).eta());
      h_jets_phi->Fill(jets_.at(nJ).phi());
    }



  //Fill jet histos for jets who passes the criteria
  if (jets_.size() > 0)
    {
      sort(jets_.begin(), jets_.end(), highestJetPt());
      if(debugme)
	{
	  for (uint k=0; k<jets_.size(); k++)
	    {
	      cout << "jets_ " << k << " has pt of " << jets_.at(k).pt() <<  endl;
	    }
	}
    }


  if (jets_.size() > 0)
    {
      h_jet1_pt->Fill(jets_.at(0).pt());
      h_jet1_eta->Fill(jets_.at(0).eta());    
      h_jet1_phi->Fill(jets_.at(0).phi());
      h_jet1_mass->Fill(jets_.at(0).mass());

      if (oneMu==1)
	{
	  h_deltaR_jet1muon1->Fill(deltaR(jets_.at(0).eta(), jets_.at(0).phi(), looseMuons_.at(0).eta(), looseMuons_.at(0).phi()));
	}
      if (twoMu ==1)
	{
	  h_deltaR_jet1muon2->Fill(deltaR(jets_.at(0).eta(), jets_.at(0).phi(), looseMuons_.at(1).eta(), looseMuons_.at(1).phi()));
	  }
      
      if (jets_.size() > 1)
	{
	  h_jet2_pt->Fill(jets_.at(1).pt());
	  h_jet2_eta->Fill(jets_.at(1).eta());    
	  h_jet2_phi->Fill(jets_.at(1).phi());
	  h_jet2_mass->Fill(jets_.at(1).mass());
	  if (oneMu==1)
	    {
	      h_deltaR_jet2muon1->Fill(deltaR(jets_.at(1).eta(), jets_.at(1).phi(), looseMuons_.at(0).eta(), looseMuons_.at(0).phi()));
	    }
	  if (twoMu ==1)
	    {
	      h_deltaR_jet2muon2->Fill(deltaR(jets_.at(1).eta(), jets_.at(1).phi(), looseMuons_.at(1).eta(), looseMuons_.at(1).phi()));
	    }
	  
	}
      

    }



  if (debugme){
    printf("    Contains: %i jets_(s)\n",
	   (int)jets_.size());
  }

  int passedNLeptons = 0;
  int passedNJets = 0;

  if (looseMuons_.size() > minNLeptons)
    {
      passedNLeptons = 1;
    }

  if (jets_.size() > minNJets)
    {
      passedNJets = 1;
    }

  if (debugme){
    if (passedNLeptons==1)
      cout << "Passed #leptons cut" << endl;
    else 
      cout << "Didnt pass #leptons cut" << endl;

    if (passedNJets==1)
      cout << "Passed #jets cut" << endl;
    else
      cout << "Didnt pass #jets cut" << endl;

  }
  
  if(!(passedNLeptons==1 && passedNJets==1)){
  
    if(debugme) 
      cout << "No leptons or jets. Bad bad event, returning now..." << endl;
    return;
  }

  if (debugme) 
    cout << "Im just before boson cands" << endl;

  // Make a Z candidate out of the muons. 
  int oppCharge=0;

  for (size_t i=0; i<looseMuons_.size(); i++){
    if (debugme) 
      cout << "looseMuons_ has charge "  << looseMuons_.at(i).charge() << endl;
    for (size_t j=i+1; j<looseMuons_.size(); j++){
      if (looseMuons_[i].charge() != looseMuons_[j].charge())
	oppCharge=1;
    }
  }

  if (oppCharge==1){
    //  zCand = getZCands(looseMuons_).front();
    //    zCand = getZCands(looseMuons_).at(0);

    ZCandV zCand = getZCands(looseMuons_, 25.0);
    
    if(zCand.size()>0)
      { 
	zCand_ =  zCand[0];
	if (debugme)
	  {
	    cout << "Daugh 1 pT is " << zCand[0].daughter(0)->pt() << endl;
	    cout << "Daugh 2 pT is " << zCand[0].daughter(1)->pt() << endl;
	  }
	
	/*	const TeVMuon m1 = FindMuon(*zCand_.daughter(0));
	const TeVMuon m2 = FindMuon(*zCand_.daughter(1));
	h_Zmuon1_pt->Fill(m1.pt());
	h_Zmuon1_eta->Fill(m1.eta());
	h_Zmuon1_phi->Fill(m1.phi());
	h_Zmuon2_pt->Fill(m2.pt());
	h_Zmuon2_eta->Fill(m2.eta());
	h_Zmuon2_phi->Fill(m2.phi());	  
	h_deltaR_muon1muon2->Fill(deltaR(m1.eta(), m1.phi(), m2.eta(), m2.phi()));*/
	    
	
      }

    else                           
      zCand_ = ZCandidate();
   
    // zCand_ = zCand.size() ? zCand[0] : ZCandidate();

    if (debugme)
      cout << "Passed zCand" << endl;
  }
 
  // Make a W candidate out of the jets.

  wCand = getWCand(jets_);

  if (debugme)
    cout << "Passed wCand" << endl;


  bool validZ    = PassValidZCut();
  if (debugme && validZ)
    cout << "Valid Z from muon" << endl;

  bool validHadV = PassValidHadVCut();
  if(debugme && validHadV)
    cout << "Valid Z from jet" << endl;


  bool goodZ     = PassZMassCut() && PassZptCut();
  if (debugme && goodZ)
    cout << "Good Z from muon" << endl;


  bool goodHadV  = PassHadVMassCut() && PassHadVptCut();
  if (debugme && goodHadV)
    cout << "Good Z from jet" << endl;



  if (goodZ)
    {
      const TeVMuon m1 = FindMuon(*zCand_.daughter(0));
      const TeVMuon m2 = FindMuon(*zCand_.daughter(1));
      h_Zmuon1_pt->Fill(m1.pt());
      h_Zmuon1_eta->Fill(m1.eta());
      h_Zmuon1_phi->Fill(m1.phi());
      h_Zmuon2_pt->Fill(m2.pt());
      h_Zmuon2_eta->Fill(m2.eta());
      h_Zmuon2_phi->Fill(m2.phi());	  
      h_deltaR_muon1muon2->Fill(deltaR(m1.eta(), m1.phi(), m2.eta(), m2.phi()));
    }

  if (goodHadV)
    {
      h_jet_HadV_pt->Fill(wCand.pt());
      h_jet_HadV_eta->Fill(wCand.eta());
      h_jet_HadV_phi->Fill(wCand.phi());
    }

  
  if(validZ && validHadV && goodZ && goodHadV) {
    reco::CompositeCandidate hadVZ;
    hadVZ.addDaughter(zCand_);
    hadVZ.addDaughter(wCand);
    AddFourMomenta addP4;
    addP4.set(hadVZ);
    h_HadVZMass->Fill(hadVZ.mass());
    h_HadVZpt->Fill(hadVZ.pt());
    h_HadVZeta->Fill(hadVZ.eta());
    h_HadVZphi->Fill(hadVZ.phi());


    const TeVMuon VZm1 = FindMuon(*zCand_.daughter(0));
    const TeVMuon VZm2 = FindMuon(*zCand_.daughter(1));
    h_Zmuon1_VZCut_pt->Fill(VZm1.pt());
    h_Zmuon1_VZCut_eta->Fill(VZm1.eta());
    h_Zmuon1_VZCut_phi->Fill(VZm1.phi());
    h_Zmuon2_VZCut_pt->Fill(VZm2.pt());
    h_Zmuon2_VZCut_eta->Fill(VZm2.eta());
    h_Zmuon2_VZCut_phi->Fill(VZm2.phi());
    
    h_jet_VZCut_pt->Fill(wCand.pt());
    h_jet_VZCut_eta->Fill(wCand.eta());
    h_jet_VZCut_phi->Fill(wCand.phi());

    h_deltaR_HadVmuon1->Fill(deltaR(wCand.eta(), wCand.phi(), VZm1.eta(), VZm1.phi()));
    h_deltaR_HadVmuon2->Fill(deltaR(wCand.eta(), wCand.phi(), VZm2.eta(), VZm2.phi()));


  }
  else{
    if (debugme)
      cout << "Didnt manage to make a HadVZ Candidate" << endl;
  }

  
}


/////////////////Accessors///////////////////////

/////////////////Modifies///////////////////////
/*void HadronicVZAnalyzer::CheckStream(ofstream& stream, string s){
  if(!stream) { 
    cout << "Cannot open file " << s << endl; 
    abort();
  } 
}
*/
void HadronicVZAnalyzer::SetCandEvtFile(string s){
  outCandEvt.open(s.c_str());
  WPrimeUtil::CheckStream(outCandEvt, s);
}

void HadronicVZAnalyzer::SetLogFile(string s){
  outLogFile.open(s.c_str());      
  WPrimeUtil::CheckStream(outLogFile, s); 
}

/////////////////Cuts///////////////////////

//Always true
bool HadronicVZAnalyzer::PassNoCut()
{
  return true;
}

//Trigger requirements
//-----------------------------------------------------------
bool HadronicVZAnalyzer::PassTriggersCut()
{
//-----------------------------------------------------------
  if(debugme) cout<<"Trigger requirements"<<endl;
  // To be implemented
  return true;
}//--- PassTriggersCut()

bool
HadronicVZAnalyzer::PassNLeptonsCut(){
  return (looseMuons_.size()) > minNLeptons;
}

bool
HadronicVZAnalyzer::PassNJetsCut(){
  return ((jets_.size() > minNJets) && (jets_.size() < maxNJets));
}

bool
HadronicVZAnalyzer::PassValidHadVCut(){
  return wCand && wCand.mass()>0.;
}

bool
HadronicVZAnalyzer::PassValidZCut(){
  return zCand_ && zCand_.mass()>0.;
}

bool
HadronicVZAnalyzer::PassValidHadVZCandCut(){
  // To be implemented
  return true;
}

bool
HadronicVZAnalyzer::PassNumberOfZsCut(){
  return numZs < maxNumZs;
}

bool
HadronicVZAnalyzer::PassLeadingLeptonPtCut(){
  return LeadPt > minLeadPt;
}

////////////////////////////////
/////////Check Z Properties/////
////////////////////////////////
bool
HadronicVZAnalyzer::PassZMassCut(){
  return (zCand_.mass() > minZmass) && (zCand_.mass() < maxZmass);  
}

bool
HadronicVZAnalyzer::PassZptCut(){
  return zCand_.pt() > minZpt;
}
////////////////////////////////
///////Check Had. V Properties//
////////////////////////////////

bool HadronicVZAnalyzer::PassHadVMassCut(){
  return (wCand.mass() > minHadVmass) && (wCand.mass() < maxHadVmass);
}

bool
HadronicVZAnalyzer::PassHadVptCut(){
  return wCand.pt() > minHadVpt;
}

////////////////////////////////
//////Check Muon Properties/////
////////////////////////////////
/// Big function to check all cuts.
bool HadronicVZAnalyzer::PassMuonCut(const TeVMuon* mu){
  for(uint i=0; i<MuonCutFns_.size(); ++i){
    bool result = CALL_MEMBER_FN(*this,MuonCutFns_[i])(mu);
    if(result==false) return false;
  }
  return true;
}

bool HadronicVZAnalyzer::PassMuonLoosePtCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Loose Pt Cut"<<endl;
  return (mu->pt() > minMuonLoosePt);
}

bool HadronicVZAnalyzer::PassMuonTightPtCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Tight Pt Cut"<<endl;
  return (mu->pt() > minMuonTightPt);
}//--- PassMuonPtCut

bool HadronicVZAnalyzer::PassMuonGlobalCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Global Cut"<<endl;
  return (mu->isGlobalMuon()); 
}//--- PassMuonGlobalCut

bool HadronicVZAnalyzer::PassMuonNpixhitCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon NpixhitCut"<<endl;
  return (mu->globalTrack()->hitPattern().numberOfValidPixelHits() > minMuonNPixHit);
}//--- PassMuonNpixhitCut

bool HadronicVZAnalyzer::PassMuonNtrkhitCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon NtrkhitCut"<<endl;
  return (mu->globalTrack()->hitPattern().numberOfValidTrackerHits() > minMuonNTrkHit);
}//--- PassMuonNtrkhitCut

bool HadronicVZAnalyzer::PassMuonNormChi2Cut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Chi2 Cut"<<endl;
  return (mu->globalTrack()->normalizedChi2() < maxMuonNormChi2);
}//--- PassMuonChi2Cut

bool HadronicVZAnalyzer::PassMuonHitsUsedCut(const TeVMuon* mu){
  //Num Valid Muon Hits
  if(debugme) cout<<"Check Muon Hits Used Cut"<<endl;
  return (mu->globalTrack()->hitPattern().numberOfValidMuonHits() > minMuonHitsUsed);
}//--- PassMuonHits Used Cut

bool HadronicVZAnalyzer::PassMuonStationsCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Stations Cut"<<endl;
  return (mu->numberOfMatches() > minMuonStations);
}//--- PassMuonStationsCut

bool HadronicVZAnalyzer::PassMuonEtaCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Eta Cut"<<endl;
  return (fabs(mu->eta()) < maxMuonEta);
}//--- PassMuonEta Cut

/*
bool HadronicVZAnalyzer::PassMuonCombRelIsoCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon CombRelIso Cut"<<endl;
  return (Calc_MuonRelIso(mu) < maxWmunuCombRelIso);
}//--- PassMuonCombRelIsoCut
*/

bool HadronicVZAnalyzer::PassMuonDxyCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Dxy Cut"<<endl;
  return (fabs(mu->userFloat("d0")) < maxMuonDxy);
}//--- PassMuonDxyCut

////////////////////////////////
/////Check jet Properties //////
////////////////////////////////
bool HadronicVZAnalyzer::PassJetCut(const pat::Jet* jet){
  bool jetKinStatus = PassJetPtCut(jet) && PassJetEtaCut(jet);
  bool jetIDStatus = PassJetIDCut(jet);
  return jetKinStatus && jetIDStatus;
}

bool HadronicVZAnalyzer::PassJetPtCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet Pt Cut"<<endl;
  return (jet->pt() > minJetPt);
}

bool HadronicVZAnalyzer::PassJetEtaCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet Eta Cut"<<endl;
  return (fabs(jet->eta()) < maxJetEta);
}

bool HadronicVZAnalyzer::PassJetNHFCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet NHF Cut"<<endl;
  return (jet->neutralHadronEnergyFraction() < maxJetNHF);
}

bool HadronicVZAnalyzer::PassJetNEFCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet NEF Cut"<<endl;
  return (jet->neutralEmEnergyFraction() < maxJetNEF);
}

bool HadronicVZAnalyzer::PassJetNConstCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet NConst Cut"<<endl;
  return (jet->numberOfDaughters() > minJetnumConst);
}

bool HadronicVZAnalyzer::PassJetCHFCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet CHF Cut"<<endl;
  return (jet->chargedHadronEnergyFraction() > minJetCHF);
}

bool HadronicVZAnalyzer::PassJetCMultCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet CMult Cut"<<endl;
  return ((unsigned int)jet->chargedMultiplicity() > minJetcMult);
}

bool HadronicVZAnalyzer::PassJetCEFCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet CEF Cut"<<endl;
  return (jet->chargedEmEnergyFraction() < maxJetCEF);
}

bool HadronicVZAnalyzer::PassJetIDCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Pass JetID Cut"<<endl;
  bool neutralStatus = 
    PassJetNHFCut(jet) &&
    PassJetNEFCut(jet) &&
    PassJetNConstCut(jet);
  bool chargedStatus = (fabs(jet->eta()) > 2.4) ||
    (PassJetCHFCut(jet) &&
     PassJetNHFCut(jet) &&
     PassJetCMultCut(jet)
     );
  return (neutralStatus && chargedStatus);
}

////////////////////////////////
/////Check TeV Properties/////
////////////////////////////////

/*float
HadronicVZAnalyzer::Calc_MuonRelIso(const TeVMuon* mu){
  return (mu->isolationR03().emEt + mu->isolationR03().hadEt + mu->isolationR03().sumPt)
    / mu->pt();
}
*/
///////////////Utilities//////////////////

void
HadronicVZAnalyzer::ClearEvtVariables(){
  muons_.clear();
  jets_.clear();
  looseMuons_.clear();
  tightMuons_.clear();
  zCand_ = ZCandidate();
  wCand = WCandidate();
}

void 
HadronicVZAnalyzer::reportProgress(int eventNum) {
  if (eventNum % intOptions_["report"] == 0) {
    printf("\rWe've processed %i events so far...", eventNum);
    cout.flush();
    printf("\n");
  }
  printf("Event number: %i", ++eventNum);
}

/// Print to screen (like printf), but only if --verbose option is on
void HadronicVZAnalyzer::verbose(const char *string, ...)
{
  if (intOptions_["verbose"]) {
    va_list ap;
    va_start(ap, string);
    vprintf(string, ap);
    va_end(ap);
    cout << endl;
  }
}

//Cory: This should really be beginSample!!!
void HadronicVZAnalyzer::beginFile(std::vector<wprime::InputFile>::const_iterator fi){
  TFileDirectory dir = wprimeUtil_->getFileService()->mkdir(fi->samplename); 
  Declare_Histos(dir);
  //ResetCounters();
}

// operations to be done when closing input file 
// (e.g. print summary)
void HadronicVZAnalyzer::endFile(std::vector<wprime::InputFile>::const_iterator fi,
                         ofstream & out){
  //ScaleHistos();//Already scaled
  //printSummary(fi->samplename);  
  //deleteHistos();
  //listOfHists.clear();
}

void HadronicVZAnalyzer::endAnalysis(ofstream & out){
}


TeVMuon &
HadronicVZAnalyzer::FindMuon(reco::Candidate & p){
  for(uint i=0; i<muons_.size(); ++i){
    if(Match(muons_[i], p)) return muons_[i];
  }
  cout<<"Didn't find match for muon!!!, returning random one\n";
  return muons_[0];
}


bool
HadronicVZAnalyzer::Match(TeVMuon & p1, reco::Candidate & p2){
  float tolerance = 0.0001;
  if (p1.pdgId() == p2.pdgId() &&
      fabs(p1.eta() - p2.eta()) < tolerance &&
      fabs(p1.phi() - p2.phi()) < tolerance
    )
    return true;
  return false;
}


/*
inline bool HadronicVZAnalyzer::inEE(const TeVMuon& mu){
  return fabs(mu.eta()) >= 1.05;
}

inline float
  HadronicVZAnalyzer::MuonPU(const TeVMuon & m){
  return rhoFastJet_*effectiveMuonArea_[inEE(m)];
}
*/
/* FROM WZUtilities.h
struct highestPt {
  bool operator() (const reco::Candidate * a, const reco::Candidate * b){

    return a->pt() > b->pt();
  }
};
*/


