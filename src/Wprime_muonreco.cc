#include "UserCode/CMGWPrimeGroup/interface/Wprime_muonreco.h"

// system include files
#include <memory>
#include <cmath>

#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/CaloJet.h"


#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

#include <iostream>
#include <map>

using std::cout; using std::endl; using std::string;

using namespace Wprime_muonreco_histo;

/// Constructor
Wprime_muonreco::Wprime_muonreco(const edm::ParameterSet& iConfig):
  muonTag_(iConfig.getParameter<edm::InputTag> ("MuonTag")),
  metTag_(iConfig.getParameter<edm::InputTag> ("MetTag")),
  HLTTag_(iConfig.getParameter<edm::InputTag>( "HLTriggerResults" ) ),
  //  isoTag_(iConfig.getParameter<edm::InputTag> ("IsolationTag")),
  jetTag_(iConfig.getParameter<edm::InputTag> ("JetTag")),
  tkIsoMapTag_(iConfig.getParameter<edm::InputTag> ("TkIsoMapTag")),
  ecalIsoMapTag_(iConfig.getParameter<edm::InputTag> ("EcalIsoMapTag")),
  hcalIsoMapTag_(iConfig.getParameter<edm::InputTag> ("HcalIsoMapTag")),
  eJetMin_(iConfig.getParameter<double>("EtJetCut")),
  detmu_acceptance(iConfig.getParameter<double>("Detmu_acceptance")),
  muHLT_20x(iConfig.getParameter<string>("SingleMuHLT_20x")),
  muHLT_21x(iConfig.getParameter<string>("SingleMuHLT_21x")),
  muL1(iConfig.getParameter<string>("SingleMuL1"))
{
   //now do what ever initialization is needed

}


/// Destructor
Wprime_muonreco::~Wprime_muonreco()
{
 
}

// get the generator info, populate gen_muons, set genmu_acceptance flag
void Wprime_muonreco::getGenMuons(const edm::Event & iEvent)
{
  if(realData)return;

  edm::Handle<edm::HepMCProduct> genEvt;
  iEvent.getByLabel("source", genEvt);
  if(!genEvt.isValid())
    {
      cerr << " Did not find generator-level info using label source! " << endl;
      abort();
    }
  const HepMC::GenEvent * myGenEvent = genEvt->GetEvent();
  
  // Et of single generator-level muon in event (if applicable)
  float Et1mu=-999;

  for(HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin(); p != myGenEvent->particles_end(); ++p ) 
    { // loop over pythia particles
      if ( !( abs((*p)->pdg_id())==13 && (*p)->status()==1 ) )  continue;

      gen_muons.push_back(*(*p));
      float eta   =(*p)->momentum().eta();
      if(abs(eta) < detmu_acceptance)
	genmu_acceptance = true;

      float Et = (*p)->momentum().e() * sin((*p)->momentum().theta());
      hPtGen->Fill(Et);
      Et1mu = Et;
    } //loop over pythia particles

  if(gen_muons.size() == 1)
    hPtGenOne->Fill(Et1mu);

}

// get trigger info, update MuTrig/genMuTrig, set muL1/HLT_acceptance flag
void Wprime_muonreco::getTriggers(const edm::Event & iEvent)
{

  edm::Handle<edm::TriggerResults> hltresults;
  iEvent.getByLabel(HLTTag_,hltresults);
  
  if (! hltresults.isValid() ) 
    { 
      cerr << "  -- No HLTRESULTS" << endl; 
      abort(); // what is the proper way of handling this??? skip run?
    }
  
  // shouldn't this be done at begin-job instead???
  if(triggerNames.size() == 0)
    {
      init_trigger(hltresults);
      const edm::Provenance & prov = iEvent.getProvenance(hltresults.id());
      HLTversion = prov.releaseVersion();
      cout << " CMSSW release used to produce trigger decisions : "
	   << HLTversion << endl;
    }
  
  // loop over triggers: extract event counts for trigger efficiencies
  for (unsigned int itrig = 0; itrig != N_triggers; ++itrig)
    {
      string trigName = triggerNames.triggerName(itrig);
      bool accept = hltresults->accept(itrig);
      
      if(accept)
	{  // trigger <trigName> has fired in this event
	  It it;
	  if(genmu_acceptance)
	    {
	      it = genMuTrig.trigger_count.find(trigName);
	      if(it != genMuTrig.trigger_count.end())
		(genMuTrig.trigger_count[trigName])++; 
	    }
	  it = MuTrig.trigger_count.find(trigName);
	  if(it != MuTrig.trigger_count.end())
	    (MuTrig.trigger_count[trigName])++; 
	  
	  if( (is20x(HLTversion) && trigName == muHLT_20x) ||
	      (is21x(HLTversion) && trigName == muHLT_21x) )
	    muHLT_acceptance = true;

	  if(trigName == muL1)
	    muL1_acceptance = true;
	}  // trigger <trigName> has fired in this event
      
    }


  for(unsigned imc = 0; imc != gen_muons.size(); ++imc)
    { // loop over pythia muons
      
      if(genmu_acceptance)
	{
	  float etag = gen_muons[imc].momentum().eta();
	  float EtGen = gen_muons[imc].momentum().e() * 
	    sin(gen_muons[imc].momentum().theta());
	  
	  h_mcmu_pt->Fill(EtGen);
	  h_mcmu_eta->Fill(etag);
	  if(muHLT_acceptance)
	    {
	      h_mcmu_pt_hlt->Fill(EtGen);
	      h_mcmu_eta_hlt->Fill(etag);
	    }
	  if(muL1_acceptance)
	    {
	      h_mcmu_eta_l1->Fill(etag);
	    }
	}
      
    } // loop over pythia muons

}

// get Calo-MET, initialize MET
void Wprime_muonreco::getCaloMET(const edm::Event & iEvent)
{
  // Get the MET collection from the event
  edm::Handle<reco::CaloMETCollection> metCollection;
  iEvent.getByLabel(metTag_, metCollection);
  
  CaloMETCollection::const_iterator caloMET = metCollection->begin();
  met_x += caloMET->px();
  met_y += caloMET->py();
}

// get Jets
void Wprime_muonreco::getJets(const edm::Event & iEvent)
{
  // Get the Jet collection from the event
  edm::Handle<reco::CaloJetCollection> jetCollection;
  iEvent.getByLabel(jetTag_, jetCollection);

  unsigned nJetsAboveThres[nBinNJets_veto] = {0};

  TAxis * yaxis = hEthresNjet->GetYaxis();

  // # of jets above <eJetMin_> threshold
  unsigned njets = 0;
  for (CaloJetCollection::const_iterator jet = jetCollection->begin(); 
       jet!=jetCollection->end(); ++jet) { // loop over jets
  
    hJetEt->Fill(jet->et());
    if (jet->et() > eJetMin_) ++njets;

    for(unsigned ietj = 1; ietj <= nBinEtJets_veto; ++ietj)
      { // loop over jet-Et bins
	float Ethresh = yaxis->GetBinCenter(ietj);
	if (jet->et() > Ethresh)   
	  ++(nJetsAboveThres[ietj - 1]);   
      } // loop over jet-Et bins

  } // loop over jets

  for(unsigned ietj = 1; ietj <= nBinEtJets_veto; ++ietj)
    { // loop over jet-Et bins
      float Ethresh = yaxis->GetBinCenter(ietj);
      hEthresNjet->Fill(nJetsAboveThres[ietj - 1], Ethresh);
    } // loop over jet-Et bins
   
   hNAlljets->Fill(jetCollection->size());
   hNjetsGtEJetMin->Fill(njets);

}  

// get Isolation
void Wprime_muonreco::getIsolation(const edm::Event & iEvent)
{
  // take iso deposits for tracks 
  // (contains (eta,phi, pt) of tracks in R<X (1.0) around each muon)
  iEvent.getByLabel(tkIsoMapTag_, tkMapH);

  // take iso deposits for ECAL 
  // (contains (eta,phi, pt) of ecal in R<X (1.0) around each muon)
  iEvent.getByLabel(ecalIsoMapTag_, ecalMapH);

  // take iso deposits for HCAL 
  // (contains (eta,phi, pt) of hcal in R<X (1.0) around each muon)
  iEvent.getByLabel(hcalIsoMapTag_, hcalMapH);
}

// get TeV muons
void Wprime_muonreco::getTeVMuons(const edm::Event & iEvent)
{
  edm::Handle<reco::TrackToTrackMap> tevMapH_default;
  edm::Handle<reco::TrackToTrackMap> tevMapH_1stHit;
  edm::Handle<reco::TrackToTrackMap> tevMapH_picky;

  iEvent.getByLabel("tevMuons", "default", tevMapH_default);
  tevMap_default = tevMapH_default.product();

  iEvent.getByLabel("tevMuons", "firstHit", tevMapH_1stHit);
  tevMap_1stHit = tevMapH_1stHit.product();

  iEvent.getByLabel("tevMuons", "picky", tevMapH_picky);
  tevMap_picky = tevMapH_picky.product();
}


// do MC matching
void Wprime_muonreco::doMCmatching()
{      
  for(unsigned i = 0; i != N_muons; ++i) 
    { // loop over reco muons 
      MuonRef mu(muonCollection,i);

      bool matched = false;

      for(unsigned imc = 0; imc != gen_muons.size(); ++imc)
	{ // loop over pythia muons
	  float etag = gen_muons[imc].momentum().eta();

    	  float deta = etag - mu->eta();
	  float dphi = gen_muons[imc].momentum().phi()-mu->phi();
	  if(dphi > 2*M_PI) dphi -= 2*M_PI;
	  if(dphi > M_PI) dphi = 2*M_PI - dphi;
	  float dr = sqrt(dphi*dphi+deta*deta);
	  
	  if( dr < 0.15)
	    {
	      matched = true;
	      float EtGen = gen_muons[imc].momentum().e() * 
		sin(gen_muons[imc].momentum().theta());
	      hPtRecoOverPtGen->Fill(mu->pt()/EtGen) ;
	      hPtRecoVsPtGen ->Fill(mu->pt(),EtGen); 
	    }

	} // loop over reco muons
      if (!matched) {hPtMuUnMatched->Fill(mu->pt()); }

    } // loop over pythia muons

}

// get muons, update MET
void Wprime_muonreco::getMuons(const edm::Event & iEvent)
{
  // Get the Muon Track collection from the event
  iEvent.getByLabel(muonTag_, muonCollection);

  if(RECOversion == INVALID_RELEASE)
    { // need a better way of doing this - in beginRun???
      const edm::Provenance& prov=iEvent.getProvenance(muonCollection.id());
      RECOversion = prov.releaseVersion();
      cout << " CMSSW release used to produce muon RECO : "
	   << RECOversion << endl;
      
    }

  // # of reconstructed muons in event  
  N_muons = muonCollection->size();
  hNMu->Fill(N_muons);

  N_muons_tot += N_muons;
  if(1 == N_muons)
    {
      ++(MuTrig.Nev_1mu);
      if(genmu_acceptance)
	++(genMuTrig.Nev_1mu);
    }

  for(unsigned i = 0; i != N_muons; ++i) 
    { // loop over reco muons 
      MuonRef mu(muonCollection,i);
      met_x -= mu->px();
      met_y -= mu->py();
    } // loop over reco muons 
  met = sqrt(met_x*met_x +met_y*met_y);
  hMET->Fill(met);
  
  for(unsigned i = 0; i != N_muons; ++i) 
    { // loop over reco muons 
      
      MuonRef mu(muonCollection,i);
      // pt
      hPtMu->Fill(mu->pt());
      // eta
      hEtaMu->Fill(mu->eta());
      
      // acoplanarity
      Geom::Phi<double> deltaphi(mu->phi()-atan2(met_y,met_x));
      double acop = deltaphi.value();
      if (acop<0) acop = - acop;
      acop = M_PI - acop;
      hAcop->Fill(acop);
       
      // transverse mass
      double w_et = mu->pt() + met;
      double w_px = mu->px() + met_x;
      double w_py = mu->py() + met_y;
      double massT = w_et*w_et - w_px*w_px - w_py*w_py;
      massT = (massT>0) ? sqrt(massT) : 0;
      hTMass->Fill(massT);
            
    } // loop over reco muons
  
}


// ------------ method called to for each event  ------------
void
Wprime_muonreco::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  init_event();

  // should this be done only at begin-job/run instead???
  realData = iEvent.isRealData();

  getGenMuons(iEvent);

  ++(MuTrig.Nev); // # of processed events
  if(genmu_acceptance)
    // # of processed events with gen-muon within det.acceptance
    ++(genMuTrig.Nev);

  getTriggers(iEvent);

  getCaloMET(iEvent);

  getJets(iEvent);

  getIsolation(iEvent);
  
  getMuons(iEvent);

  if(!realData)
    doMCmatching();

}

// initialize histograms
void Wprime_muonreco::init_histograms()
{
  hPtGen = fs->make<TH1F>("ptGenMu","Pt gen mu",nBinPtMu,minPtMu,maxPtMu);
  hPtGenJetVeto = fs->make<TH1F>("ptGenMuJetVeto","Pt gen mu Jet Veto",nBinPtMu,minPtMu,maxPtMu);
  hPtGenOne = fs->make<TH1F>("ptGenOne","Pt gen one mc mu",nBinPtMu,minPtMu,maxPtMu);
  hPtGenOneMatch= fs->make<TH1F>("ptGenOneMatch","Pt gen one mu reco",nBinPtMu,minPtMu,maxPtMu);
  hPtMuOneMatch= fs->make<TH1F>("ptMuOneMatch","Pt mu one mu reco",nBinPtMu,minPtMu,maxPtMu);
  hPtMuOneMatchJetVeto= fs->make<TH1F>("ptMuOneMatchJetVeto","Pt mu one mu reco JetVeto",nBinPtMu,minPtMu,maxPtMu);
  hPtGenOnePtlt5Match= fs->make<TH1F>("ptGenOnePtlt5Match","Pt gen mu One reco mu Ptlt5 Matched",nBinPtMu,minPtMu,maxPtMu);
  hPtRecoOverPtGen = fs->make<TH1F>("ptRecoOverPtGen","Pt reco over Pt Gen",200,0,30);
  hPtRecoOverPtGenJetVeto = fs->make<TH1F>("ptRecoOverPtGenJetVeto","Pt reco over Pt Gen Jet Veto",200,0,30);
  h1PtGen1PtReco = fs->make<TH1F>("OneptRecoOnePtGen","(1/Pt reco - 1/ Pt Gen)*PtGen",400,-1,1);
  h1PtGen1PtRecoVsPtGen = fs->make<TH2F>("OneptRecoOnePtGenVsPtGen","(1/Pt reco - 1/ Pt Gen)*PtGen vs pt gen",nBinPtMu,minPtMu,maxPtMu,400,-1,1);
  hPtRecoVsPtGen = fs->make<TH2F>("ptRecoVsPtGen","Pt reco vs pt Gen",nBinPtMu,minPtMu,maxPtMu,nBinPtMu,minPtMu,maxPtMu);

  hPtMuUnMatched = fs->make<TH1F>("ptMuUnMatched","Pt mu unmatched",nBinPtMu,minPtMu,maxPtMu);
  hPtMuUnMatchedJetVeto = fs->make<TH1F>("ptMuUnMatchedJetVeto","Pt mu unmatched",nBinPtMu,minPtMu,maxPtMu);

  hNMu    = fs->make<TH1F>("NMu","Nb. muons in the event",10,0.,10.);
  hPtMu   = fs->make<TH1F>("ptMu","Pt mu",nBinPtMu,minPtMu,maxPtMu);
  hPtOneMu   = fs->make<TH1F>("ptOneMu","Pt mu",nBinPtMu,minPtMu,maxPtMu);
  hPtMaxMu   = fs->make<TH1F>("ptMaxMu","Pt mu",nBinPtMu,minPtMu,maxPtMu);
  hPtMuJetVeto   = fs->make<TH1F>("ptMuJetVeto","Pt mu JetVeto",nBinPtMu,minPtMu,maxPtMu);
  hPtOneMuJetVeto   = fs->make<TH1F>("ptOneMuJetVeto","Pt mu JetVeto",nBinPtMu,minPtMu,maxPtMu);
  hPtMaxMuJetVeto   = fs->make<TH1F>("ptMaxMuJetVeto","Pt mu JetVeto",nBinPtMu,minPtMu,maxPtMu);
  hEtaMu  = fs->make<TH1F>("etaMu","Eta mu",50,-2.5,2.5);
  hEtaOnePtlt5Match=fs->make<TH1F>("etaOnePtlt5Match","Eta One Pt lt 5 Match",50,-2.5,2.5);
  hMET    = fs->make<TH1F>("MET","Missing Transverse Energy (GeV)",nBinPtMu,minPtMu,maxPtMu);  
  hTMass  = fs->make<TH1F>("TMass","Rec. Transverse Mass (GeV)",3000,0,6000);
  hAcop   = fs->make<TH1F>("Acop","Mu-MET acoplanarity",100,0,6.28);
  hNAlljets  = fs->make<TH1F>("NAlljets","njets",nBinNJets, minNJets, maxNJets);
  hNjetsGtEJetMin  = fs->make<TH1F>("Njetsgt40 ","njets e gt 40",nBinNJets, minNJets, maxNJets);
  hEthresNjet = 
    fs->make<TH2F>("EthresNjet ","Ethr vs NJetsSurv", 
		   nBinNJets_veto, minNJets_veto, maxNJets_veto, 
		   nBinEtJets_veto, minEtJets_veto, maxEtJets_veto);
  hEthresNjet_norm =
    fs->make<TH2F>("EthresNjet_norm ","Norm Ethr vs NJetsSurv",
		   nBinNJets_veto, minNJets_veto, maxNJets_veto, 
		   nBinEtJets_veto, minEtJets_veto, maxEtJets_veto);
  hJetEt =  fs->make<TH1F>("jetsEt","jets et (GeV)",nBinEtJets,minEtJets,maxEtJets);

  hHitTrack =   fs->make<TH1F>("hHitTrack","number of hits in the tracker ",400,0,100);
  hHitMuon =   fs->make<TH1F>("hHitMuon","number of hits in the muon system ",400,0,100);
  hHitComb =   fs->make<TH1F>("hHitComb","number of combined hits  ",400,0,100);
  hHitCheckMu = fs->make<TH1F>("hHitChekMu","number of combined hits - number of tracker hits ",400,0,100);

  hPtSumR03  = fs->make<TH1F>("ptSumR03","Sum R03 pT (GeV)",100,0.,100.);
  hPtSumNR03 = fs->make<TH1F>("ptSumNR03","Sum pT R03/pT",100,0.,50.);
  hTMass_PtSumR03 = fs->make<TH2F>("TMass_ptSumR03","Rec. Transverse Mass (GeV) vs Sum pT R03 (GeV)",100,0.,50.,150,0.,300.);

  hPtTkIso_ConeSize = fs->make<TH2F>("TkptSum_conesize","PtSum tk (GeV) vs cone size",nBinCone,minCone,maxCone,nBinSumPt,minSumPt,maxSumPt);
  hPtTkIsoThresh_ConeSize = fs->make<TH2F>("TkptSumThresh_conesize","PtSum tk Thresh (GeV) vs cone size",nBinCone,minCone,maxCone,nBinSumPt,minSumPt,maxSumPt);
  hPtTkIsoThresh_ConeSize_norm = fs->make<TH2F>("TkptSumThresh_conesize_norm","PtSum tk Thresh Norm (GeV) vs cone size",nBinCone,minCone,maxCone,nBinSumPt,minSumPt,maxSumPt);
  hPtTkIso_mmupt_ConeSize = fs->make<TH2F>("TkptSum_mmupt_conesize","PtSum tk - mu pt (GeV) vs cone size",nBinCone,minCone,maxCone,nBinSumPt,minSumPt,maxSumPt);
  hNTkIso_ConeSize = fs->make<TH2F>("NTk_conesize","Ntk vs cone size",nBinCone,minCone,maxCone,50,0.,50.);
  hEcalIso_ConeSize = fs->make<TH2F>("EcalEtSum_conesize","Ecal EtSum (GeV) vs cone size",nBinCone,minCone,maxCone,nBinSumPt,minSumPt,maxSumPt);
  hHcalIso_ConeSize = fs->make<TH2F>("HcalEtsum_conesize","Hcal EtSum (GeV) vs cone size",nBinCone,minCone,maxCone,nBinSumPt,minSumPt,maxSumPt);

  hMuChi2 = fs->make<TH1F>("muChi2","muChi2",100,0.,50.);
  hTrackerMuChi2= fs->make<TH1F>("muTrackerChi2","muTrackerChi2",100,0.,50.);
  hSAMuChi2= fs->make<TH1F>("muSAChi2","muSAChi2",100,0.,50.);

  hMuChi2UnMatched = fs->make<TH1F>("muChi2UnMatched","muChi2UnMatched",100,0.,50.);
  hTrackerMuChi2UnMatched= fs->make<TH1F>("muTrackerChi2UnMatched","muTrackerChi2UnMatched",100,0.,50.);
  hSAMuChi2UnMatched= fs->make<TH1F>("muSAChi2UnMatched","muSAChi2UnMatched",100,0.,50.);

  hMuNdf = fs->make<TH1F>("muNdf","muNdf",200,0.,100.);
  hTrackerMuNdf= fs->make<TH1F>("muTrackerNdf","muTrackerNdf",200,0.,100.);
  hSAMuNdf= fs->make<TH1F>("muSANdf","muSANdf",200,0.,100.);

  hMuNdfUnMatched = fs->make<TH1F>("muNdfUnMatched","muNdfUnMatched",200,0.,100.);
  hTrackerMuNdfUnMatched= fs->make<TH1F>("muTrackerNdfUnMatched","muTrackerNdfUnMatched",200,0.,100.);
  hSAMuNdfUnMatched= fs->make<TH1F>("muSANdfUnMatched","muSANdfUnMatched",200,0.,100.);
  
  string TeVName[4]={"TeVstd","TeVfirstHit","TeVpicky","TeVcocktail"};
  for (int u=0; u <4; u++){

    string hname ="ptMu"; hname += TeVName[u];
    string htit ="Pt mu "; htit += TeVName[u];
    hPtTevMu[u]=fs->make<TH1F>(hname.c_str(),htit.c_str(),nBinPtMu,minPtMu,maxPtMu);

    hname ="ptMuJetVeto"; hname += TeVName[u];
    htit ="Pt mu Jet Veto"; htit += TeVName[u];
    hPtTevMuJetVeto[u]=fs->make<TH1F>(hname.c_str(),htit.c_str(),nBinPtMu,minPtMu,maxPtMu);

    hname = "ptGlobalRecoOverPtTevReco_"   ; hname += TeVName[u];
    htit  = "Global muon Pt over Tev Muon pt " ; htit  += TeVName[u];
    hPtTevMuOverPtMu[u] = fs->make<TH1F>(hname.c_str(),htit.c_str(),500,0,10);

    hname = "OneptRecoOnePtGen_"   ; hname += TeVName[u];
    htit  = "(1/Pt reco - 1/ Pt Gen)*PtGen " ; htit  += TeVName[u];
    h1PtGen1PtRecoTevMu[u]=fs->make<TH1F>(hname.c_str(),htit.c_str(),400,-1,1);

    hname = "ptRecoOverPtGen_"   ; hname += TeVName[u];
    htit  = "Pt reco over Pt Gen " ; htit  += TeVName[u];
    hPtTevMuOverPtGen[u] = fs->make<TH1F>(hname.c_str(),htit.c_str(),500,0,10);
    
    hname = "OneptRecoOnePtGenVsPtGen_"   ; hname += TeVName[u];
    htit  = "(1/Pt reco - 1/ Pt Gen)*PtGen vs pt gen " ; htit  += TeVName[u];
    h1PtGen1PtRecoVsPtGenTevMu[u] = fs->make<TH2F>(hname.c_str(),htit.c_str(),nBinPtMu,minPtMu,maxPtMu,400,-1,1);
  }




  h_mcmu_pt = new TH1F("mcmu_pt", "Gen muons pt distribution", 100, 0, 500.0);
  h_mcmu_pt_hlt = new TH1F("mcmu_pt_hlt", "Gen muons pt distribution - HLT mu accepted events", 100, 0, 500.0);
  h_mcmu_eta = new TH1F("mcmu_eta", "Gen muons eta distribution",100,-2.5,2.5);
  h_mcmu_eta_hlt = new TH1F("mcmu_eta_hlt", "Gen muons eta distribution - HLT mu accepted events",100,-2.5,2.5);

  h_mcmu_eta_l1 = new TH1F("mcmu_eta_l1", "Gen muons eta distribution - L1 mu accepted events", 100, -2.5, 2.5);
}

const string Wprime_muonreco::INVALID_RELEASE = "invalid release number";

bool Wprime_muonreco::is21x(const string & release_string)
{
  return (release_string.find("CMSSW_2_1_") != string::npos);
}

bool Wprime_muonreco::is20x(const string & release_string)
{
  return (release_string.find("CMSSW_2_0_") != string::npos);
}

// initialize job info
void Wprime_muonreco::init_run()
{
  N_triggers = 0; realData = false;
  HLTversion = RECOversion = INVALID_RELEASE;
  init_histograms();
}

// initialize event info
void Wprime_muonreco::init_event()
{
  gen_muons.clear(); N_muons = N_muons_tot = 0;
  genmu_acceptance = muL1_acceptance = muHLT_acceptance = false;
  met_x = met_y = met = 0.0;
}


// ------------ method called once each job just before starting event loop  ----
void Wprime_muonreco::beginJob(const edm::EventSetup&)
{
  init_run();
}

// ------------ method called once each job just after ending the event loop  ------------
void Wprime_muonreco::endJob() {

  for (unsigned ietj = 0; ietj != nBinEtJets; ++ietj)
    {
      for(unsigned jnj = 0; jnj != nBinNJets; ++jnj)
	{
	  double this_eff= hEthresNjet->GetBinContent(ietj,jnj)/MuTrig.Nev;
	  hEthresNjet_norm->SetBinContent(ietj,jnj,this_eff);
	}
    }
  
  for(unsigned i = 1; i != 24; ++i) 
    {
      for (unsigned isumptt = 0; isumptt != 600; ++isumptt)
	{
	  double this_eff =hPtTkIsoThresh_ConeSize ->GetBinContent(i,isumptt)/MuTrig.Nev_1mu;
	  hPtTkIsoThresh_ConeSize_norm->SetBinContent(i,isumptt,this_eff);
	}
    }

    printSummary();
}

// print summary info over full job
void Wprime_muonreco::printSummary() const
{
  cout << " Processed " << MuTrig.Nev << " events" << endl;
  cout << " " << genMuTrig.Nev << " of those contained at least one MC muon"
       << " within |eta| < " << detmu_acceptance  << endl;
  
  std::ostringstream description; 
  description << " Trigger efficiencies for events\n w/ gen-muon within det. acceptance (" << genMuTrig.Nev << " events processed) :";
  if(!realData)
    printSummary2(genMuTrig, description.str());

  description.str(""); // interesting way to clear the description content...
  description << " Trigger efficiencies for all events (" << MuTrig.Nev 
		<< " events processed) :";
  printSummary2(MuTrig, description.str());
}

// print summary info for real
void Wprime_muonreco::printSummary2(const trigEff & trig, 
				    const string & description) const
{
  string longline("=============================================================");
  cout << longline << endl;
  cout<< description << endl << endl;
  
  for(It it = trig.trigger_count.begin(); it != trig.trigger_count.end();
      ++it)
    {
      float eff = 1.*it->second/(trig.Nev);
      float d_eff = sqrt(eff*(1-eff)/trig.Nev);
      cout << " TriggerName " << it->first;
      cout << ", Accepted events: " << it->second
	   << " (" << 100*eff << " +- " << 100*d_eff << " %) " << endl;
    }
  cout << endl;
}

// initialize trigger structure
void Wprime_muonreco::init_trigger(const edm::Handle<edm::TriggerResults> & 
				   hltresults)
{
  N_triggers = hltresults->size();
  triggerNames.init(*hltresults);

  for (unsigned int itrig = 0; itrig != N_triggers; ++itrig)
    {
      string trigName = triggerNames.triggerName(itrig);
      // the trigger class name (i.e. Muon) should be release-dependent ???
      if(trigName.find("Mu") != string::npos) // in 21x
	//      if(trigName.find("HLT1Muon") != string::npos) // in 20x
	{
	  if(!realData)genMuTrig.trigger_count[trigName] = 0;
	  MuTrig.trigger_count[trigName] = 0;
	}
    }

}


//define this as a plug-in
DEFINE_FWK_MODULE(Wprime_muonreco);
