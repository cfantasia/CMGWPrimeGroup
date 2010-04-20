#include "UserCode/CMGWPrimeGroup/interface/Wprime_muonreco.h"
//

// system include files
#include <memory>
#include <cmath>

#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <fstream>
#include <map>

using std::cout; using std::endl; using std::string; using std::ifstream;

using namespace Wprime_muonreco_histo;

/// Constructor
Wprime_muonreco::Wprime_muonreco(const edm::ParameterSet& iConfig):
  muonTag_(iConfig.getParameter<edm::InputTag> ("MuonTag")),
  pfmetTag_(iConfig.getParameter<edm::InputTag> ("pfMetTag")),
  HLTTag_(iConfig.getParameter<edm::InputTag>( "HLTriggerResults" ) ),
  //  isoTag_(iConfig.getParameter<edm::InputTag> ("IsolationTag")),
  jetTag_(iConfig.getParameter<edm::InputTag> ("JetTag")),
  tkIsoMapTag_(iConfig.getParameter<edm::InputTag> ("TkIsoMapTag")),
  ecalIsoMapTag_(iConfig.getParameter<edm::InputTag> ("EcalIsoMapTag")),
  hcalIsoMapTag_(iConfig.getParameter<edm::InputTag> ("HcalIsoMapTag")),
  detmu_acceptance(iConfig.getParameter<double>("Detmu_acceptance")),
  muHLT_20x(iConfig.getParameter<string>("SingleMuHLT_20x")),
  muHLT_21x(iConfig.getParameter<string>("SingleMuHLT_21x")),
  muL1(iConfig.getParameter<string>("SingleMuL1")),
  sample_description(iConfig.getParameter<string>("description")),
  Nprod_evt(iConfig.getParameter<int>("Nprod_evt"))
{
   //now do what ever initialization is needed
  tree_job = tree_event = 0;

  evt = new wprime::Event(); job = new wprime::JobInfo();
  job->sample = sample_description;
  job->Nprod_evt = Nprod_evt;
  software_version = "V00-00-00";
}


/// Destructor
Wprime_muonreco::~Wprime_muonreco()
{
  if(evt) delete evt; if(job) delete job;
}

// get the generator info, populate gen_muons, set genmu_acceptance flag
void Wprime_muonreco::getGenParticles(const edm::Event & iEvent)
{
  if(realData)return;

  TClonesArray & mcmu = *(evt->mu_mc);
  TClonesArray & mcneu = *(evt->neu_mc);
  TClonesArray & mcw = *(evt->w_mc);
  TClonesArray & mcwp = *(evt->wp_mc);

  int Nm = 0; int Nn = 0; int Nw = 0; int Nwp = 0;

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);
  for(size_t i = 0; i != genParticles->size(); ++i) {
    const GenParticle & p = (*genParticles)[i];
    int id = p.pdgId(); 
    if( abs(id)!=13 && abs(id) != 14 && abs(id) != 24 && abs(id) != 34)
      continue; // keep only muons, neutrino, W or W'

    int st = p.status(); int ch = p.charge();
    int momid = -9999;

    TLorentzVector p4(p.px(), p.py(), p.pz(),  p.energy());
    const Candidate * mom = p.mother();
    if(mom)
      momid = mom->pdgId();

    switch(abs(id))
      {
      case 13: // muon
	new(mcmu[Nm]) wprime::MCParticle(p4, ch, momid, st);
	++Nm;
	
	if(st == 1)
	  {
	    gen_muons.push_back(p);
	    if(abs(p.eta())  < detmu_acceptance)
	      genmu_acceptance = true;
	    
	  }

	break;
	
      case 14: // muon-neutrino
	new(mcneu[Nn]) wprime::MCParticle(p4, ch, momid, st);
	++Nn;

	break;
	
      case 24: // W
	new(mcw[Nw]) wprime::MCParticle(p4, ch, momid, st);
	++Nw;
      
	break;
	
      case 34: // W'
	new(mcwp[Nwp]) wprime::MCParticle(p4, ch, momid, st);
	++Nwp;
	break;

      default:
	{
	  ; // do nothing for now
	}

      } // switch (abs(id))


  } // loop over genParticles

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

      job->HLTversion = HLTversion;
    }
  
  // loop over triggers: extract event counts for trigger efficiencies
  for (unsigned int itrig = 0; itrig != N_triggers; ++itrig)
    {
      string trigName = triggerNames.triggerName(itrig);
      bool accept = hltresults->accept(itrig);
      
      if(accept)
	{  // trigger <trigName> has fired in this event
	  tIt it;
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

	  if(trigName == "HLT_L1MuOpen")
	    evt->HLT_L1MuOpen = true;
	  else if (trigName == "HLT_L1Mu")
	    evt->HLT_L1Mu = true;
	  else if (trigName == "HLT_Mu3")
	    evt->HLT_Mu3 = true;
	  else if (trigName == "HLT_Mu5")
	    evt->HLT_Mu5 = true;
	  else if (trigName == "HLT_Mu9")
	    evt->HLT_Mu9 = true;

	  if(trigName == muL1)
	    muL1_acceptance = true;
	}  // trigger <trigName> has fired in this event
      
    }

}

// get Particle-Flow MET
void Wprime_muonreco::getPFMET(const edm::Event & iEvent)
{
  // Get the MET collection from the event
  edm::Handle<reco::PFMETCollection> pfmetCollection;
  iEvent.getByLabel(pfmetTag_, pfmetCollection);
  
  PFMETCollection::const_iterator pfMET = pfmetCollection->begin();
  met_x += pfMET->px();
  met_y += pfMET->py();
  evt->pfmet.Set(met_x, met_y);
}

// get Jets
void Wprime_muonreco::getJets(const edm::Event & iEvent)
{
  // Get the Jet collection from the event
  edm::Handle<reco::CaloJetCollection> jetCollection;
  iEvent.getByLabel(jetTag_, jetCollection);

  TClonesArray & jet = *(evt->jet);

  int j = 0;
  for (CaloJetCollection::const_iterator jet_ = jetCollection->begin(); 
       jet_ !=jetCollection->end(); ++jet_) { // loop over jets

    new(jet[j]) TLorentzVector(jet_->px(), jet_->py(), jet_->pz(), 
			       jet_->energy());
    ++j;
  } // loop over jets

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

// get muons
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
      
      job->RECOversion = RECOversion;
    }

  // # of reconstructed muons in event (including standalone-only)
  N_all_muons = muonCollection->size();
}

// copy tracking info from reco::Track to wprime::Track
void Wprime_muonreco::getTracking(wprime::Track & track, const reco::Track & p)
{
  TVector3 p3(p.px(), p.py(), p.pz());
  track.p.SetVectM(p3, wprime::MUON_MASS);
  track.q = p.charge();
  track.chi2 = p.chi2();
  track.dpt = p.ptError();
  track.dq_over_p = p.qoverpError();
  track.ndof = int(p.ndof());
  track.Ntot_hits = p.numberOfValidHits();
  track.Ntrk_hits = p.hitPattern().numberOfValidTrackerHits();
}

// do muon analysis
void Wprime_muonreco::doMuons()
{
  TClonesArray & recomu = *(evt->mu);

  for(unsigned i = 0; i != N_all_muons; ++i) 
    { // loop over reco muons 
      
      MuonRef mu(muonCollection, i);
      if(!(mu->isGlobalMuon()) )
	continue; // keep only global muons
      
      new(recomu[N_muons]) wprime::Muon(); 
      wprime::Muon * wpmu = (wprime::Muon *) recomu[N_muons];
      wpmu->Nmu_hits = mu->standAloneMuon()->recHitsSize(); 

      getTracking(wpmu->tracker, *(mu->track()));
      getTracking(wpmu->global, *(mu->combinedMuon()) );

      ++N_muons;
      
      doIsolation(mu, wpmu);

      doTeVanalysis(mu, wpmu);
      
    } // loop over reco muons

  if(1 == N_muons)
    {
      ++(MuTrig.Nev_1mu);
      if(genmu_acceptance)
	++(genMuTrig.Nev_1mu);
    
    }
  
}

// do TeV-muon analysis
void Wprime_muonreco::doTeVanalysis(reco::MuonRef mu, wprime::Muon * wpmu)
{
  if(!(mu->isGlobalMuon()) ) return; // keep only global muons

  TrackToTrackMap::const_iterator iTeV_default;
  TrackToTrackMap::const_iterator iTeV_1stHit;
  TrackToTrackMap::const_iterator iTeV_picky;

  if(is21x(RECOversion) || is22x(RECOversion) || is31x(RECOversion))
    {
      iTeV_default = tevMap_default->find(mu->globalTrack());
      iTeV_1stHit = tevMap_1stHit->find(mu->globalTrack());
      iTeV_picky = tevMap_picky->find(mu->globalTrack());
    }
  else if(is20x(RECOversion))
    {
      iTeV_default = tevMap_default->find(mu->combinedMuon());
      iTeV_1stHit = tevMap_1stHit->find(mu->combinedMuon());
      iTeV_picky = tevMap_picky->find(mu->combinedMuon());
    }
  else
    {
      cerr << " RECOversion " << RECOversion << " not supported. Sorry! " 
	   << endl;
            abort();
    }

  if(iTeV_default == tevMap_default->end() 
     || iTeV_1stHit == tevMap_1stHit->end() 
     || iTeV_picky == tevMap_picky->end())
    {
      cout << "-Wprime_muonreco- Warning: No Tev muons found for this event !! "
	   << endl; 
      return;
    }

  getTracking(wpmu->tev_1st, *(iTeV_1stHit->val) );

}



// do isolation
void Wprime_muonreco::doIsolation(MuonRef mu,  wprime::Muon * wpmu)
{
  // General Isolation
  const reco::IsoDeposit tkDep((*tkMapH)[mu]);
  const reco::IsoDeposit ecalDep((*ecalMapH)[mu]);
  const reco::IsoDeposit hcalDep((*hcalMapH)[mu]);

  for(unsigned i = 0; i != nBinCone; ++i) 
    { // loop over cone sizes
      float coneSize = minCone+i*(maxCone-minCone)/nBinCone;
      double ptsum_cone=tkDep.depositWithin(coneSize); 
      wpmu->SumPtIso[i] = ptsum_cone;

      double ntk_cone=tkDep.depositAndCountWithin(coneSize).second;
      wpmu->NtrkIso[i] = ntk_cone;

      double ecetsum_cone=ecalDep.depositWithin(coneSize); 
      wpmu->ECALIso[i] = ecetsum_cone;

      double hcetsum_cone=hcalDep.depositWithin(coneSize); 
      wpmu->HCALIso[i] = hcetsum_cone;
    } // loop over cone sizes

}

// ------------ method called to for each event  ------------
void
Wprime_muonreco::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  init_event();
  evt->evt_no = iEvent.id().event();
  evt->run_no = iEvent.id().run();


  // should this be done only at begin-job/run instead???
  realData = iEvent.isRealData();

  getGenParticles(iEvent);

  ++(MuTrig.Nev); // # of processed events
  if(genmu_acceptance)
    // # of processed events with gen-muon within det.acceptance
    ++(genMuTrig.Nev);

  getTriggers(iEvent);

  getPFMET(iEvent);

  getJets(iEvent);

  getIsolation(iEvent);

  getMuons(iEvent);
  getTeVMuons(iEvent);
  doMuons();

  tree_event->Fill();
}

// initialize histograms
void Wprime_muonreco::init_histograms()
{
  string common_desc = " - UserCode/CMGWPrimeGroup version " + software_version;
  string tree_title = "Job/file info" + common_desc;
  tree_job = fs->make<TTree>("jobinfo", tree_title.c_str());
  tree_title = "Wprime kinematic info per event" + common_desc;
  tree_event = fs->make<TTree>("wprime", tree_title.c_str());
  tree_job->Branch("job", "wprime::JobInfo", &job, 8000, 2);
  tree_event->Branch("wp", "wprime::Event", &evt, 8000, 2);
}

const string Wprime_muonreco::INVALID_RELEASE = "invalid release number";

bool Wprime_muonreco::is21x(const string & release_string)
{
  return (release_string.find("CMSSW_2_1_") != string::npos);
}

bool Wprime_muonreco::is31x(const string & release_string)
{
  return (release_string.find("CMSSW_3_1_") != string::npos);
}

bool Wprime_muonreco::is22x(const string & release_string)
{
  return (release_string.find("CMSSW_2_2_") != string::npos);
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
  gen_muons.clear();
  N_muons = N_all_muons = 0;
  genmu_acceptance = muL1_acceptance = muHLT_acceptance = false;
  met_x = met_y = met = 0.0;

  evt->mu_mc->Clear(); evt->neu_mc->Clear(); 
  evt->w_mc->Clear(); evt->wp_mc->Clear(); 
  evt->mu->Clear(); evt->jet->Clear();

  evt->HLT_L1MuOpen = evt->HLT_L1Mu = evt->HLT_Mu3 = evt->HLT_Mu5 =
    evt->HLT_Mu9 = false;
}


// ------------ method called once each job just before starting event loop  ----
void Wprime_muonreco::beginJob(const edm::EventSetup&)
{
  const char * _path = getenv("CMSSW_BASE");
  const char * _filename = "/src/UserCode/CMGWPrimeGroup/VERSION";

  if(_path == 0)
    {
      cerr << " *** Wprime_muonreco::beginJob *** " << endl;
      cerr << " *** Error! can't locate working directory ! ***\n";
      cerr << " (maybe the UserCode/CMGWPrimeGroup directory is missing?)\n"
	   << endl;
    }
  else
    {
      char verFile[1024];
      sprintf(verFile, "%s%s", _path, _filename);
      cout << " Found version file: " << verFile << endl;

      ifstream cin; cin.open(verFile);
      while(cin >> software_version)
	{
	  ; // write version # into string
	}
      cout << " Extracted version: " << software_version << endl;
    }

  init_run();

}

// ------------ method called once each job just after ending the event loop  ------------
void Wprime_muonreco::endJob() {

    printSummary();    
    tree_job->Fill();
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
  
  for(tIt it = trig.trigger_count.begin(); it != trig.trigger_count.end();
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
