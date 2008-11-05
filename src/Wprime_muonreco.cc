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
// using namespace edm; using namespace reco;
using std::cout; using std::endl; using std::string;

/// Constructor
Wprime_muonreco::Wprime_muonreco(const edm::ParameterSet& iConfig):
  muonTag_(iConfig.getParameter<edm::InputTag> ("MuonTag")),
  metTag_(iConfig.getParameter<edm::InputTag> ("MetTag")),
  HLTTag_(iConfig.getParameter<edm::InputTag>( "HLTriggerResults" ) ),
  //  isoTag_(iConfig.getParameter<edm::InputTag> ("IsolationTag")),
  jetTag_(iConfig.getParameter<edm::InputTag> ("JetTag")),
  eJetMin_(iConfig.getParameter<double>("EJetMin")),
  nBinNMu(iConfig.getUntrackedParameter<unsigned int>("NBinNMu")),
  minNMu(iConfig.getUntrackedParameter<double>("MinNMu")),
  maxNMu(iConfig.getUntrackedParameter<double>("MaxNMu")),
  nBinPtMu(iConfig.getUntrackedParameter<unsigned int>("NBinPtMu")),
  minPtMu(iConfig.getUntrackedParameter<double>("MinPtMu")),
  maxPtMu(iConfig.getUntrackedParameter<double>("MaxPtMu")),
  nBinEtaMu(iConfig.getUntrackedParameter<unsigned int>("NBinEtaMu")),
  minEtaMu(iConfig.getUntrackedParameter<double>("MinEtaMu")),
  maxEtaMu(iConfig.getUntrackedParameter<double>("MaxEtaMu")),
  nBinMET(iConfig.getUntrackedParameter<unsigned int>("NBinMET")),
  minMET(iConfig.getUntrackedParameter<double>("MinMET")),
  maxMET(iConfig.getUntrackedParameter<double>("MaxMET")),
  nBinTMass(iConfig.getUntrackedParameter<unsigned int>("NBinTMass")),
  minTMass(iConfig.getUntrackedParameter<double>("MinTMass")),
  maxTMass(iConfig.getUntrackedParameter<double>("MaxTMass")),
  nBinAcop(iConfig.getUntrackedParameter<unsigned int>("NBinAcop")),
  minAcop(iConfig.getUntrackedParameter<double>("MinAcop")),
  maxAcop(iConfig.getUntrackedParameter<double>("MaxAcop")),
  nBinNJets(iConfig.getUntrackedParameter<unsigned int>("NBinNJets")),
  minNJets(iConfig.getUntrackedParameter<double>("MinNJets")),
  maxNJets(iConfig.getUntrackedParameter<double>("MaxNJets")),
  detmu_acceptance(iConfig.getParameter<double>("Detmu_acceptance")),
  muHLT_20x(iConfig.getParameter<string>("SingleMuHLT_20x")),
  muHLT_21x(iConfig.getParameter<string>("SingleMuHLT_21x")),
  muL1(iConfig.getParameter<string>("SingleMuL1")),
  outputFileName(iConfig.getUntrackedParameter<string>("OutputFileName"))
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
  const HepMC::GenEvent * myGenEvent = genEvt->GetEvent();
  
  for(HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin(); p != myGenEvent->particles_end(); ++p ) 
    { // loop over pythia particles
      if ( !( abs((*p)->pdg_id())==13 && (*p)->status()==1 ) )  continue;

      gen_muons.push_back(*(*p));
      float eta   =(*p)->momentum().eta();
      if(abs(eta) < detmu_acceptance)
	genmu_acceptance = true;

      float Et = (*p)->momentum().e() * sin((*p)->momentum().theta());
      hPtGen->Fill(Et);
    } //loop over pythia particles

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
  
  CaloJetCollection::const_iterator jet = jetCollection->begin();
  unsigned njets = 0;
  for (jet=jetCollection->begin(); jet!=jetCollection->end(); ++jet) {
    if (jet->et()>eJetMin_) 
      ++njets;
  }
  hNjets->Fill(njets);
}  

// do MC matching
void Wprime_muonreco::doMCmatching()
{      
  for(unsigned i = 0; i != N_muons; ++i) 
    { // loop over reco muons 
      TrackRef mu(muonCollection,i);

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

  for(unsigned i = 0; i != N_muons; ++i) 
    { // loop over reco muons 
      TrackRef mu(muonCollection,i);
      met_x -= mu->px();
      met_y -= mu->py();
    } // loop over reco muons 
  met = sqrt(met_x*met_x +met_y*met_y);
  hMET->Fill(met);
  
  for(unsigned i = 0; i != N_muons; ++i) 
    { // loop over reco muons 
      
      TrackRef mu(muonCollection,i);
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
  
  getMuons(iEvent);

  if(!realData)
    doMCmatching();

}

// initialize histograms
void Wprime_muonreco::init_histograms()
{
  hPtGen = new TH1F("ptGenMu","Pt gen mu",nBinPtMu,minPtMu,maxPtMu);
  hPtRecoOverPtGen = new TH1F("ptRecoOverPtGen","Pt reco over Pt Gen",80,0,10);
  //  hPtRecoVsPtGen = new TH2F("ptRecoVsPtGen","Pt reco vs pt Gen",nBinPtMu,minPtMu,maxPtMu,nBinPtMu,minPtMu,maxPtMu);
  hPtRecoVsPtGen = new TH2F("ptRecoVsPtGen","Pt reco vs pt Gen",nBinPtMu,minPtMu, 800.0,nBinPtMu,minPtMu, 800.0);
  //  hPtMuUnMatched = new TH1F("ptMuUnMatched","Pt mu unmatched",nBinPtMu,minPtMu,maxPtMu);
  hPtMuUnMatched = new TH1F("ptMuUnMatched","Pt mu unmatched",80,minPtMu,200.0);

  // hNMu    = new TH1F("NMu","Nb. muons in the event",10,0.,10.);
  hNMu    = new TH1F("NMu","Nb. muons in the event",nBinNMu,minNMu,maxNMu);
  //  hPtMu   = new TH1F("ptMu","Pt mu",100,0.,100.);
  hPtMu   = new TH1F("ptMu","Pt mu",nBinPtMu,minPtMu,maxPtMu);
  //  hEtaMu  = new TH1F("etaMu","Eta mu",50,-2.5,2.5);
  hEtaMu  = new TH1F("etaMu","Eta mu",nBinEtaMu,minEtaMu,maxEtaMu);
  //  hMET    = new TH1F("MET","Missing Transverse Energy (GeV)", 100,0.,100.);
  hMET    = new TH1F("MET","Missing Transverse Energy (GeV)",nBinMET,minMET,maxMET);  
  //  hTMass  = new TH1F("TMass","Rec. Transverse Mass (GeV)",150,0.,300.);
  //  hTMass  = new TH1F("TMass","Rec. Transverse Mass (GeV)",100,0.,200.);
  hTMass  = new TH1F("TMass","Rec. Transverse Mass (GeV)",nBinTMass,minTMass,maxTMass);
  //  hAcop   = new TH1F("Acop","Mu-MET acoplanarity",50,0.,M_PI);
  //  hNjets  = new TH1F("Njets","njets",25,0.,25.);
  hAcop   = new TH1F("Acop","Mu-MET acoplanarity",nBinAcop,minAcop,maxAcop);
  hNjets  = new TH1F("Njets","njets",nBinNJets,minNJets,maxNJets);

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
  gen_muons.clear(); N_muons = 0;
  genmu_acceptance = muL1_acceptance = muHLT_acceptance = false;
  met_x = met_y = met = 0.0;
}


// ------------ method called once each job just before starting event loop  ----
void Wprime_muonreco::beginJob(const edm::EventSetup&)
{
  // Create output files

  outputRootFile = TFile::Open(outputFileName.c_str(), "RECREATE");
  outputRootFile->cd();

  init_run();
}

// ------------ method called once each job just after ending the event loop  ------------
void Wprime_muonreco::endJob() {
  // Write the histos to file
  outputRootFile->cd();


  // do we really need to explicitly call h_i->Write() for every i-th histogram???
  hPtGen->Write(); 
  hPtRecoOverPtGen->Write(); 
  hPtRecoVsPtGen->Write(); 
  hPtMuUnMatched->Write(); 

  hNMu->Write();
  hPtMu->Write();
  hEtaMu->Write();
  hMET->Write();
  hTMass->Write();
  hAcop->Write();
  hNjets->Write();

  h_mcmu_pt_hlt->Write(); h_mcmu_pt->Write();
  h_mcmu_eta_hlt->Write(); h_mcmu_eta_l1->Write(); h_mcmu_eta->Write();


  outputRootFile->Close();

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
