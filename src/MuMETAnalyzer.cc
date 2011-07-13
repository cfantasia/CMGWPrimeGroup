#include "UserCode/CMGWPrimeGroup/interface/MuMETAnalyzer.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"

#include <TH1F.h>
#include <TLorentzVector.h>

using std::string; using std::cout; using std::endl;

MuMETAnalyzer::MuMETAnalyzer(const edm::ParameterSet& cfg,WPrimeUtil * wprimeUtil)
{
  wprimeUtil_ = wprimeUtil;
  assert(wprimeUtil_);

  muons_       = cfg.getParameter<edm::InputTag>("muons"  );
  met_       = cfg.getParameter<edm::InputTag>("met"  );
  particleFlow_ = cfg.getParameter<edm::InputTag>("particleFlow" );
  muReconstructor_   = cfg.getParameter<int>("muonReconstructor");
  muonPtThreshold_   = cfg.getParameter<double>("muonPtThreshold");
  chi2Cut_           = cfg.getParameter<double>("chi2Cut");
  muonEtaCut_        = cfg.getParameter<double>("muonEtaCut");
  oneMuPtTrackCut_   = cfg.getParameter<double>("oneMuPtTrackCut");
  relIsoCut_        = cfg.getParameter<double>("relIsoCut");
  highestPtMuonOnly_ = cfg.getParameter<bool>("highestPtMuonOnly");
  dumpHighPtMuons_   = cfg.getParameter<bool>("dumpHighPtMuons");
  dumpHighPtMuonThreshold_ = cfg.getParameter<double>("dumpHighPtMuonThreshold");

  assert(muReconstructor_ < Num_MuTeVtrkAlgos);

  setupCutOrder();

  reconstructors.push_back(kGLOBAL);
  reconstructors.push_back(kINNER);
  reconstructors.push_back(kCOCKTAIL);
  reconstructors.push_back(kDYT);
 
}

MuMETAnalyzer::~MuMETAnalyzer()
{
}

void MuMETAnalyzer::defineHistos(TFileDirectory & dir)
{
  for(int i = 0; i != Num_mumet_cuts; ++i)
    hPT[i] = hETA[i] = hPHI[i] = /*hMJDPHI[i] =*/ hTM[i] = 0;

  defineHistos_MuonPt(dir);
  defineHistos_MuonEta(dir);
  defineHistos_MuonPhi(dir);
  //  defineHistos_MuonJetDPhi(dir);
  defineHistos_MuonIso(dir);
  defineHistos_TMass(dir);
}

// Get the hardest muon (based on tracker-pt) in event
// (returns index in pat::MuonCollection)
int MuMETAnalyzer::getTheHardestMuon()
{
  int nmuons = muons->size();
  float temp_muPT = -999;
  int ret = -1;
  for(int j = 0; j != nmuons; ++j)
    {
      TeVMuon muon((*muons)[j], muReconstructor_);
      if(muon.innerTrack().isNull())continue;
      float current_muPT = muon.innerTrack()->pt();
      if (current_muPT > temp_muPT) 
	{
	  temp_muPT = current_muPT;
	  ret = j;
	}
    }
  return ret;
}


void MuMETAnalyzer::eventLoop(edm::EventBase const & event)
{
  event.getByLabel(muons_, muons);
  event.getByLabel(met_, met);

  // switch to help us keep track of whether a muon has already
  // been found in current event surviving the ith-cut;
  // this will ensure that we increase Num_surv_cut maximum once per evet
  // whereas we nevertheless fill the histograms 
  // for every muon surviving the i-th cut
  bool accountMe[Num_mumet_cuts];
  for(int cut = 0; cut != Num_mumet_cuts; ++cut)
    accountMe[cut] = true;

  unsigned iMuMin = 0; unsigned iMuMax = muons->size();
  if(highestPtMuonOnly_) // if true, will only consider highest-pt muon in event
    {
      int iMu = getTheHardestMuon();
      if(iMu >=0) // there is at least one muon in the event
	{
	  iMuMin = iMu;
	  iMuMax = iMu + 1;
	}
    }

  //loop over muons
  for (unsigned theMu = iMuMin; theMu != iMuMax; ++theMu){//loop over muons
    bool fill_entry = true; // if true, will histogram muon
    isInvalidMuon_ = false;
    pfMETwithoutMuCalculated_ = false;

    TeVMuon muon((*muons)[theMu], muReconstructor_);
    mu4D = muon.p4(isInvalidMuon_);
    if(isInvalidMuon_)continue;
    if(mu4D.Pt() < muonPtThreshold_) continue;
    for(int cut_index = 0; cut_index != Num_mumet_cuts; ++cut_index)
      { // loop over selection cuts
	string arg = mumet_cuts_desc_short[cut_index];
	// call to function [as implemented in setupCutOder]
	// don't get me started about the syntax!
	bool survived_cut=(this->*cuts[arg])(&fill_entry, &muon, event);
	if(!survived_cut)break; // skip rest of selection cuts
	
	if(fill_entry)
	  tabulateMe(cut_index, accountMe, event, &muon);
	
	if(dumpHighPtMuons_ && fill_entry 
	   && cut_index == Num_mumet_cuts-1
	   && wprimeUtil_->runningOnData() && 
	   muon.innerTrack()->pt() > dumpHighPtMuonThreshold_)
	  printHighPtMuon(event, muon);
      
      } // loop over selection cuts


  } // loop over muons

}


// fill histograms for muon if fill_entry=true; update book-keeping 
// (via private member: stats); make sure stats gets updated maximum once per event
void MuMETAnalyzer::tabulateMe(int cut_index, bool accountMe[], 
			       edm::EventBase const& event, const TeVMuon * muon)
{
  // if the accountMe switch is on, increase the number of events passing the cuts
  // and turn the switch off so we don't count more than once per event
  if(accountMe[cut_index])
    {
      wprime::SampleStat::iterator it;
      it = stats.find(wprimeUtil_->getSampleName());
      if(it == stats.end())abort();
      ++((it->second)[cut_index].Nsurv_evt_cut);
      accountMe[cut_index] = false;
    }
  float weight = wprimeUtil_->getWeight();
  // fill the histograms
  hPT[cut_index]->Fill(mu4D.Pt(), weight);
  hETA[cut_index]->Fill(mu4D.Eta(), weight);
  hPHI[cut_index]->Fill(mu4D.Phi(), weight);
  //  hMJDPHI[cut_index]->Fill( XJetDPhi(mu4D, event), weight);
  hTM[cut_index]->Fill(WPrimeUtil::TMass(mu4D, getNewMET(event, mu4D)), weight);
  hISO[cut_index]->Fill(muon->trkRelIsolation(),weight);
}

// operations to be done when changing input file (e.g. create new histograms)
void MuMETAnalyzer::beginFile(std::vector<wprime::InputFile>::const_iterator fi)
{
  wprime::FilterEff tmp; 
  for(int cut = 0; cut != Num_mumet_cuts; ++cut)
    stats[fi->samplename].push_back(tmp);

  // add channel/analysis name here?
  TFileDirectory dir= wprimeUtil_->getFileService()->mkdir(fi->samplename); 
  defineHistos(dir); // one set of histograms per input file

}

// operations to be done when closing input file 
// (e.g. print summary)
void MuMETAnalyzer::endFile(std::vector<wprime::InputFile>::const_iterator fi,
			    ofstream & out)
{
  printFileSummary(fi, out);
}

// print summary of efficiencies
void MuMETAnalyzer::printFileSummary(std::vector<wprime::InputFile>::const_iterator fi, ofstream & out)
{
  string sample = fi->samplename;
  float weight = fi->weight;
  wprime::SampleStat::iterator it;
  it = stats.find(sample);
  if(it == stats.end())abort();  

  float eff, deff;
  out << "\n Sample: " << sample << endl;
  out << " Total # of produced events for " << wprimeUtil_->getLumi_ipb() 
       << " ipb = " << fi->Nprod_evt*weight << endl;
  out << " Total # of events after pre-selection for " 
       << wprimeUtil_->getLumi_ipb() << " ipb = " << fi->Nact_evt*weight 
       << endl;

  WPrimeUtil::getEff(eff, deff, fi->Nact_evt, fi->Nprod_evt);
  out << " Preselection efficiency = " << eff << " +- " << deff << endl;

  for(int cut_index = 0; cut_index != Num_mumet_cuts; ++cut_index){
    
    (it->second)[cut_index].Nsurv_evt_cut_w = 
      (it->second)[cut_index].Nsurv_evt_cut*weight;
    
    out << " Cut # " << cut_index << ": " << mumet_cuts_desc_long[cut_index] 
	 <<", expected # of evts = " 
	 << (it->second)[cut_index].Nsurv_evt_cut_w;
    
    //calculate efficiencies
    if(cut_index == 0)
      WPrimeUtil::getEff(eff, deff, (it->second)[cut_index].Nsurv_evt_cut,
	     fi->Nact_evt);
    else
      WPrimeUtil::getEff(eff, deff, (it->second)[cut_index].Nsurv_evt_cut,
	     (it->second)[cut_index-1].Nsurv_evt_cut);
    out << ", Relative eff = "<<eff << " +- " << deff;
    WPrimeUtil::getEff(eff, deff, (it->second)[cut_index].Nsurv_evt_cut, 
	   fi->Nprod_evt);
    out << ", Absolute eff = "<< eff << " +- " << deff
	 << endl;
    
    (it->second)[cut_index].eff = eff;
    (it->second)[cut_index].deff = deff;
    
  } // loop over different cuts
	
}

// e.g. print summmary of expected events for all samples
void MuMETAnalyzer::endAnalysis(ofstream & out)
{
  float N_SM = 0; 
  std::map<std::string, std::vector<wprime::FilterEff> >::const_iterator it;

  int index = Num_mumet_cuts-1; // get # of events after last cut

  out << endl;
  for(it = stats.begin(); it != stats.end(); ++it)
    { // loop over samples
      string sample = it->first;
      float N_evt = (it->second)[index].Nsurv_evt_cut_w;
      out<< " "<< sample << ": " << N_evt
	 << " evts (eff = " << 100.*(it->second)[index].eff
	 << " +- " << 100.*(it->second)[index].deff
	 << " %) " << endl;
      
      if(!wprimeUtil_->runningOnData() &&
	 sample.find("wprime") == string::npos)
	N_SM += N_evt;
      
    } // loop over samples

  out << " Total # of SM events: " 
      << N_SM << endl;
  
}


void MuMETAnalyzer::defineHistos_MuonPt(TFileDirectory & dir)
{
  for(int cut = 0; cut != Num_mumet_cuts; ++cut)
    {
      string name = "hPT" + algo_desc_short[muReconstructor_] + "_" + 
	mumet_cuts_desc_short[cut];
      string title = algo_desc_long[muReconstructor_]+ " muon p_{T} with " + 
	mumet_cuts_desc_long[cut];
      
      hPT[cut] = dir.make<TH1F>(name.c_str(), title.c_str(), nBinPtMu,
			    minPtMu, maxPtMu);
    }
}

void MuMETAnalyzer::defineHistos_MuonEta(TFileDirectory & dir)
{
  for(int cut = 0; cut != Num_mumet_cuts; ++cut)
    {
      string name = "hETA" + algo_desc_short[muReconstructor_] + "_" 
	+ mumet_cuts_desc_short[cut];
      string title = algo_desc_long[muReconstructor_] + " muon #eta with " 
	+ mumet_cuts_desc_long[cut];
      hETA[cut] = dir.make<TH1F>(name.c_str(), title.c_str(), nBinEtaMu,
			   minEtaMu,maxEtaMu);
    }

}

void MuMETAnalyzer::defineHistos_MuonPhi(TFileDirectory & dir)
{
  for(int cut = 0; cut != Num_mumet_cuts; ++cut)
    {
      string name = "hPHI" + algo_desc_short[muReconstructor_] + "_" + 
	mumet_cuts_desc_short[cut];
      string title = algo_desc_long[muReconstructor_] + " muon #phi with " + 
	mumet_cuts_desc_long[cut];
      hPHI[cut] = dir.make<TH1F>(name.c_str(), title.c_str(), nBinPhiMu,
				 minPhiMu,maxPhiMu);
    }

}

#if 0
void MuMETAnalyzer::defineHistos_MuonJetDPhi(TFileDirectory & dir)
{
  for(int cut = 0; cut != Num_mumet_cuts; ++cut)
    {
      string name = "hMJDPHI" + algo_desc_short[muReconstructor_] + "_" + 
	mumet_cuts_desc_short[cut];
      string title = algo_desc_long[muReconstructor_] + " muon-jet #Delta#phi with " + 
	mumet_cuts_desc_long[cut];
      hMJDPHI[cut] = dir.make<TH1F>(name.c_str(), title.c_str(), 
			      nBinDPhiMu,minDPhiMu,maxDPhiMu);
    }

}
#endif

void MuMETAnalyzer::defineHistos_MuonIso(TFileDirectory & dir)
{
  for(int cut = 0; cut != Num_mumet_cuts; ++cut)
    {
      string name = "hISO" + algo_desc_short[muReconstructor_] + "_" + 
	mumet_cuts_desc_short[cut];
      string title = algo_desc_long[muReconstructor_] + " muon isol with " + 
	mumet_cuts_desc_long[cut];
      hISO[cut] = dir.make<TH1F>(name.c_str(), title.c_str(), nBinIsoMu,
			   minIsoMu,maxIsoMu);
    }

}

void MuMETAnalyzer::defineHistos_TMass(TFileDirectory & dir)
{
  for(int cut = 0; cut != Num_mumet_cuts; ++cut)
    {
      string name = "hTM" + algo_desc_short[muReconstructor_] + "_" 
	+ mumet_cuts_desc_short[cut];
      string title = algo_desc_long[muReconstructor_]+ " Transv. Mass with " 
	+ mumet_cuts_desc_long[cut];
      hTM[cut] = dir.make<TH1F>(name.c_str(), title.c_str(), nBinTmMu,
			  minTmMu,maxTmMu);
    }
  
}

void MuMETAnalyzer::setupCutOrder()
{
  cuts.clear();
#if debugmeMuMet
  cout << "\n Mu+MET cuts will be applied in this order: " << endl;
#endif

  for(int cut_i = 0; cut_i != Num_mumet_cuts; ++cut_i)
    { // loop over selection cuts
      string arg = mumet_cuts_desc_short[cut_i];
#if debugmeMuMet
      cout << " Cut #" << (cut_i+1) << ": " << mumet_cuts_desc_long[cut_i]
	   << " (" << arg << ") " << endl;
#endif
      if(arg == "hlt")cuts[arg] = &MuMETAnalyzer::passedHLT;
      else if(arg == "qual")cuts[arg] = &MuMETAnalyzer::goodQualityMuon;
      else if(arg == "1mu")cuts[arg] = &MuMETAnalyzer::onlyOneHighTrackPtMuon;
      else if(arg == "iso")cuts[arg] = &MuMETAnalyzer::isolatedMuon;
      else if(arg == "met")cuts[arg] = &MuMETAnalyzer::kinematicCuts;
      else
	{
	  cout << " Oops! Don't understand how to prepare for cut nicknamed as "
	       << arg << endl;
	  abort();
	}
    } // loop over selection cuts

#if debugmeMuMet
  cout << endl;
#endif

}

// dump on screen info about high-pt muon
void MuMETAnalyzer::printHighPtMuon(edm::EventBase const & event, TeVMuon & muon) 
{
  cout << "\n Run # = " << event.id().run() << " Event # = " 
       << event.id().event() << " LS = " << event.id().luminosityBlock() 
       << endl;
  cout << " Muon eta = " << mu4D.Eta() << "  phi = " << mu4D.Phi() << endl;
  pat::METCollection::const_iterator oldMET = met->begin();
  TVector2 oldMETv(oldMET->px(), oldMET->py());
  //  cout << " default pfMET = " << oldMET->pt() << " GeV ";

  typedef std::vector<unsigned>::iterator It;

  for(It it = reconstructors.begin(); it != reconstructors.end(); ++it)
    {
      bool isInvalid = true;
      unsigned rec_i(*it);
      TLorentzVector p4 = muon.p4(rec_i, isInvalid);
      if(isInvalid)continue;
      TVector2 newMET = getNewMET(event, p4);
      cout << " " << algo_desc_long[rec_i] << " pt = "
	   << muon.GetTrack(rec_i)->pt() << " +- " 
	   << muon.GetTrack(rec_i)->ptError()
	   << " GeV, charge = " << muon.GetTrack(rec_i)->charge() 
	   << ", TM = " << WPrimeUtil::TMass(p4, newMET) << " GeV " << endl;
    }
      

}

// Get new MET: there are two corrections to be made:
// (a) the hadronic MET component (that needs to be corrected 
// if applyMETCorrection=true)from Z data; this will be done according to hadronic 
// activity from Z->mumu reconstructed events
// (b) the muon-pt component that needs to be updated if we switch to one
// of the dedicated high-pt muon reconstructors
TVector2 MuMETAnalyzer::getNewMET(edm::EventBase const & event, const TLorentzVector & mu_p)
{
  if(wprimeUtil_->shouldApplyMETCorrection())
    return wprimeUtil_->getHadronicMET(event) - TVector2(mu_p.Px(), mu_p.Py());

  return getPFMETwithoutMu(event) - TVector2(mu_p.Px(), mu_p.Py());
}

// get (PF) MET without the default-pt for the running muon in event (mu4D);
// this is done so that we can adjust the muon-pt component of the MET by 
// switching to one of the dedicated high-pt muon reconstructors
TVector2 MuMETAnalyzer::getPFMETwithoutMu(edm::EventBase const & event)
{
  if(pfMETwithoutMuCalculated_)
    return pfMETwithoutMuCached_;

  pat::METCollection::const_iterator oldMET = met->begin();
  float met_x = oldMET->px(); float met_y = oldMET->py();

  //----Loop over PFCandidates to get the hardest muon and put its
  //pT back to the MET TVector2; create a new object
  edm::Handle<pat::PFParticleCollection > PFCandidates;
  event.getByLabel(particleFlow_,PFCandidates);
  int nmuon = 0;
  pat::PFParticleCollection::const_iterator iParticle;
  pat::PFParticleCollection::const_iterator icorrespondingPFMu;
  for( iParticle = (PFCandidates.product())->begin() ; 
       iParticle != (PFCandidates.product())->end() ; ++iParticle ){
    const reco::Candidate* candidate = &(*iParticle);
    if (candidate) {
      const reco::PFCandidate* pfCandidate = 
	dynamic_cast<const reco::PFCandidate*> (candidate);
      if (pfCandidate){
	//check that it is a muon
	if (!(pfCandidate->particleId() == 3)) continue;
	double deta = mu4D.Eta() - pfCandidate->eta();
	double dphi = reco::deltaPhi(mu4D.Phi(), pfCandidate->phi());
	double dR = TMath::Sqrt(deta*deta + dphi*dphi);
	if(dR < 0.01)
	  {
	    icorrespondingPFMu = iParticle;
	    ++nmuon;
	    break;
	  }
         
      }//pfCandidate
    }//if candidate
  }//loop over PFCandidates
  
  //Add back the pT for the hardest muon if found.
  if (nmuon > 0){
    const reco::Candidate* mycandidate = &(*icorrespondingPFMu);
    const reco::PFCandidate* mypfCandidate = 
      dynamic_cast<const reco::PFCandidate*> (mycandidate);
    met_x += mypfCandidate->px();
    met_y += mypfCandidate->py();
  }//if nmuon >0
  else { // if no match found, use mu4D value
    met_x += mu4D.Px();
    met_y += mu4D.Py();
  }
  pfMETwithoutMuCached_.Set(met_x,met_y);
  pfMETwithoutMuCalculated_ = true;
  return pfMETwithoutMuCached_;
}


// whether HLT accepted the event
bool MuMETAnalyzer::passedHLT(bool *, const TeVMuon *, edm::EventBase const &)
{
  // needs implementation
  return true;
}

// check if muon satisfies quality requirements
// fill goodQual; always returns true
bool MuMETAnalyzer::goodQualityMuon(bool * goodQual, const TeVMuon * muon, edm::EventBase const &)
{
  *goodQual = muon->goodQualityMuon(chi2Cut_, muonEtaCut_);
  return true;
}

// true if only one muon with track pt > the threshold
bool MuMETAnalyzer::onlyOneHighTrackPtMuon(bool *, const TeVMuon *, edm::EventBase const &)
{
  return (nMuAboveThresh(oneMuPtTrackCut_) == 1);
}

// returns # of (global) muons with tracker-pt above <tracker_muon_pt>
unsigned MuMETAnalyzer::nMuAboveThresh(float tracker_muon_pt)
{
  unsigned N = 0;
  int iMuMin = 0; int iMuMax = muons->size();
  //loop over muons
  for (int theMu = iMuMin; theMu != iMuMax; ++theMu){//loop over muons
    const pat::Muon * pMuon = &(*muons)[theMu];
    
    // consider only global muons
    if(!pMuon->isGood("AllGlobalMuons"))continue;
    
    reco::TrackRef trk = pMuon->innerTrack();
    if(trk.isNull())continue;
    if(trk->pt() > tracker_muon_pt)
      ++N;
  }
 
  return N;
}

// set bool flag to true if muon isolated
// always returns true
bool MuMETAnalyzer::isolatedMuon(bool * goodQual, const TeVMuon * muon, 
				 edm::EventBase const &)
{
  if(muon->trkRelIsolation() > relIsoCut_)
    *goodQual = false;

  return true;
}

// check if muon, MET pass kinematic cuts, updated goodQual
// always returns true
bool MuMETAnalyzer::kinematicCuts(bool * goodQual, const TeVMuon *, 
				  edm::EventBase const & event)
{
  TVector2 MET = getNewMET(event, mu4D);
  float ratio = mu4D.Pt()/MET.Mod();
  
  TVector2 muon_T(mu4D.Px(), mu4D.Py());
  float delta_phi = MET.DeltaPhi(muon_T);
  
  if(ratio < 0.4 || ratio > 1.5 || TMath::Abs(delta_phi) < 2.5)
    *goodQual = false;

  return true;
}
