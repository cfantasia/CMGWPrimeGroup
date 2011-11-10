#include "UserCode/CMGWPrimeGroup/interface/MuMETAnalyzer.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"

#include <TH1F.h>
#include <TH2F.h>

#include <stdio.h>

using std::string; using std::cout; using std::endl;

MuMETAnalyzer::MuMETAnalyzer(const edm::ParameterSet& cfg,int fileToRun) :
  AnalyzerBase(cfg, fileToRun){

  muonPtThreshold_   = cfg.getParameter<double>("muonPtThreshold");
  chi2Cut_           = cfg.getParameter<double>("chi2Cut");
  muonEtaCut_        = cfg.getParameter<double>("muonEtaCut");
  oneMuPtTrackCut_   = cfg.getParameter<double>("oneMuPtTrackCut");
  relIsoCut_        = cfg.getParameter<double>("relIsoCut");
  highestPtMuonOnly_ = cfg.getParameter<bool>("highestPtMuonOnly");
  dumpHighPtMuons_   = cfg.getParameter<bool>("dumpHighPtMuons");
  dumpHighPtMuonThreshold_ = cfg.getParameter<double>("dumpHighPtMuonThreshold");
  dumpHighMtMuonThreshold_ = cfg.getParameter<double>("dumpHighMtMuonThreshold");
  
  muMetRes.assign(NumSigSamples, NULL);
  setupCutOrder();

}

MuMETAnalyzer::~MuMETAnalyzer()
{
}

void MuMETAnalyzer::defineResolutionHistos(const TFileDirectory & dir, float Mass)
{
  int index = getSignalSampleIndex();
  createResolutionHist(dir, Mass, "mu", muMetRes[index]);
}


void MuMETAnalyzer::defineHistos(const TFileDirectory & dir)
{
  AnalyzerBase::defineHistos(dir);
  defineHistoset("hPT"+algo_desc_short[muReconstructor_], 
                 algo_desc_long[muReconstructor_]+ " muon p_{T} with ",
                 "p_{T} (GeV)", nBinPtMu, minPtMu, maxPtMu, "GeV", hPT,dir);
  defineHistoset("hPTgen"+algo_desc_short[muReconstructor_], 
                 algo_desc_long[muReconstructor_]+ " muon p_{T} with ",
                 "p_{T} (GeV)", nBinPtMu, minPtMu, maxPtMu, "GeV", hPTgen,dir);
  defineHistoset("hETA"+algo_desc_short[muReconstructor_], 
                 algo_desc_long[muReconstructor_]+ " muon #eta with ",
                 "#eta", nBinEtaMu, minEtaMu, maxEtaMu, "NONE", hETA,dir);
  defineHistoset("hPHI"+algo_desc_short[muReconstructor_], 
                 algo_desc_long[muReconstructor_]+ " muon #phi with ",
                 "#phi", nBinPhiMu, minPhiMu, maxPhiMu, "NONE", hPHI,dir);
  defineHistoset("hISO"+algo_desc_short[muReconstructor_], 
                 algo_desc_long[muReconstructor_]+ " muon isol with ",
                 "Iso", nBinIsoMu, minIsoMu, maxIsoMu, "NONE", hISO,dir);
  defineHistoset("hTM"+algo_desc_short[muReconstructor_], 
                 algo_desc_long[muReconstructor_]+ " Transv. Mass with ",
                 "m_{T} (GeV)", nBinTmMu, minTmMu, maxTmMu, "GeV", hTM,dir);
  defineHistos_TMvPT(dir);
}

// get the hardest muon (based on tracker-pt) in event
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
  event.getByLabel(muonsLabel_, muons);
  event.getByLabel(metLabel_, defMet);
  // switch to help us keep track of whether a muon has already
  // been found in current event surviving the ith-cut;
  // this will ensure that we increase Num_surv_cut maximum once per evet
  // whereas we nevertheless fill the histograms 
  // for every muon surviving the i-th cut
  bool* accountMe = new bool[NCuts_];
  for(int cut = 0; cut != NCuts_; ++cut)
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

    TeVMuon muon((*muons)[theMu], muReconstructor_);
    if(!muon.isValid())continue;
    if(muon.pt() < muonPtThreshold_) continue;
        
    Wcand = wprimeUtil_->getNewMETandW(event, muon, met, metLabel_);

    for(int cut_index = 0; cut_index != NCuts_; ++cut_index)
      { // loop over selection cuts
	
	string arg = CutNames_[cut_index];
	// call to function [as implemented in setupCutOder]
	// don't get me started about the syntax!
	bool survived_cut=(this->*cuts[arg])(&fill_entry, &muon, event);
	if(!survived_cut)break; // skip rest of selection cuts
	
	if(fill_entry)
	  tabulateMe(cut_index, accountMe, event, &muon);
	
	if(dumpHighPtMuons_ && fill_entry 
	   && cut_index == NCuts_-1
	   && wprimeUtil_->runningOnData() && 
	   ((muon.innerTrack()->pt() > dumpHighPtMuonThreshold_) ||
	    Wcand.mt() > dumpHighMtMuonThreshold_))
	  printHighPtMuon(event, (*muons)[theMu]);
      
      } // loop over selection cuts


  } // loop over muons
  delete[] accountMe;
}


// fill histograms for muon if fill_entry=true; update book-keeping 
// (via private member: stats); make sure stats gets updated maximum once per event
void MuMETAnalyzer::tabulateMe(int cut_index, bool accountMe[], 
			       edm::EventBase const& event, const TeVMuon * muon)
{
  // if the accountMe switch is on, increase the number of events passing the cuts
  // and turn the switch off so we don't count more than once per event
  float weight = wprimeUtil_->getWeight();
  if(accountMe[cut_index])
    {
      std::vector<wprime::InputFile>::iterator file = wprimeUtil_->getCurrentSample();
      file->results[cut_index].Nsurv_evt_cut_w += weight;
      file->results[cut_index].Nsurv_evt_cut++;
      accountMe[cut_index] = false;
    }
  // fill the histograms
  double genpt = -999;
  for(unsigned int i = 0 ; i != muon->genParticleRefs().size() ; ++i ){
    
    switch( muon->genParticle(i)->status() ){
    case 1:
      //		  histContainer_["DR_status1Match"]->Fill( ROOT::Math::VectorUtil::DeltaR(muon->p4() , muon->genParticle(i)->p4() ) ); 
      genpt = muon->genParticle(i)->pt();
      break;
    case 3:
      //		  histContainer_["DR_status3Match"]->Fill( ROOT::Math::VectorUtil::DeltaR(muon->p4() , muon->genParticle(i)->p4() ) ); 
      break;
    default:
      //		  histContainer_["DR_defaultMatch"]->Fill( ROOT::Math::VectorUtil::DeltaR(muon->p4() , muon->genParticle(i)->p4() ) ); 
      break;
    }
  }

  if(genpt > 0)
    hPTgen[cut_index]->Fill(genpt, weight);
  hPT[cut_index]->Fill(muon->pt(), weight);
  hETA[cut_index]->Fill(muon->eta(), weight);
  hPHI[cut_index]->Fill(muon->phi(), weight);
  hTM[cut_index]->Fill(Wcand.mt(), weight);
  hISO[cut_index]->Fill(muon->trkRelIsolation(),weight);
  hTMvPT[cut_index]->Fill(muon->pt(), Wcand.mt(), weight);

  if(wprimeUtil_->isSignalSample() && cut_index == NCuts_-1){
    float genMt = wprimeUtil_->getGenWprimeMt(event, PDG_ID_MUON, PDG_ID_MUONNEU, 
					      muon);
    int index = getSignalSampleIndex();
    muMetRes[index]->Fill(Wcand.mt() - genMt);
  }

}

void MuMETAnalyzer::defineHistos_TMvPT(const TFileDirectory & dir)
{
  for(int cut = 0; cut != NCuts_; ++cut)
    {
      string name = "hTMvPT" + algo_desc_short[muReconstructor_] + "_" 
	+ CutNames_[cut];
      string title = algo_desc_long[muReconstructor_]+ 
	" Transv. Mass vs p_{T} with " + CutDescs_[cut];
      hTMvPT[cut] = dir.make<TH2F>(name.c_str(), title.c_str(), 
				   nBinPtMu,minPtMu, maxPtMu,
				   nBinTmMu, minTmMu, maxTmMu);
    }
  
}

void MuMETAnalyzer::setupCutOrder()
{
  CutDescs_.resize(NCuts_);
  cuts.clear();
#if debugmeMuMet
  cout << "\n Mu+MET cuts will be applied in this order: " << endl;
#endif

  for(int cut_i = 0; cut_i != NCuts_; ++cut_i)
    { // loop over selection cuts
      string arg = CutNames_[cut_i];
      
      if(arg == "hlt")CutDescs_[cut_i] = "Single-muon HLT"; 
      else if(arg == "qual")CutDescs_[cut_i] = "Quality";
      else if(arg == "1mu")CutDescs_[cut_i] = "1 muon only";
      else if(arg == "iso")CutDescs_[cut_i] = "Isolation"; 
      else if(arg == "met")CutDescs_[cut_i] = "MET kinematic cuts";

#if debugmeMuMet
      cout << " Cut #" << (cut_i+1) << ": " << CutDescs_[cut_i]
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

void MuMETAnalyzer::printHighPtMuon(edm::EventBase const & event, const pat::Muon & muon) 
{
  cout << "\n Run # = " << event.id().run() << " Event # = " 
       << event.id().event() << " LS = " << event.id().luminosityBlock() 
       << endl;
  
  cout << " Muon eta = " << muon.eta() << "  phi = " << muon.phi() << endl;
  pat::METCollection::const_iterator oldMET = defMet->begin();
  TVector2 oldMETv(oldMET->px(), oldMET->py());
  cout << " default pfMET = " << oldMET->pt() << " GeV " << endl;

  typedef std::vector<unsigned>::iterator It;
  for(It it = muon_reconstructors.begin(); it != muon_reconstructors.end(); ++it)
    {
      TeVMuon mu(muon, *it);
      if(!mu.isValid())continue;

      pat::MET myMet;    
      WCandidate w = wprimeUtil_->getNewMETandW(event, mu, myMet, metLabel_);
 
      mu.printPtInfo(*it);
      cout << " MET = " << myMet.et() << " GeV, TM = " << w.mt() << " GeV " << endl;
    
      if(*it == muon_reconstructors.size()-1)
	mu.printTrackerInfo();
    }

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
bool MuMETAnalyzer::kinematicCuts(bool * goodQual, const TeVMuon * muon, 
				  edm::EventBase const & event)
{
  float ratio = muon->pt()/met.et();
  float delta_phi = Wcand.calcDPhi();

  if(ratio < 0.4 || ratio > 1.5 || TMath::Abs(delta_phi) < 2.5)
    *goodQual = false;

  return true;
}

void MuMETAnalyzer::fillHistos(const int& index, const float& weight){
}
