#include "UserCode/CMGWPrimeGroup/interface/EleMETAnalyzer.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"

#include <TH1F.h>
#include <TLorentzVector.h>

using std::string; using std::cout; using std::endl;

EleMETAnalyzer::EleMETAnalyzer(const edm::ParameterSet& cfg,WPrimeUtil * wprimeUtil) :
  cuts_(cfg)
{
  wprimeUtil_ = wprimeUtil;
  assert(wprimeUtil_);

  electrons_       = cfg.getParameter<edm::InputTag>("electrons"  );
  met_       = cfg.getParameter<edm::InputTag>("met"  );
  oneEleEtCut_   = cfg.getParameter<double>("oneEleEtCut");
  highestEtElectronOnly_ = cfg.getParameter<bool>("highestEtElectronOnly");
  dumpHighEtElectrons_   = cfg.getParameter<bool>("dumpHighEtElectrons");
  dumpHighEtElectronThreshold_ = cfg.getParameter<double>("dumpHighEtElectronThreshold");

  setupCutOrder();

}

EleMETAnalyzer::~EleMETAnalyzer()
{
}

void EleMETAnalyzer::defineHistos(TFileDirectory & dir)
{
  for(int i = 0; i != Num_elmet_cuts; ++i)
    hPT[i] = hETA[i] = hPHI[i] = /*hMJDPHI[i] =*/ hTM[i] = 0;

  defineHistos_ElectronEt(dir);
  defineHistos_ElectronEta(dir);
  defineHistos_ElectronPhi(dir);
  defineHistos_TMass(dir);
}

// Get the hardest muon (based on HEEP Et) in event
// (returns index in pat::ElectronCollection)
int EleMETAnalyzer::getTheHardestElectron()
{
  int nelectrons = electrons->size();
  float temp_elEt = -999;
  int ret = -1;
  for(int j = 0; j != nelectrons; ++j)
    {
      heep::Ele el((*electrons)[j]);
      float current_elEt = el.et();
      if (current_elEt > temp_elEt) 
	{
	  temp_elEt = current_elEt;
	  ret = j;
	}
    }
  return ret;
}


void EleMETAnalyzer::eventLoop(edm::EventBase const & event)
{
  event.getByLabel(electrons_, electrons);
  event.getByLabel(met_, met);

  // switch to help us keep track of whether a electron has already
  // been found in current event surviving the ith-cut;
  // this will ensure that we increase Num_surv_cut maximum once per evet
  // whereas we nevertheless fill the histograms 
  // for every electron surviving the i-th cut
  bool accountMe[Num_elmet_cuts];
  for(int cut = 0; cut != Num_elmet_cuts; ++cut)
    accountMe[cut] = true;

  int iEleMin = 0; int iEleMax = electrons->size();
  if(highestEtElectronOnly_) // if true, will only consider highest-pt muon in event
    {
      int iEle = getTheHardestElectron();
      if(iEle >=0) // there is at least one muon in the event
	{
	  iEleMin = iEle;
	  iEleMax = iEle + 1;
	}
    }

  //loop over electrons
  for (int theEle = iEleMin; theEle != iEleMax; ++theEle){//loop over electrons
    bool fill_entry = true; // if true, will histogram muon

    heep::Ele el((*electrons)[theEle]);
    cutCode = -1; // this means HEEP cuts have not been run for this electron
    setElectronMomentum(el);
    for(int cut_index = 0; cut_index != Num_elmet_cuts; ++cut_index)
      { // loop over selection cuts
	string arg = elmet_cuts_desc_short[cut_index];
	// call to function [as implemented in setupCutOder]
	// don't get me started about the syntax!
	bool survived_cut=(this->*cuts[arg])(&fill_entry, el, event);
	if(!survived_cut)break; // skip rest of selection cuts
	
	if(fill_entry)
	  tabulateMe(cut_index, accountMe, event, theEle);
	
	if(dumpHighEtElectrons_ && fill_entry 
	   && cut_index == Num_elmet_cuts-1
	   && wprimeUtil_->getSampleName()=="data"
	   && el.et() > dumpHighEtElectronThreshold_)
	  printHighEtElectron(event);
      
      } // loop over selection cuts


  } // loop over electrons

}

// set electron 4-d momentum (sets el4D)
void EleMETAnalyzer::setElectronMomentum(const heep::Ele & el)
{
  TVector3 p3(el.p4().px(), el.p4().py(), el.p4().pz());
  el4D.SetVectM(p3, wprime::ELECTRON_MASS);
}


// fill histograms for muon if fill_entry=true; update book-keeping 
// (via private member: stats); make sure stats gets updated maximum once per event
void EleMETAnalyzer::tabulateMe(int cut_index, bool accountMe[], 
				edm::EventBase const& event, int theEle)
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
  hPT[cut_index]->Fill(el4D.Pt(), weight);
  hETA[cut_index]->Fill(el4D.Eta(), weight);
  hPHI[cut_index]->Fill(el4D.Phi(), weight);
  hTM[cut_index]->Fill(WPrimeUtil::TMass(el4D, getNewMET(event, el4D)), weight);
}

// operations to be done when changing input file (e.g. create new histograms)
void EleMETAnalyzer::beginFile(std::vector<wprime::InputFile>::const_iterator fi)
{
  wprime::FilterEff tmp; 
  for(int cut = 0; cut != Num_elmet_cuts; ++cut)
    stats[fi->samplename].push_back(tmp);

  // add channel/analysis name here?
  TFileDirectory dir= wprimeUtil_->getFileService()->mkdir(fi->samplename); 
  defineHistos(dir); // one set of histograms per input file

}

// operations to be done when closing input file 
// (e.g. print summary)
void EleMETAnalyzer::endFile(std::vector<wprime::InputFile>::const_iterator fi,
			    ofstream & out)
{
  printFileSummary(fi, out);
}

// print summary of efficiencies
void EleMETAnalyzer::printFileSummary(std::vector<wprime::InputFile>::const_iterator fi, ofstream & out)
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

  for(int cut_index = 0; cut_index != Num_elmet_cuts; ++cut_index){
    
    (it->second)[cut_index].Nsurv_evt_cut_w = 
      (it->second)[cut_index].Nsurv_evt_cut*weight;
    
    out << " Cut # " << cut_index << ": " << elmet_cuts_desc_long[cut_index] 
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
void EleMETAnalyzer::endAnalysis(ofstream & out)
{
  float N_SM = 0; 
  std::map<std::string, std::vector<wprime::FilterEff> >::const_iterator it;

  int index = Num_elmet_cuts-1; // get # of events after last cut

  out << endl;
  for(it = stats.begin(); it != stats.end(); ++it)
    { // loop over samples
      string sample = it->first;
      float N_evt = (it->second)[index].Nsurv_evt_cut_w;
      out<< " "<< sample << ": " << N_evt
	 << " evts (eff = " << 100.*(it->second)[index].eff
	 << " +- " << 100.*(it->second)[index].deff
	 << " %) " << endl;
      
      if(sample == "W" || sample == "Wlowpt" || sample == "QCD" 
	 || sample == "Z" || sample == "Top")
	N_SM += N_evt;
      
    } // loop over samples

  out << " Total # of SM (W + QCD + Z + Top) events: " 
      << N_SM << endl;
  
}


void EleMETAnalyzer::defineHistos_ElectronEt(TFileDirectory & dir)
{
  for(int cut = 0; cut != Num_elmet_cuts; ++cut)
    {
      string name = "hET_" + elmet_cuts_desc_short[cut];
      string title = " Electron E_{T} with " + elmet_cuts_desc_long[cut];
      
      hPT[cut] = dir.make<TH1F>(name.c_str(), title.c_str(), nBinEtEle,
			    minEtEle, maxEtEle);
    }
}

void EleMETAnalyzer::defineHistos_ElectronEta(TFileDirectory & dir)
{
  for(int cut = 0; cut != Num_elmet_cuts; ++cut)
    {
      string name = "hETA_" + elmet_cuts_desc_short[cut];
      string title = " Electron #eta with " + elmet_cuts_desc_long[cut];
      hETA[cut] = dir.make<TH1F>(name.c_str(), title.c_str(), nBinEtaEle,
			   minEtaEle,maxEtaEle);
    }

}

void EleMETAnalyzer::defineHistos_ElectronPhi(TFileDirectory & dir)
{
  for(int cut = 0; cut != Num_elmet_cuts; ++cut)
    {
      string name = "hPHI_" + elmet_cuts_desc_short[cut];
      string title = " Electron #phi with " + elmet_cuts_desc_long[cut];
      hPHI[cut] = dir.make<TH1F>(name.c_str(), title.c_str(), nBinPhiEle,
				 minPhiEle,maxPhiEle);
    }

}

void EleMETAnalyzer::defineHistos_TMass(TFileDirectory & dir)
{
  for(int cut = 0; cut != Num_elmet_cuts; ++cut)
    {
      string name = "hTM_" + elmet_cuts_desc_short[cut];
      string title = " Transv. Mass with " + elmet_cuts_desc_long[cut];
      hTM[cut] = dir.make<TH1F>(name.c_str(), title.c_str(), nBinTmEle,
			  minTmEle,maxTmEle);
    }
  
}

void EleMETAnalyzer::setupCutOrder()
{
  cuts.clear();
#if debugmeElMet
  cout << "\n Ele+MET cuts will be applied in this order: " << endl;
#endif

  for(int cut_i = 0; cut_i != Num_elmet_cuts; ++cut_i)
    { // loop over selection cuts
      string arg = elmet_cuts_desc_short[cut_i];
#if debugmeElMet
      cout << " Cut #" << (cut_i+1) << ": " << elmet_cuts_desc_long[cut_i]
	   << " (" << arg << ") " << endl;
#endif
      if(arg == "hlt")cuts[arg] = &EleMETAnalyzer::passedHLT;
      else if(arg == "qual")cuts[arg] = &EleMETAnalyzer::goodQualityElectron;
      else if(arg == "1el")cuts[arg] = &EleMETAnalyzer::onlyOneHighEtElectron;
      else if(arg == "iso")cuts[arg] = &EleMETAnalyzer::isolatedElectron;
      else if(arg == "met")cuts[arg] = &EleMETAnalyzer::kinematicCuts;
      else
	{
	  cout << " Oops! Don't understand how to prepare for cut nicknamed as "
	       << arg << endl;
	  abort();
	}
    } // loop over selection cuts

#if debugmeElMet
  cout << endl;
#endif

}

// dump on screen info about high-Et electron
void EleMETAnalyzer::printHighEtElectron(edm::EventBase const & event)
{
  cout << " Run # = " << event.id().run() << " Event # = " 
       << event.id().event() << " LS = " << event.id().luminosityBlock() 
       << endl;

  cout << " Electron eta = " << el4D.Eta() << "  phi = " << el4D.Phi()
       << " Et = " << el4D.Et() << endl;

  pat::METCollection::const_iterator oldMET = met->begin();
  TVector2 oldMETv(oldMET->px(), oldMET->py());
  TVector2 newMET = getNewMET(event, el4D);
  cout << " default pfMET = " << oldMET->pt();
  cout << " hadronic-recoil-adjusted pfMET = " << newMET.Mod() << endl;

  cout << " default TM = " << WPrimeUtil::TMass(el4D, oldMETv);
  cout << " hadronic-recoil-adjusted TM = " << WPrimeUtil::TMass(el4D, newMET) << endl;

#if 0
  cout << " P = " << 

  printMe("  pt = ", theEle);
  printMe(" dpt = ", theEle);
  cout << " # of layers: (strip, pixel) " << endl;
  printMe("layers", theEle);
  cout << " # of layers w/o measurement: (strip, pixel) " << endl;
  printMe("layersNoMeas", theEle);
  cout << " # of hits (strip, pixel, muon) " << endl;
  printMe("hits", theEle);
  cout << " Chi2/Ndof " << endl;
  printMe("chi2", theEle);
  cout << " Tracker eta =  " << theEle->tracker.p.Eta()
       << ", Tracker phi = " << theEle->tracker.p.Phi() << endl;
  cout << " # of standalone muon hits = " << theEle->Nmu_hits << endl;
  cout << " Global: " << theEle->AllGlobalElectrons
       << " Tracker: " << theEle->AllTrackerElectrons
       << " Standalone: " << theEle->AllStandAloneElectrons
       << " Global prompt tight: " << theEle->GlobalElectronPromptTight << endl;
#endif
}

// Get new MET: here (unlike the muon case) there is only one correction to be made:
// the hadronic MET component (that needs to be corrected 
// if applyCorrection=true) from Z data; this will be done according to hadronic 
// activity from Z->mumu (not ee?) reconstructed events
TVector2 EleMETAnalyzer::getNewMET(edm::EventBase const & event, const TLorentzVector & el_p)
{
  if(wprimeUtil_->shouldApplyMETCorrection())
    return wprimeUtil_->getHadronicMET(event) - TVector2(el_p.Px(), el_p.Py());

  pat::METCollection::const_iterator pfMET = met->begin();
  return TVector2 (pfMET->px(), pfMET->py());
}

// whether HLT accepted the event
bool EleMETAnalyzer::passedHLT(bool *, const heep::Ele &, edm::EventBase const &)
{
  // needs implementation
  return true;
}

int EleMETAnalyzer::ignoreIsolationMask = (~heep::CutCodes::ISOLEMHADDEPTH1)
  & (~heep::CutCodes::ISOLHADDEPTH2) & (~heep::CutCodes::ISOLPTTRKS);

int EleMETAnalyzer::useOnlyIsolationMask = heep::CutCodes::ISOLEMHADDEPTH1
	  | heep::CutCodes::ISOLHADDEPTH2 | heep::CutCodes::ISOLPTTRKS;

// make sure HEEP cuts are calculated once per electron
void EleMETAnalyzer::runHEEPcuts(const heep::Ele & el)
{
  if(cutCode < 0)
    cutCode = cuts_.getCutCode(el);
}

// check if electron satisfies quality requirements
// fill goodQual; always returns true
bool EleMETAnalyzer::goodQualityElectron(bool * goodQual, const heep::Ele & el, edm::EventBase const &)
{
  runHEEPcuts(el);
  if(cutCode & ignoreIsolationMask)
    *goodQual = false;
  return true;
}

// true if only one electron with Et > the threshold
bool EleMETAnalyzer::onlyOneHighEtElectron(bool *, const heep::Ele &, edm::EventBase const &)
{
  return (nEleAboveThresh(oneEleEtCut_) == 1);
}

// returns # of (global) electrons with with Et above <Et_thresh>
unsigned EleMETAnalyzer::nEleAboveThresh(float Et_thresh)
{
  unsigned N = 0;
  int iEleMin = 0; int iEleMax = electrons->size();
  //loop over electrons
  for (int theEle = iEleMin; theEle != iEleMax; ++theEle){//loop over electrons
    heep::Ele el((*electrons)[theEle]);
    if(el.et() > Et_thresh)
      ++N;
  }

  return N;
}

// set bool flag to true if electron isolated
// always returns true
bool EleMETAnalyzer::isolatedElectron(bool * goodQual, const heep::Ele & el, 
				 edm::EventBase const &)
{
  runHEEPcuts(el);
  if(cutCode & useOnlyIsolationMask)
    *goodQual = false;
  return true;
}

// check if electron, MET pass kinematic cuts, updated goodQual
// always returns true
bool EleMETAnalyzer::kinematicCuts(bool * goodQual, const heep::Ele &, 
				  edm::EventBase const & event)
{
  TVector2 MET = getNewMET(event, el4D);
  float ratio = el4D.Et()/MET.Mod();
  
  TVector2 electron_T(el4D.Px(), el4D.Py());
  float delta_phi = MET.DeltaPhi(electron_T);
  
  if(ratio < 0.4 || ratio > 1.5 || TMath::Abs(delta_phi) < 2.5)
    *goodQual = false;

  return true;
}

