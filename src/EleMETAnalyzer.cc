#include "UserCode/CMGWPrimeGroup/interface/EleMETAnalyzer.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"

#include <TH1F.h>
#include <TLorentzVector.h>

using std::string; using std::cout; using std::endl;

EleMETAnalyzer::EleMETAnalyzer(const edm::ParameterSet& cfg,WPrimeUtil * wprimeUtil) :
  AnalyzerBase(cfg, wprimeUtil), cuts_(cfg)
{
  electronPtThreshold_   = cfg.getParameter<double>("electronPtThreshold");
  oneEleEtCut_   = cfg.getParameter<double>("oneEleEtCut");
  highestEtElectronOnly_ = cfg.getParameter<bool>("highestEtElectronOnly");
  dumpHighEtElectrons_   = cfg.getParameter<bool>("dumpHighEtElectrons");
  dumpHighEtElectronThreshold_ = cfg.getParameter<double>("dumpHighEtElectronThreshold");
  dumpHighMtElectronThreshold_ = cfg.getParameter<double>("dumpHighMtElectronThreshold");
  
  triggerResults_ = cfg.getParameter<edm::InputTag>("triggerResults");
  HLTPathsByName_= cfg.getParameter<std::vector<std::string > >("hltPaths");
  HLTPathsByIndex_.resize(HLTPathsByName_.size());
  
  mkTuple_ = cfg.getParameter<bool>("mkTuple");

  doEoP_ = cfg.getParameter<bool>("doEoP");
  electronEoverPthreshold_  = cfg.getParameter<double>("electronEoverPthreshold");

  setupCutOrder();
  
  analysis = "eleMET";
}

EleMETAnalyzer::~EleMETAnalyzer()
{
}

void EleMETAnalyzer::defineHistos(const TFileDirectory & dir)
{
  AnalyzerBase::defineHistos(dir);
  for(int i = 0; i != Num_elmet_cuts; ++i)
    {
      hPT[i] = hETA[i] = hPHI[i] = /*hMJDPHI[i] =*/ hTM[i] = 0;
      cloneTrees[i] = 0;
    }
   
  if(mkTuple_) { NTuple = 0; defineTrees(dir); }
  defineHistos_ElectronEt(dir);
  defineHistos_ElectronEta(dir);
  defineHistos_ElectronPhi(dir);
  defineHistos_TMass(dir);
}

// get the hardest muon (based on HEEP Et) in event
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
  if(mkTuple_)  ClearWprimeVariables(vars, analysis);

  event.getByLabel(electronsLabel_, electrons);
  event.getByLabel(metLabel_, defMet);

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
    setElectronMomentum(el); // this is needed for the ntuple-making

    Wcand = wprimeUtil_->getNewMETandW(event, el, met, metLabel_);

    if(mkTuple_) {FillNtuple(theEle,iEleMin,event);  NTuple->Fill();}

    if(el.p4().pt() < electronPtThreshold_) continue;
    if(doEoP_ && Wcand.elec()->epIn() > electronEoverPthreshold_ ) continue;

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
	   && wprimeUtil_->getSampleName().find("data") != string::npos
	   && ((el.et() > dumpHighEtElectronThreshold_) ||
	       Wcand.mt() > dumpHighMtElectronThreshold_))
	  printHighEtElectron(event);
      
      } // loop over selection cuts


  } // loop over electrons

}
//Fill Ntuple
void EleMETAnalyzer::FillNtuple(int & theEle, int & iEleMin, edm::EventBase const & event){

  const reco::Candidate *WpMet =  Wcand.met();

  if( theEle==iEleMin ){
    vars.runId   = event.id().run();
    vars.lumiId  = event.id().luminosityBlock();
    vars.eventId = event.id().event();
    TLorentzVector v(WpMet->pt(),WpMet->eta(),WpMet->phi(),WpMet->energy() );
    vars.met = v;
    vars.p_met = &vars.met;
  }  

  vars.ele = el4D;
  vars.p_ele = &vars.ele;
  const heep::Ele *ele = Wcand.elec();

  vars.ele_charge        = ele->charge();
  vars.ele_e1x5          = ele->e1x5();
  vars.ele_e2x5          = ele->e2x5Max();
  vars.ele_e5x5          = ele->e5x5();
  vars.ele_tkIso         = ele->isolPtTrks();
  vars.ele_emIso         = ele->isolEm();
  vars.ele_hadIso        = ele->isolHad();
  vars.ele_hadIso_d1     = ele->isolHadDepth1();
  vars.ele_hadIso_d2     = ele->isolHadDepth2();
  vars.ele_isEB          = ele->isEB();
  vars.ele_isEE          = ele->isEE();
  vars.ele_isEcalDriven  = ele->isEcalDriven();
  vars.ele_sigmaIetaIeta = ele->sigmaIEtaIEta();
  vars.ele_DphiIn        = ele->dPhiIn();
  vars.ele_DetaIn        = ele->dEtaIn();
  vars.ele_HOverE        = ele->hOverE();
  vars.ele_fbrem         = ele->fbrem();
  vars.ele_EOverP        = ele->epIn();
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
  hTM[cut_index]->Fill(Wcand.mt(), weight);
 
  //fill the TTrees
  if(mkTuple_){
    int iEleMin =0;
    FillNtuple(theEle,iEleMin,event);
    cloneTrees[cut_index] -> Fill();
  }
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

  if(mkTuple_) { 
    vars.totEvents = fi->Nprod_evt;
    vars.crossSection = fi->x_sect;
    NTuple->Fill();
  }

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
      
      if(sample.find("data") == string::npos && 
	 sample.find("wprime") == string::npos)
	N_SM += N_evt;
      
    } // loop over samples

  out << " Total # of SM events: " 
      << N_SM << endl;
  
}

void EleMETAnalyzer::defineTrees(const TFileDirectory & dir)
{

  string name = "ntuple";// + elmet_cuts_desc_short[cut];
  NTuple = dir.make<TTree>(name.c_str(), name.c_str());
  InitializeTree(vars, NTuple, analysis);

  for(int cut = 0; cut != Num_elmet_cuts; ++cut)
    {
      string name = "ntu_" + elmet_cuts_desc_short[cut];
      cloneTrees[cut] = dir.make<TTree>(name.c_str(), name.c_str());
      InitializeTree(vars, cloneTrees[cut], analysis);
    }
}


void EleMETAnalyzer::defineHistos_ElectronEt(const TFileDirectory & dir)
{
  for(int cut = 0; cut != Num_elmet_cuts; ++cut)
    {
      string name = "hET_" + elmet_cuts_desc_short[cut];
      string title = " Electron E_{T} with " + elmet_cuts_desc_long[cut];
      
      hPT[cut] = dir.make<TH1F>(name.c_str(), title.c_str(), nBinEtEle,
			    minEtEle, maxEtEle);
    }
}


void EleMETAnalyzer::defineHistos_ElectronEta(const TFileDirectory & dir)
{
  for(int cut = 0; cut != Num_elmet_cuts; ++cut)
    {
      string name = "hETA_" + elmet_cuts_desc_short[cut];
      string title = " Electron #eta with " + elmet_cuts_desc_long[cut];
      hETA[cut] = dir.make<TH1F>(name.c_str(), title.c_str(), nBinEtaEle,
			   minEtaEle,maxEtaEle);
    }

}

void EleMETAnalyzer::defineHistos_ElectronPhi(const TFileDirectory & dir)
{
  for(int cut = 0; cut != Num_elmet_cuts; ++cut)
    {
      string name = "hPHI_" + elmet_cuts_desc_short[cut];
      string title = " Electron #phi with " + elmet_cuts_desc_long[cut];
      hPHI[cut] = dir.make<TH1F>(name.c_str(), title.c_str(), nBinPhiEle,
				 minPhiEle,maxPhiEle);
    }

}

void EleMETAnalyzer::defineHistos_TMass(const TFileDirectory & dir)
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
  cout << "\n Run # = " << event.id().run() << " Event # = " 
       << event.id().event() << " LS = " << event.id().luminosityBlock() 
       << endl;

  cout << " Electron eta = " << el4D.Eta() << "  phi = " << el4D.Phi()
       << " Et = " << el4D.Et() << " GeV " << endl;

  pat::METCollection::const_iterator oldMET = defMet->begin();
  TVector2 oldMETv(oldMET->px(), oldMET->py());
  cout << " default pfMET = " << oldMET->pt() << " GeV" << " corrected MET = " 
       << met.et() << " GeV " << endl << " TM = " << Wcand.mt() << " GeV " << endl;

}

//##Trigger
inline bool
EleMETAnalyzer::SameTrigger(string & A, string & B){
  return (B.find("*") == string::npos) ? !A.compare(B) : !A.compare(0, A.size()-1, B, 0, B.size()-1);
}


// whether HLT accepted the event
bool EleMETAnalyzer::passedHLT(bool *, const heep::Ele &, edm::EventBase const & event)
{
  // retrieve TriggerResults from the event
  edm::Handle<edm::TriggerResults> triggerResults ;
  event.getByLabel(triggerResults_,triggerResults) ;
  
  bool accept = false ;
  
  if (triggerResults.isValid())
    {
      // get trigger names
      const edm::TriggerNames & triggerNames = event.triggerNames(*triggerResults);
      
      unsigned int n = HLTPathsByName_.size() ;
      for (unsigned int i=0; i!=n; i++)
	{
	  HLTPathsByIndex_[i]=triggerNames.triggerIndex(HLTPathsByName_[i]) ;
	}
      
      // empty input vectors (n==0) means any trigger paths
      if (n==0)
	{
	  n=triggerResults->size() ;
	  HLTPathsByName_.resize(n) ;
	  HLTPathsByIndex_.resize(n) ;
	  for ( unsigned int i=0 ; i!=n ; i++)
	    {
              HLTPathsByName_[i]=triggerNames.triggerName(i) ;
              HLTPathsByIndex_[i]=i ;
	    }
	}
      
      // count number of requested HLT paths which have fired
      unsigned int fired=0 ;
      for ( unsigned int i=0 ; i!=n ; i++ )
	{
	  if (HLTPathsByIndex_[i]<triggerResults->size())
	    {
	      if (triggerResults->accept(HLTPathsByIndex_[i]))
		{
		  fired++ ;
		  accept = true ;
		}
	    }
	}
      
    }
  
  return accept;
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
  float ratio = el4D.Et()/met.et();
  float delta_phi = Wcand.calcDPhi();
  
  if(ratio < 0.4 || ratio > 1.5 || TMath::Abs(delta_phi) < 2.5)
    *goodQual = false;

  return true;
}

void EleMETAnalyzer::fillHistos(const int& index, const float& weight){
}
