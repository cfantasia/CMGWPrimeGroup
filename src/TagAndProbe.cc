#include "UserCode/CMGWPrimeGroup/interface/TagAndProbe.h"

using namespace std;

TagAndProbe::TagAndProbe(){}
TagAndProbe::TagAndProbe(const edm::ParameterSet & cfg, int fileToRun) :
  AnalyzerBase(cfg, fileToRun){

  setupCutOrder();
  if(debug_) printf("Using %i cuts\n",NCuts_);

  removeTauEvents_ = cfg.getUntrackedParameter<bool>("removeTauEvents", false);
  adjustMETPhi_ = cfg.getUntrackedParameter<bool>("adjustMETPhi", false);
  elScaleFactor_ = cfg.getParameter<double>("elScaleFactor");
  muScaleFactor_ = cfg.getParameter<double>("muScaleFactor");
  cout<<"removeTauEvents:"<<removeTauEvents_<<endl
      <<"adjustMETPhi:"<<adjustMETPhi_<<endl
      <<"elScaleFactor:"<<elScaleFactor_<<endl
      <<"muScaleFactor:"<<muScaleFactor_<<endl;

  rhoFastJetLabel_ = cfg.getParameter<edm::InputTag>("rhoFastJet");
  effectiveElecArea_ = cfg.getParameter<vector<double> >("effectiveElecArea");
  effectiveMuonArea_ = cfg.getParameter<vector<double> >("effectiveMuonArea");
  
// +++++++++++++++++++General Cut values
  maxNumZs_ = cfg.getParameter<uint>("maxNumZs");
  minLeadPt_ = cfg.getParameter<double>("minLeadPt");
  minNLeptons_ = cfg.getUntrackedParameter<uint>("minNLeptons", 0);
  maxNVLLeptons_ = cfg.getUntrackedParameter<uint>("maxNVLLeptons", 999);
  minNTightLeptons_ = cfg.getUntrackedParameter<uint>("minNTightLeptons", 0);
  maxNJets_ = cfg.getUntrackedParameter<uint>("maxNJets", 999);
 
  minMET_ = cfg.getUntrackedParameter<double>("minMET", 0.);

// +++++++++++++++++++Z Cuts
  minZeePt1_ = cfg.getParameter<double>("minZeePt1");
  minZeePt2_ = cfg.getParameter<double>("minZeePt2");
  minZmmPt1_ = cfg.getParameter<double>("minZmmPt1");
  minZmmPt2_ = cfg.getParameter<double>("minZmmPt2");
  minZmass_ = cfg.getUntrackedParameter<double>("minZmass", 0.);
  maxZmass_ = cfg.getUntrackedParameter<double>("maxZmass", 9e9);

  maxZMassDiff_ = cfg.getParameter<double>("maxZMassDiff");

  //Selectors
  Pset eSelectorPset = cfg.getParameter<Pset>("electronSelectors");
  string countElectronType = cfg.getUntrackedParameter<string>("CountElectronType", "wp95");
  cout<<"Using "<<countElectronType<<" for count electrons "<<endl;
  countElectron_ = ElectronSelector(eSelectorPset, countElectronType);
  string TagElectronType = cfg.getUntrackedParameter<string>("TagElectronType", "wp95");
  cout<<"Using "<<TagElectronType<<" for Tag electrons "<<endl;
  TagElectron_ = ElectronSelector(eSelectorPset, TagElectronType);
  string looseProbeElectronType = cfg.getUntrackedParameter<string>("LooseProbeElectronType", "wp95");
  cout<<"Using "<<looseProbeElectronType<<" for loose Probe electrons "<<endl;
  looseProbeElectron_ = ElectronSelector(eSelectorPset, looseProbeElectronType);
  string tightProbeElectronType = cfg.getUntrackedParameter<string>("TightProbeElectronType", "wp95");
  cout<<"Using "<<tightProbeElectronType<<" for tight Probe electrons "<<endl;
  tightProbeElectron_ = ElectronSelector(eSelectorPset, tightProbeElectronType);

  Pset mSelectorPset = cfg.getParameter<Pset>("muonSelectors");
  string countMuonType = cfg.getUntrackedParameter<string>("CountMuonType", "exotica");
  cout<<"Using "<<countMuonType<<" for count muons"<<endl;
  countMuon_ = MuonSelector(mSelectorPset, countMuonType);
  string TagMuonType = cfg.getUntrackedParameter<string>("TagMuonType", "exotica");
  cout<<"Using "<<TagMuonType<<" for Tag muons"<<endl;
  TagMuon_ = MuonSelector(mSelectorPset, TagMuonType);
  string looseProbeMuonType = cfg.getUntrackedParameter<string>("LooseProbeMuonType", "exotica");
  cout<<"Using "<<looseProbeMuonType<<" for loose Probe muons"<<endl;
  looseProbeMuon_ = MuonSelector(mSelectorPset, looseProbeMuonType);
  string tightProbeMuonType = cfg.getUntrackedParameter<string>("TightProbeMuonType", "exotica");
  cout<<"Using "<<tightProbeMuonType<<" for tight Probe muons"<<endl;
  tightProbeMuon_ = MuonSelector(mSelectorPset, tightProbeMuonType);

  countElectrons_.reserve(20);
  TagElectrons_.reserve(20);
  looseProbeElectrons_.reserve(20);
  tightProbeElectrons_.reserve(20);
  countMuons_.reserve(20);
  TagMuons_.reserve(20);
  looseProbeMuons_.reserve(20);
  tightProbeMuons_.reserve(20);

  Pset jSelectorPset = cfg.getParameter<Pset>("jetSelectors");
  string looseJetType = cfg.getUntrackedParameter<string>("LooseJetType", "Base");
  cout<<"Using "<<looseJetType<<" for jets\n";
  looseJet_ = JetSelector(jSelectorPset, looseJetType);

  looseJets_.reserve(30);
  

}

TagAndProbe::~TagAndProbe(){
}

void TagAndProbe::setupCutOrder(){
  map<string,fnCut > mFnPtrs;
  mFnPtrs["NoCuts"] = boost::bind(&TagAndProbe::passNoCut, this);
  mFnPtrs["HLT"] = boost::bind(&TagAndProbe::passTriggersCut, this);
  mFnPtrs["MinNLeptons"] = boost::bind(&TagAndProbe::passMinNLeptonsCut, this, boost::cref(TagElectrons_), boost::cref(TagMuons_), boost::cref(minNLeptons_));
  mFnPtrs["MaxNVLLeptons"] = boost::bind(&TagAndProbe::passMaxNLeptonsCut, this, boost::cref(countElectrons_), boost::cref(countMuons_), boost::cref(maxNVLLeptons_));
  mFnPtrs["ValidZ"] = boost::bind(&TagAndProbe::passValidZCut, this, boost::ref(zCand_));
  mFnPtrs["NumZs"] = boost::bind(&TagAndProbe::passNumberOfZsCut, this);
  mFnPtrs["ZMass"] = boost::bind(&TagAndProbe::passZMassCut, this, boost::cref(zCand_), boost::cref(minZmass_), boost::cref(maxZmass_));
  mFnPtrs["MET"] = boost::bind(&TagAndProbe::passMinMETCut, this, boost::cref(met_), boost::cref(minMET_));

  fillCuts(mFnPtrs);

}

void TagAndProbe::defineHistos(const TFileDirectory & dir){
  if(debug_) printf("Declare histos\n");
    tEvts.assign(NCuts_,NULL);
    for(int i=0; i<NCuts_; ++i){
      string name = "tEvts_" + CutNames_[i];
      string title = "Analysis Variables (After " + CutDescs_[i] + " Cut)";
      tEvts[i] = dir.make<TTree>(name.c_str(), title.c_str());
      tEvts[i]->Branch("Run", &runNumber_);
      tEvts[i]->Branch("Lumi", &lumiNumber_);
      tEvts[i]->Branch("Event", &evtNumber_);
      tEvts[i]->Branch("EvtType", &evtType_);
      tEvts[i]->Branch("ZMass", &ZMass_);
      tEvts[i]->Branch("MET", &MET_);
      tEvts[i]->Branch("NVtxs", &NVtxs_);
      tEvts[i]->Branch("weight", &weight_);
      tEvts[i]->Branch("ZLep1Pt", &ZLep1Pt_);
      tEvts[i]->Branch("ZLep1Eta", &ZLep1Eta_);
      tEvts[i]->Branch("ZLep1Phi", &ZLep1Phi_);
      tEvts[i]->Branch("ZLep2Pt", &ZLep2Pt_);
      tEvts[i]->Branch("ZLep2Eta", &ZLep2Eta_);
      tEvts[i]->Branch("ZLep2Phi", &ZLep2Phi_);
      tEvts[i]->Branch("ZLep1PtGen", &ZLep1PtGen_);
      tEvts[i]->Branch("ZLep1EtaGen", &ZLep1EtaGen_);
      tEvts[i]->Branch("ZLep1PhiGen", &ZLep1PhiGen_);
      tEvts[i]->Branch("ZLep2PtGen", &ZLep2PtGen_);
      tEvts[i]->Branch("ZLep2EtaGen", &ZLep2EtaGen_);
      tEvts[i]->Branch("ZLep2PhiGen", &ZLep2PhiGen_);
      tEvts[i]->Branch("ZTightCode", &ZTightCode_);

    }
}

//fill Histograms
void TagAndProbe::fillHistos(const int& index, const float& weight){
  if(debug_) printf("filling Histos\n");
  if(index > 0) tEvts[index]->Fill();//trying to keep the file size down
}

void 
TagAndProbe::eventLoop(edm::EventBase const & event){
  clearEvtVariables();
  runNumber_ = event.id().run();
  lumiNumber_ = event.id().luminosityBlock();
  evtNumber_ = event.id().event();

/*
  if(!wprimeUtil_->runningOnData()){//Don't do this for data
    GenParticleV genParticles = getProduct<GenParticleV>(event, "genParticles");
    for (size_t i = 0; i < genParticles.size(); i++){
      if (abs(genParticles[i].pdgId()) == PDGTAU){
        return;
      }
    }
  }
*/  
  if(debug_){
    WPrimeUtil::printEvent(event, cout);
  }

  if(removeTauEvents_ && !wprimeUtil_->runningOnData()){//Don't do this for data
    GenParticleV genParticles = getProduct<GenParticleV>(event, "prunedGenParticles");
      for (size_t i = 0; i < genParticles.size(); i++){
        if (abs(genParticles[i].pdgId()) == PDG_ID_TAU ||
            abs(genParticles[i].pdgId()) == PDG_ID_TAUNEU){
          if(debug_) cout<<"Removing Tau event."<<endl;
          return;
        }
      }
  }//Remove Tau Events
 
  // get leptons
  event.getByLabel(electronsLabel_,patElectronsH_);
  event.getByLabel(muonsLabel_,patMuonsH_);
  event.getByLabel(metLabel_, metH_);

  if(useAdjustedMET_) event.getByLabel(pfCandsLabel_, pfCandidatesH_);
  WPrimeUtil::getLeptonsMET(patElectronsH_, allElectrons_,
                            patMuonsH_, muReconstructor_, allMuons_,
                            metH_, useAdjustedMET_, met_,
                            pfCandidatesH_);
  
  //get Vertex
  event.getByLabel(vertexLabel_, verticesH_);

  if(debug_) cout<<"num of vertices is "<<verticesH_->size()<<endl;
  const reco::Vertex & primaryVertex =  verticesH_.isValid() && !verticesH_->empty() ? verticesH_->at(0) : reco::Vertex();
  if(!verticesH_.isValid() || verticesH_->empty()) cout<<" No valid PV"<<endl;

  if(adjustMETPhi_){//Correct for met phi shift
    TVector2 subtract; 
    int Nvtx = 0;
    for(unsigned ivtx=0; ivtx<verticesH_->size(); ++ivtx){
      if(verticesH_.isValid() && 
         verticesH_->at(ivtx).ndof() >= 4 && 
         verticesH_->at(ivtx).chi2() > 0 && 
         verticesH_->at(ivtx).tracksSize() > 0 && 
         fabs(verticesH_->at(ivtx).z()) < 24 && 
         fabs(verticesH_->at(ivtx).position().Rho()) < 2.)  Nvtx++;

      if(debug_) cout<<"checking vertix "<<ivtx<<endl;
    }
    if(debug_) cout<<"Found "<<Nvtx<<" vertices for met phi shift"<<endl;
    if(wprimeUtil_->runningOnData()){
      subtract.Set(+3.54233e-01 + 2.65299e-01*Nvtx,
                   +1.88923e-01 - 1.66425e-01*Nvtx);
      //pfMEtSysShiftCorrParameters_2012runAvsNvtx_data = cms.PSet(
      //  px = cms.string("+3.54233e-01 + 2.65299e-01*Nvtx"),
      //  py = cms.string("+1.88923e-01 - 1.66425e-01*Nvtx")
    }else{//MC
      subtract.Set(-2.99576e-02 - 6.61932e-02*Nvtx,
                   +3.70819e-01 - 1.48617e-01*Nvtx);
      //pfMEtSysShiftCorrParameters_2012runAvsNvtx_mc = cms.PSet(
      //  px = cms.string("-2.99576e-02 - 6.61932e-02*Nvtx"),
      //  py = cms.string("+3.70819e-01 - 1.48617e-01*Nvtx")
    }
    TVector2 newmet(met_.px()-subtract.Px(), met_.py()-subtract.Py());
    if(debug_) cout<<"Before shifting phi:"<<met_.p4()<<endl;
    met_.setP4(LorentzVector(newmet.Px(), newmet.Py(), 0., newmet.Mod()));
    if(debug_) cout<<"After shifting phi:"<<met_.p4()<<endl;
  }

  event.getByLabel(rhoFastJetLabel_, rhoFastJetH_);
  if(debug_){
    cout<<" rho = "<<(*rhoFastJetH_)<<endl;
    printf("    Contains: %lu electron(s), %lu muon(s)\n",
           allElectrons_.size(), allMuons_.size());
    print(allElectrons_, "All Electrons");
    print(allMuons_,     "All Muons");
  }

  // Make vectors of leptons passing various criteria

  PatElectronV scaledPatElectrons;  
  if(elScaleFactor_) scaledPatElectrons = *patElectronsH_;
    
  for (size_t i = 0; i < allElectrons_.size(); i++) {
    if(elScaleFactor_){
      pat::Electron & newEl = scaledPatElectrons[i];
      newEl.setP4(elScaleFactor_*newEl.p4()); 
      allElectrons_[i] = heep::Ele(newEl);
    }
    const pat::Electron & p = allElectrons_[i].patEle();
    if(Overlap(p, *patMuonsH_.product(), 0.01)) continue;
    const float pu = ElecPU(allElectrons_[i]);

    if (countElectron_(p, pu, primaryVertex))
      countElectrons_.push_back(allElectrons_[i]);

    if (TagElectron_(p, pu, primaryVertex))
      TagElectrons_.push_back(allElectrons_[i]);

    if (looseProbeElectron_(p, pu, primaryVertex))
      looseProbeElectrons_.push_back(allElectrons_[i]);

    if (tightProbeElectron_(p, pu, primaryVertex))
      tightProbeElectrons_.push_back(allElectrons_[i]);

  }

  for (size_t i = 0; i < allMuons_.size(); i++) {
    TeVMuon & p = allMuons_[i];
    if(muScaleFactor_) p.setP4(muScaleFactor_*p.p4());

    const float pu = MuonPU(p);
    if(debug_) cout<<"muon number "<<i<<" has pu "<<pu<<endl;

    if (countMuon_(p, pu, primaryVertex))
      countMuons_.push_back(p);

    if (TagMuon_(p, pu, primaryVertex))
      TagMuons_.push_back(p);

    if (looseProbeMuon_(p, pu, primaryVertex))
      looseProbeMuons_.push_back(p);

    if (tightProbeMuon_(p, pu, primaryVertex))
      tightProbeMuons_.push_back(p);
  }

  for(size_t i = 0; i < TagMuons_.size(); i++) {
    const TeVMuon & m = TagMuons_[i];
    if(!WPrimeUtil::Contains(m, tightProbeMuons_)){
      cout<<"This tag muon isn't a probe muon!"<<endl;
      print(TagMuons_,     "Tag Muons");
      print(tightProbeMuons_,     "Tight Probe Muons");
      abort();
    }
  }

  if(debug_){
    printf("    Contains: %lu count electron(s), %lu count muon(s)\n",
           countElectrons_.size(), countMuons_.size());
    printf("    Contains: %lu tag electron(s), %lu tag muon(s)\n",
           TagElectrons_.size(), TagMuons_.size());
    printf("    Contains: %lu loose Probe electron(s), %lu loose Probe muon(s)\n",
           looseProbeElectrons_.size(), looseProbeMuons_.size());
    printf("    Contains: %lu tight Probe electron(s), %lu tight Probe muon(s)\n",
           tightProbeElectrons_.size(), tightProbeMuons_.size());
  }

  if(debug_){
    cout<<"Tag Leptons: \n";
    print(TagElectrons_, "Tag Electrons");
    print(TagMuons_,     "Tag Muons");
    cout<<"Loose Probe Leptons: \n";
    print(looseProbeElectrons_, "Loose Probe Electrons");
    print(looseProbeMuons_,     "Loose Probe Muons");
    cout<<"Tight Probe Leptons: \n";
    print(tightProbeElectrons_, "Tight Probe Electrons");
    print(tightProbeMuons_,     "Tight Probe Muons");
  }

  // Preselection - skip events that don't look promising
  if (doPreselect_){
    if(debug_) cout<<"Testing Preselection...\n";
    if(TagElectrons_.size() + TagMuons_.size() == 0) return;
  }

  //get Jets
  if(debug_) cout<<"Fill Jet Vectors\n";
  event.getByLabel(jetsLabel_, patJetsH_);
  const JetV & allJets  = *patJetsH_;
  for (size_t i = 0; i < allJets.size(); i++) {
    if (looseJet_(allJets[i]) && !Overlap(allJets[i], TagMuons_, 1.0) && !Overlap(allJets[i], TagElectrons_, 1.0))
      looseJets_.push_back(allJets[i]);
  }

  if(debug_) cout<<"Get Trigger\n";
  //get Trigger 
  event.getByLabel(hltEventLabel_, triggerEventH_);
  //if(debug_) printDebugEvent();

  MET_ = met_.et();
  NVtxs_ = (*verticesH_).size();


  if(!wprimeUtil_->runningOnData()){//Don't do this for data
    GenParticleV genParticles = getProduct<GenParticleV>(event, genLabel_);
    for (size_t i = 0; i < genParticles.size(); i++){
      const reco::Candidate * gen = &genParticles[i];
      if(gen->status() != 3) continue; //I think this is before fsr
      int pdgid = abs(gen->pdgId());
      if(pdgid == PDG_ID_ELEC || pdgid == PDG_ID_MUON ||
         pdgid == PDG_ID_ELECNEU || pdgid == PDG_ID_MUONNEU){
        if(debug_) cout<<"Found particle with pdgid, status, pt = "<<gen->pdgId()<<", "<<gen->status()<<", "<<gen->pt()<<endl;
        const reco::Candidate * genMother = findMother(gen);
        if(!genMother) continue;
        if(debug_) cout<<"Found mother with pdgid, status, pt = "<<genMother->pdgId()<<", "<<genMother->status()<<", "<<genMother->pt()<<endl;
        if      (abs(genMother->pdgId()) == PDG_ID_Z && gen->charge() > 0){
          ZLep1PtGen_  = gen->pt();
          ZLep1EtaGen_ = gen->eta();
          ZLep1PhiGen_ = gen->phi();
        }else if(abs(genMother->pdgId()) == PDG_ID_Z && gen->charge() < 0){
          ZLep2PtGen_  = gen->pt();
          ZLep2EtaGen_ = gen->eta();
          ZLep2PhiGen_ = gen->phi();
        }
      }
    }
  }//MC Only If


  weight_ = wprimeUtil_->getWeight();
  if(!passCuts(weight_)) return;
  if(wprimeUtil_->runningOnData()){
    //cout<<" The following data event passed All Cuts!!!\n";
    //printPassingEvent(event);
    //if(debug_) printEventLeptons();
    //cout<<" ------------------\n";
  }
  ///if(debug_) printEventLeptons();

}

////////////////////////////////
/////////Check Z Properties/////
////////////////////////////////

inline bool
TagAndProbe::passValidZCut(ZCandidate& z){
  calcZVariables();//Cory: These should take args
  return AnalyzerBase::passValidZCut(z);
}

inline bool TagAndProbe::passZLepTightCut(bool firstDaughter){
  if(zCand_.flavor() == PDG_ID_ELEC){
    const heep::Ele & e = firstDaughter ? *zCand_.elec1() : *zCand_.elec2();
    return WPrimeUtil::Contains(e, tightProbeElectrons_);
  }else if(zCand_.flavor() == PDG_ID_MUON){
    const TeVMuon & m = firstDaughter ? *zCand_.muon1() : *zCand_.muon2();
    return WPrimeUtil::Contains(m, tightProbeMuons_);
  }
  return false;
}

inline bool
TagAndProbe::passNumberOfZsCut() const{
  return numZs_ <= maxNumZs_;
}

int
TagAndProbe::countZCands(ZCandV & Zs) const{
  int count =0;
  for(uint i=0; i<Zs.size(); ++i) 
    if( (Zs[i].mass() > minZmass_) && (Zs[i].mass() < maxZmass_)) 
      count++;
  return count;
}

void
TagAndProbe::calcZVariables(){
  if (debug_) cout<<"In calc Z Variables"<<endl;
  // Reconstruct the Z
  ZCandV zeeCands = ZCandsAsym(TagElectrons_, tightProbeElectrons_);
  if(zeeCands.size() == 0){
    if(debug_) cout<<" Making Z out of loose probe electrons"<<endl;
    zeeCands = ZCandsAsym(TagElectrons_, looseProbeElectrons_);
  }
  removeLowLepPtCands(zeeCands, minZeePt1_, minZeePt2_);

  ZCandV zmmCands = ZCandsAsym(TagMuons_, tightProbeMuons_);
  if(zmmCands.size() == 0){
    if(debug_) cout<<" Making Z out of loose probe muons"<<endl;
    zmmCands = ZCandsAsym(TagMuons_,     looseProbeMuons_);
  }
  removeLowLepPtCands(zmmCands, minZmmPt1_, minZmmPt2_);

  if(debug_){
    printf("    Contains: %lu Zee and %lu Zmm preliminary candidate(s)\n", zeeCands.size(), zmmCands.size());
    for (ZCandV::iterator i = zeeCands.begin(); i != zeeCands.end(); ++i){
      cout<<"Electron Z Cand with mass: "<<i->mass()<<endl;
    }
    for (ZCandV::iterator i = zmmCands.begin(); i != zmmCands.end(); ++i){
      cout<<"Muon Z Cand with mass: "<<i->mass()<<endl;
    }
  }
  
  ZCandV zCands;
  zCands.insert(zCands.end(), zeeCands.begin(), zeeCands.end());
  zCands.insert(zCands.end(), zmmCands.begin(), zmmCands.end());

  removeWorstCands(zCands, minZmass_, maxZmass_);
  sort(zCands.begin(), zCands.end(), closestToZMass());
  removeOverlapping(zCands);

  if(debug_) printf("    Contains: %lu Z final candidate(s)\n", zCands.size());
  zCand_ = zCands.size() ? zCands[0] : ZCandidate();

  ZCandV zCandsAll;
  if(zCand_){
    if(debug_) printf("Good Cand: Flavor: %i Mass: %.2f\n", zCand_.flavor(), zCand_.mass());
    
    ZCandV zeeCandsAll = getZCands(looseProbeElectrons_);
    ZCandV zmmCandsAll = getZCands(looseProbeMuons_    );
    if(debug_) printf("    Contains: %lu loose preliminary Zee candidate(s) and %lu loose preliminary Zmm candidate(s)\n", zeeCandsAll.size(), zmmCandsAll.size());

    zCandsAll.insert(zCandsAll.end(), zeeCandsAll.begin(), zeeCandsAll.end());
    zCandsAll.insert(zCandsAll.end(), zmmCandsAll.begin(), zmmCandsAll.end());
    removeWorstCands(zCandsAll, minZmass_, maxZmass_);

    bool tight1=false, tight2=false;
    tight1 = passZLepTightCut(true);
    tight2 = passZLepTightCut(false);
    ZTightCode_ = (2*tight1) + tight2;
  }
  if(zCand_) zCandsAll.insert(zCandsAll.begin(), zCand_);
  removeOverlapping(zCandsAll);
  numZs_ = countZCands(zCandsAll); 
  ZMass_ = zCand_.mass();
  if(zCand_){
    evtType_ = 2 * (zCand_.flavor() != 11);
    ZLep1Pt_  = zCand_.daughter(0)->pt();
    ZLep1Eta_ = zCand_.daughter(0)->eta();
    ZLep1Phi_ = zCand_.daughter(0)->phi();
    ZLep2Pt_  = zCand_.daughter(1)->pt();
    ZLep2Eta_ = zCand_.daughter(1)->eta();
    ZLep2Phi_ = zCand_.daughter(1)->phi();
  }
  if(debug_){
    printf("    Contains: %lu loose Z candidate(s)\n", zCandsAll.size());
    for (uint i=0; i<zCandsAll.size(); ++i){
      printf(" Cand: %u Flavor: %i Mass: %.2f\n", i, zCandsAll[i].flavor(), zCandsAll[i].mass());
    }
    //printEventLeptons();
    //printEventDetails();
  }
}
  
void TagAndProbe::printDebugEvent() const{
  WPrimeUtil::printPassingTriggers(*triggerEventH_,triggersToUse_);
  printEventDetails();
  //printEventLeptons();
}

void TagAndProbe::printEventDetails() const{
  if(zCand_){
    cout<<" Z Flavor: "<<zCand_.flavor()
        <<" Z Mass: "<<ZMass_
        <<" Z Tight Code: "<<ZTightCode_
        <<" Z lep1 pt "<<zCand_.daughter(0)->pt()
        <<" Z lep1 eta "<<zCand_.daughter(0)->eta()
        <<" Z lep1 phi "<<zCand_.daughter(0)->phi()
        <<" Z lep2 pt "<<zCand_.daughter(1)->pt()
        <<" Z lep2 eta "<<zCand_.daughter(1)->eta()
        <<" Z lep2 phi "<<zCand_.daughter(1)->phi()
        <<endl;
  }
}

inline void
TagAndProbe::clearEvtVariables(){
  looseJets_.clear();
  allElectrons_.clear(); countElectrons_.clear();
  TagElectrons_.clear(); 
  looseProbeElectrons_.clear(); tightProbeElectrons_.clear();
  allMuons_.clear(); countMuons_.clear();
  TagMuons_.clear(); 
  looseProbeMuons_.clear(); tightProbeMuons_.clear();
  met_ = pat::MET();
  zCand_ = ZCandidate();
  LeadPt_ = -999;
  evtType_ = -999;
  LeadElecPt_ = -999;
  LeadMuonPt_ = -999;
  ZMass_ = -999;
  ZTightCode_ = 0;
  runNumber_ = 0;
  lumiNumber_ = 0;
  evtNumber_ = 0;
  MET_ = 0;
  NVtxs_ = 0;
  ZLep1Pt_ = ZLep1Eta_ = ZLep1Phi_ = ZLep2Pt_ = ZLep2Eta_ = ZLep2Phi_ = 0;
  ZLep1PtGen_ = ZLep1EtaGen_ = ZLep1PhiGen_ = ZLep2PtGen_ = ZLep2EtaGen_ = ZLep2PhiGen_ = 0.;
  weight_ = 0;
}

inline float
  TagAndProbe::ElecPU(const heep::Ele & e) const{
  float effArea = 0.;
  float abs_eta = fabs(e.patEle().eta());
  if     (                 abs_eta<1.0  )   effArea = 0.13 ; //+/- 0.001
  else if(abs_eta>1.0   && abs_eta<1.479)   effArea = 0.14 ; //+/- 0.002
  else if(abs_eta>1.479 && abs_eta<2.0  )   effArea = 0.07 ; //+/- 0.001
  else if(abs_eta>2.0   && abs_eta<2.2  )   effArea = 0.09 ; //+/- 0.001
  else if(abs_eta>2.2   && abs_eta<2.3  )   effArea = 0.11 ; //+/- 0.002
  else if(abs_eta>2.3   && abs_eta<2.4  )   effArea = 0.11 ; //+/- 0.003
  else if(abs_eta>2.4                   )   effArea = 0.14 ; //+/- 0.004 
  else{ cout<<"Electron with eta out of range: "<<abs_eta<<endl; abort();}
  return (*rhoFastJetH_)*effArea;

  //return (*rhoFastJetH_)*effectiveElecArea_[e.patEle().isEE()];
}

inline float
TagAndProbe::MuonPU(const TeVMuon & m) const{
    float pu = 0;
    for (size_t j = 0; j < allMuons_.size(); j++) {
      if(areIdentical(m, allMuons_[j])) continue;//same muon
      if(deltaR(m, allMuons_[j]) < 0.4) pu += allMuons_[j].pt(); 
    }
    return pu;
}

