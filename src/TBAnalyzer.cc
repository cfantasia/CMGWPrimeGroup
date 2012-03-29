#include "UserCode/CMGWPrimeGroup/interface/TBAnalyzer.h"

using namespace std;

TBAnalyzer::TBAnalyzer() : MAX_MBL(-1){}
TBAnalyzer::TBAnalyzer(const edm::ParameterSet & cfg, int fileToRun) :
  AnalyzerBase(cfg, fileToRun),
  MAX_MBL(sqrt(TMASS*TMASS - WMASS*WMASS))
{
  setupCutOrder();
  if(debug_) printf("Using %i cuts\n",NCuts_);

  
// +++++++++++++++++++General Cut values
  minNLeptons_ = cfg.getUntrackedParameter<uint>("minNLeptons", 0);
  minNJets_ = cfg.getUntrackedParameter<uint>("minNJets", 0);  
  minNBJets_ = cfg.getUntrackedParameter<uint>("minNBJets", 0);  
  minMET_ = cfg.getUntrackedParameter<double>("minMET", 0.);

  minWtransMass_ = cfg.getUntrackedParameter<double>("minWtransMass", 0.);
  minWpt_ = cfg.getUntrackedParameter<double>("minWpt", 0.);
  nuAlgo_ = (NuAlgos)cfg.getUntrackedParameter<int>("NuAlgo", kMinPz);

  minTpt_ = cfg.getUntrackedParameter<double>("minTpt", 0.);
  minTMass_ = cfg.getUntrackedParameter<double>("minTMass", 0.);
  maxTMass_ = cfg.getUntrackedParameter<double>("maxTMass", 9e9);

  minBDisc_ = cfg.getUntrackedParameter<double>("minBDisc", 0.);
  BDisc_ = cfg.getUntrackedParameter<string>("BDisc", "trackCountingHighEffBJetTags");
  maxBMass_ = cfg.getUntrackedParameter<double>("maxBMass", 9e9);
  minBpt_ = cfg.getUntrackedParameter<double>("minBpt", 0);
  
  //Selectors
  Pset eSelectorPset = cfg.getParameter<Pset>("electronSelectors");
  string looseElectronType = cfg.getUntrackedParameter<string>("LooseElectronType", "wp95");
  looseElectron_ = ElectronSelector(eSelectorPset, looseElectronType);
  if(debug_) cout<<"Using "<<looseElectronType<<" for loose electrons\n";

  Pset mSelectorPset = cfg.getParameter<Pset>("muonSelectors");
  string looseMuonType = cfg.getUntrackedParameter<string>("LooseMuonType", "exotica");
  looseMuon_ = MuonSelector(mSelectorPset, looseMuonType);
  if(debug_) cout<<"Using "<<looseMuonType<<" for loose muons\n";

  Pset jSelectorPset = cfg.getParameter<Pset>("jetSelectors");
  string looseJetType = cfg.getUntrackedParameter<string>("LooseJetType", "Base");
  looseJet_ = JetSelector(jSelectorPset, looseJetType);
  if(debug_) cout<<"Using "<<looseJetType<<" for jets\n";

}

TBAnalyzer::~TBAnalyzer(){
}

void TBAnalyzer::setupCutOrder(){
  map<string,fnCut > mFnPtrs;

  mFnPtrs["NoCuts"] = boost::bind(&TBAnalyzer::passNoCut, this);
  mFnPtrs["HLT"] = boost::bind(&TBAnalyzer::passTriggersCut, this);
  mFnPtrs["MinNLeptons"] = boost::bind(&TBAnalyzer::passMinNLeptonsCut, this, boost::cref(looseElectrons_), boost::cref(looseMuons_), boost::cref(minNLeptons_));
  mFnPtrs["MinNJets"] = boost::bind(&TBAnalyzer::passMinNJetsCut, this, boost::cref(looseJets_), boost::cref(minNJets_));
  mFnPtrs["MinNBJets"] = boost::bind(&TBAnalyzer::passMinNJetsCut, this, boost::cref(looseBJets_), boost::cref(minNBJets_));
  mFnPtrs["ValidW"] = boost::bind(&TBAnalyzer::passValidWCut, this, boost::cref(wCand_));
  mFnPtrs["ValidB"] = boost::bind(&TBAnalyzer::passValidBCut, this, boost::ref(bCand1_));
  mFnPtrs["ValidT"] = boost::bind(&TBAnalyzer::passValidTCut, this, boost::ref(tCand_));
  mFnPtrs["TMass"] = boost::bind(&TBAnalyzer::passXWMassCut, this, boost::ref(tCand_), nuAlgo_, minTMass_, maxTMass_);
  mFnPtrs["TPt"] = boost::bind(&TBAnalyzer::passXWptCut, this, boost::ref(tCand_), minTpt_);
  mFnPtrs["ValidTB"] = boost::bind(&TBAnalyzer::passValidTBCut, this, boost::ref(tbCand_));
  mFnPtrs["WTransMass"] = boost::bind(&TBAnalyzer::passWtransMassCut, this, boost::cref(wCand_), boost::cref(minWtransMass_));
  mFnPtrs["MET"] = boost::bind(&TBAnalyzer::passMinMETCut, this, boost::cref(met_), boost::cref(minMET_));
  mFnPtrs["Wpt"] = boost::bind(&TBAnalyzer::passWptCut, this, boost::cref(wCand_), boost::cref(minWpt_));
  mFnPtrs["BPt1"] = boost::bind(&TBAnalyzer::passMinPtCut, this, boost::ref(bCand1_), minBpt_);
  mFnPtrs["BPt2"] = boost::bind(&TBAnalyzer::passMinPtCut, this, boost::ref(bCand2_), minBpt_);
  mFnPtrs["AllCuts"] = boost::bind(&TBAnalyzer::passNoCut, this);
  
  fillCuts(mFnPtrs); 
}

void TBAnalyzer::defineHistos(const TFileDirectory & dir){
//TB Histos
  defineHistoSet("hTBMass", "Reconstructed TB Invariant Mass",
                 "M_{tb} (GeV)", 250, 0, 2500, "GeV", hTBMass,dir);
  defineHistoSet("hTBenuMass", "Reconstructed TB(e#nu) Invariant Mass",
                 "M_{tb}^{e#nu} (GeV)", 250, 0, 2500, "GeV", hTBenuMass,dir);
  defineHistoSet("hTBmnuMass", "Reconstructed TB(#mu#nu) Invariant Mass",
                 "M_{tb}^{#mu#nu} (GeV)", 250, 0, 2500, "GeV", hTBmnuMass,dir);
    defineHistoSet("hQ", "Q=M_{TB} - M_{T} - M_{B}",
                   "Q (GeV)", 50, 0, 500, "GeV", hQ,dir);
  
    
//T Mass Histos
  defineHistoSet("hTMass" , "Reconstructed Mass of T",
                 "M_{t} (GeV)", 80, 0, 400, "GeV", hTMass,dir);
  defineHistoSet("hTenuMass","Reconstructed Mass of Te#nu",
                 "M_{t}^{e#nu} (GeV)", 80, 0, 400, "GeV", hTenuMass,dir);
  defineHistoSet("hTmnuMass","Reconstructed Mass of T#mu#nu",
                 "M_{t}^{#mu#nu} (GeV)", 80, 0, 400, "GeV", hTmnuMass,dir);
  defineHistoSet("hEvtType", "Event Type",
                 "N_{#mu}", 2, 0, 2, "NONE", hEvtType,dir);
    
//Tpt Histos
  defineHistoSet("hTpt", "p_{T}^{t}", 
                 "p_{T}^{t} (GeV)", 50, 0, 1000, "GeV", hTpt,dir);
  defineHistoSet("hTenupt", "p_{T}^{t#rightarrowe#nu}", 
                 "p_{T}^{t#rightarrowbe#nu} (GeV)", 50, 0, 1000, "GeV", hTenupt,dir);
  defineHistoSet("hTmnupt", "p_{T}^{t#rightarrow#mu#nu}", 
                 "p_{T}^{t#rightarrowb#mu#nu} (GeV)", 50, 0, 1000, "GeV", hTmnupt,dir);
//MET Histos
  defineHistoSet("hMET", "MET",
                 "#slash{E}_{T} (GeV)", 50, 0, 500, "GeV", hMET,dir);
  defineHistoSet("hMETSig", "MET Significance",
                 "#slash{E}_{T}^{Signif}", 50, 0, 50, "NONE", hMETSig,dir);
    
//W Trans Mass Histos
  defineHistoSet("hWTransMass", "Reconstructed Transverse Mass of W",
                 "M_{T}^{W} (GeV)", 20, 0, 100, "GeV", hWTransMass,dir);
  defineHistoSet("hWenuTransMass", "Reconstructed Transverse Mass of We\\nu",
                 "M_{T}^{W#rightarrowe#nu} (GeV)", 20, 0, 100, "GeV", hWenuTransMass,dir);
  defineHistoSet("hWmnuTransMass", "Reconstructed TransverseMass of W#mu\\nu",
                 "M_{T}^{W#rightarrow#mu#nu} (GeV)", 20, 0, 100, "GeV", hWmnuTransMass,dir);
    
//Wpt Histos
  defineHistoSet("hWpt", "p_{T}^{W}", 
                 "p_{T}^{W} (GeV)", 50, 0, 1000, "GeV", hWpt,dir);
    
//B Histos
  defineHistoSet("hB1Mass" , "Reconstructed Mass of B1",
                 "M_{B1} (GeV)", 100, 0, 100, "GeV", hB1Mass,dir);
  defineHistoSet("hB1Disc" , "B1 Discriminator of B1",
                 "Disc_{B1}", 100, 0, 20, "NONE", hB1Disc,dir);
  defineHistoSet("hB1pt", "p_{T}^{B1}", 
                 "p_{B1}^{T} (GeV)", 50, 0, 1000, "GeV", hB1pt,dir);

  defineHistoSet("hB2Mass" , "Reconstructed Mass of B2",
                 "M_{B2} (GeV)", 100, 0, 100, "GeV", hB2Mass,dir);
  defineHistoSet("hB2Disc" , "B2 Discriminator of B2",
                 "Disc_{B2}", 100, 0, 20, "NONE", hB2Disc,dir);
  defineHistoSet("hB2pt", "p_{T}^{B2}", 
                 "p_{B2}^{T} (GeV)", 50, 0, 1000, "GeV", hB2pt,dir);

  //M_bl
  defineHistoSet("hMbl", "M_{bl}",
                 "M_{bl} (GeV)", 20, 0., 200., "GeV", hMbl, dir);


//Other histograms
  defineHistoSet("hNLElec", "Number of Loose Electrons in Event",
                 "N_{e}^{Loose}", 10, 0, 10, "NONE", hNLElec,dir);
  defineHistoSet("hNLMuon", "Number of Loose Muons in Event",
                 "N_{#mu}^{Loose}", 10, 0, 10, "NONE", hNLMuon,dir);
  defineHistoSet("hNLLeps", "Number of Loose Leptons in Event",
                 "N_{Lep}^{Loose}", 10, 0, 10, "NONE", hNLLeps,dir);
    
  defineHistoSet("hNLJets", "Number of Loose Jets in Event",
                 "N_{Jets}^{Loose}", 30, 0, 30, "NONE", hNLJets,dir);
  defineHistoSet("hNLBJets", "Number of Loose BJets in Event",
                 "N_{B Jets}^{Loose}", 30, 0, 30, "NONE", hNLBJets,dir);
  defineHistoSet("hNTBJets", "Number of Tight BJets in Event",
                 "N_{B Jets}^{Tight}", 30, 0, 30, "NONE", hNTBJets,dir);
    
  //defineHistoSet("hNVtxs", "Number of Vertexs in Event",
  //               "N_{Vtx}", 50, 0, 50, "NONE", hNVtxs,dir);

  defineHistoSet("hWeight", "PU Weight",
                 "Weight", 40, 0, 2, "NONE", hWeight,dir);

  tTBCand = dir.make<TTree>("tTBCand", "Analysis Variables after TBCand");//Only 1 for now;
  tTBCand->Branch("Run", &runNumber_);
  tTBCand->Branch("Lumi", &lumiNumber_);
  tTBCand->Branch("Event", &evtNumber_);
  tTBCand->Branch("TBMass", &TBMass_);
  tTBCand->Branch("EvtType", &evtType_);
  tTBCand->Branch("TMass", &TMass_);
  tTBCand->Branch("Tpt", &Tpt_);
  tTBCand->Branch("BMass1", &BMass1_);
  tTBCand->Branch("BDisc1", &BDisc1_);
  tTBCand->Branch("Bpt1", &Bpt1_);
  tTBCand->Branch("BMass2", &BMass2_);
  tTBCand->Branch("BDisc2", &BDisc2_);
  tTBCand->Branch("Bpt2", &Bpt2_);
  tTBCand->Branch("WTransMass", &WTransMass_);
  tTBCand->Branch("Wpt", &Wpt_);
  tTBCand->Branch("Q", &Q_);
  tTBCand->Branch("MET", &MET_);
  tTBCand->Branch("METSig", &METSig_);
  tTBCand->Branch("NVtxs", &NVtxs_);
  tTBCand->Branch("Mbl", &Mbl_);
  tTBCand->Branch("weight", &weight_);

}//defineHistos


void TBAnalyzer::fillHistos(const int& index, const float& weight){
  if(tbCand_ && wCand_ && bCand1_.mass()){
    hTBMass[index]->Fill(tbCand_().mass(), weight);
    if     (evtType_ == 0) hTBenuMass[index]->Fill(tbCand_().mass(), weight);
    else if(evtType_ == 1) hTBmnuMass[index]->Fill(tbCand_().mass(), weight);
  }
  if(tCand_){
    hTMass[index]->Fill(tCand_().mass(), weight);
    hTpt[index]->Fill(tCand_().pt(), weight);
    hEvtType[index]->Fill(evtType_, weight);
    hQ[index]->Fill(Q_, weight);
    hMbl[index]->Fill(Mbl_, weight);
    if      (wCand_.flavor() == PDG_ID_ELEC){
      hTenuMass[index]->Fill(tCand_().mass(), weight);
      hTenupt[index]->Fill(tCand_().pt(), weight);
    }else if (wCand_.flavor() == PDG_ID_MUON){
      hTmnuMass[index]->Fill(tCand_().mass(), weight);
      hTmnupt[index]->Fill(tCand_().pt(), weight);
    }
  }
  if(wCand_){
    hWTransMass[index]->Fill(wCand_.mt(), weight);
    hWpt[index]->Fill(wCand_.pt(), weight);
    if      (wCand_.flavor() == PDG_ID_ELEC){
      hWenuTransMass[index]->Fill(wCand_.mt(), weight);
    }else if (wCand_.flavor() == PDG_ID_MUON){
      hWmnuTransMass[index]->Fill(wCand_.mt(), weight);
    }
  }  
  if(bCand1_.mass()){
    hB1Mass[index]->Fill(bCand1_.mass(), weight);
    hB1Disc[index]->Fill(bCand1_.bDiscriminator(BDisc_), weight);
    hB1pt[index]->Fill(bCand1_.pt(), weight);
  }
  if(bCand2_.mass()){
    hB2Mass[index]->Fill(bCand2_.mass(), weight);
    hB2Disc[index]->Fill(bCand2_.bDiscriminator(BDisc_), weight);
    hB2pt[index]->Fill(bCand2_.pt(), weight);
  }

  hMET[index]->Fill(met_.et(), weight);
  hMETSig[index]->Fill(met_.significance(), weight);
    
  hNLElec[index]->Fill(looseElectrons_.size(), weight);
  hNLMuon[index]->Fill(looseMuons_    .size(), weight);
  hNLLeps[index]->Fill(looseElectrons_.size()+looseMuons_.size(), weight);
    
  hNLJets[index]->Fill(looseJets_.size(), weight);
  hNLBJets[index]->Fill(looseBJets_.size(), weight);
  hNTBJets[index]->Fill(tightBJets_.size(), weight);

  //hNVtxs[index]->Fill((*verticesH_).size(), weight);
  hWeight[index]->Fill(weight_/wprimeUtil_->getSampleWeight(), 1.);//Don't weight

  if(CutNames_[index] == "ValidTB") tTBCand->Fill();

}//fillHistos


void 
TBAnalyzer::eventLoop(edm::EventBase const & event){
  clearEvtVariables();
  runNumber_ = event.id().run();
  lumiNumber_ = event.id().luminosityBlock();
  evtNumber_ = event.id().event();

  if(debug_) WPrimeUtil::printEvent(event, cout);
  if (doPreselect_){
    if(debug_) cout<<"Testing Preselection...\n";
  }
  
  event.getByLabel(electronsLabel_,patElectronsH_);
  event.getByLabel(muonsLabel_,patMuonsH_);
  event.getByLabel(metLabel_, metH_);
  met_ = (*metH_)[0];
  MET_ = met_.et();
  METSig_ = met_.significance();
  //NVtxs_ = (*verticesH_).size();

  if(patElectronsH_->size() + patMuonsH_->size() == 0){
    cout << "Not enough leptons. Bad bad event, returning now..." << endl;
    return;
  }
  for (size_t i = 0; i < patElectronsH_->size(); i++) {
    allElectrons_.push_back(heep::Ele((*patElectronsH_)[i]));   
    if(Overlap(allElectrons_[i].patEle(), *patMuonsH_.product(), 0.01)) continue;
    if (looseElectron_(allElectrons_[i].patEle()))
      looseElectrons_.push_back(allElectrons_[i]);
  }

  // Loop over muons, and see if they pass the TeVMuon criteria  
  for (size_t i = 0; i < patMuonsH_->size(); i++) {
    allMuons_.push_back(TeVMuon((*patMuonsH_)[i],muReconstructor_));   
    
    if (looseMuon_(allMuons_[i]) )
      looseMuons_.push_back(allMuons_[i]);
  }
  if(debug_){
    print(allElectrons_);
    print(allMuons_);
    printf("    Contains: %i electron(s), %i muon(s)\n",
           (int)allElectrons_.size(), (int)allMuons_.size());
    printf("    Contains: %i loose electron(s), %i loose muon(s)\n",
           (int)looseElectrons_.size(), (int)looseMuons_.size());
  }

  //get Jets
  event.getByLabel(jetsLabel_, patJetsH_);
  const JetV & allJets  = *patJetsH_;
  if(allJets.size() < 1){
    cout << "Not enough jets. Bad bad event, returning now..." << endl;
    return;
  }
  if(debug_){
    print(allJets);
    printf("    Contains: %i pat jet(s)\n",
           (int)allJets.size());
  }

  // Loop over jets, and see if they pass the jet criteria
  for (size_t i = 0; i < allJets.size(); ++i) {
    if (looseJet_(allJets[i]) 
        && !Overlap(allJets[i], looseElectrons_, 0.5)
        && !Overlap(allJets[i], looseMuons_    , 0.5)){
      looseJets_.push_back(allJets[i]);

      if(allJets[i].bDiscriminator(BDisc_) > minBDisc_){//cory: these should be loose jets too!!!!
        looseBJets_.push_back(allJets[i]);

        if(allJets[i].mass() < maxBMass_)
          tightBJets_.push_back(allJets[i]);
      }
    }
  }
  if(debug_){
    printf("    Contains: %i loose jet(s)\n",
           (int)looseJets_.size());
    printf("    Contains: %i loose B jet(s)\n",
           (int)looseBJets_.size());
    printf("    Contains: %i tight B jet(s)\n",
           (int)tightBJets_.size());
  }
  
  //////////////////////////////////////////////
  /// Start Applying Cuts///////////////////////
  //////////////////////////////////////////////
  weight_ = wprimeUtil_->getWeight();

  int iCut=0;
  if( !passNoCut() ) return;
  tabulateEvent(iCut++, weight_);

  if( !passMinNLeptonsCut(looseElectrons_, looseMuons_, minNLeptons_) ) return;
  tabulateEvent(iCut++, weight_);

  if( !passMinNJetsCut(looseJets_, minNJets_) ) return;
  tabulateEvent(iCut++, weight_);

  if( !passMinNJetsCut(looseBJets_, minNBJets_) ) return;
  tabulateEvent(iCut++, weight_);

  if( !passMinMETCut(met_, minMET_) ) return;
  tabulateEvent(iCut++, weight_);

  //////////////////////////////
  ///////  Make B Jet  /////////
  //////////////////////////////
  
  //TODO: Implement here
  /*
  for(uint i=0; i<looseJets_.size(); ++i){
    const vector<pair<string,float> > & discPair = looseJets_[i].getPairDiscri ();
    cout<<" I have this many btags "<<discPair.size()<<endl;
    for(uint j=0; j<discPair.size(); ++j){
      cout<<discPair[j].first<< " " << discPair[j].second << endl;
    }
  }
  */
  //if( !passValidBCut(bCand1_) ) return;
  //tabulateEvent(iCut++, weight_);

  //////////////////////////////
  ///////  Make Ws /////////////
  //////////////////////////////

  WCandV wCands = getWCandidates(looseElectrons_, looseMuons_, met_);

  //Clean up w's
  for (WCandV::iterator w = wCands.begin(); w != wCands.end(); ++w){
    if(!passWtransMassCut(*w, minWtransMass_) ||
       !passWptCut(*w, minWpt_)){
      wCands.erase(w);
      w--;
    }
  }

  if(!wCands.size()) return;
  tabulateEvent(iCut++, weight_);  

  //////////////////////////////
  /////  Make Top Quark  ///////
  //////////////////////////////

  //TODO: Implement here
  const JetV & bCands = looseBJets_;
  int iBestT(-1), iBestW(0), iBestB(0);
  float bestTMass = 9e9;
  XWLepV tCands;
  const int nW = wCands.size();
  const int nB = bCands.size();
  if(debug_){
    cout<<"There are "<<nW<<"W candidates"<<endl;
    cout<<"There are "<<nB<<"B candidates"<<endl;
  }
  for(int iw=0; iw<nW; ++iw){
    for(int ib=0; ib<nB; ++ib){
      //Skip if m_bl is greater than theoretical limit
      float m_bl = (bCands[ib].p4() + wCands[iw].daughter(0)->p4()).M();
      if(debug_) cout<<"m_bl: "<<m_bl<<" MAX m_bl: "<<MAX_MBL<<endl;
      if(m_bl > MAX_MBL)  continue;

      tCands.push_back(XWLeptonic(bCands[ib], wCands[iw]));
      int it = tCands.size() - 1;//iw*nB + ib;
      if(debug_) 
        cout<<"iw:ib:Best:curr = "<<iw<<":"<<ib<<"\t"
            <<bestTMass<<"\t"<<tCands[it]().mass()<<endl;
      
      if(fabs(tCands[it]().mass() - TMASS) < fabs(bestTMass - TMASS)){
        if(debug_) cout<<"Found a new best mass"<<endl;
        iBestT = it; iBestW = iw; iBestB = ib;
        bestTMass = tCands[it]().mass();
      }
    }
  }

  tCand_ = iBestT >= 0 ? tCands[iBestT] : XWLeptonic();
  wCand_ = wCands[iBestW];
  bCand1_= bCands[iBestB];

  if( !passValidTCut(tCand_) ) return;
  tabulateEvent(iCut++, weight_);
  
  evtType_ = wCand_ && wCand_.flavor()==PDG_ID_MUON;
  TMass_ = tCand_(nuAlgo_).mass();
  Tpt_ = tCand_(nuAlgo_).pt();
  BMass1_ = bCand1_.mass();
  BDisc1_ = bCand1_.bDiscriminator(BDisc_);
  Bpt1_ = bCand1_.pt();
  WTransMass_ = wCand_.mt();
  Wpt_  = wCand_.pt();
  Mbl_ = (bCand1_.p4() + wCand_.daughter(0)->p4()).M();

  if( !passXWMassCut(tCand_, nuAlgo_, minTMass_, maxTMass_) )return;
  tabulateEvent(iCut++, weight_);

  if( !passXWptCut(tCand_, minTpt_) ) return;
  tabulateEvent(iCut++, weight_);

  //////////////////////////////
  ///////  Make 2nd B  /////////
  //////////////////////////////

  //Take highest bjet that's not used.
  uint bestbjet2 = 0;
  for(uint ib=0; ib<looseJets_.size(); ++ib){
    if(areIdentical(looseJets_[ib], bCand1_)) continue;
    if(looseJets_[ib].bDiscriminator(BDisc_) > looseJets_[bestbjet2].bDiscriminator(BDisc_))
      bestbjet2 = ib;
  }
  bCand2_ = looseJets_[bestbjet2];

  BMass2_ = bCand2_.mass();
  BDisc2_ = bCand2_.bDiscriminator(BDisc_);
  Bpt2_ = bCand2_.pt();

  if( !passValidBCut(bCand2_) ) return;
  tabulateEvent(iCut++, weight_);

  if( !passMinPtCut(bCand2_, minBpt_) ) return;
  tabulateEvent(iCut++, weight_);

  //////////////////////////////
  ///////  Make TB Cand  ///////
  //////////////////////////////

  tbCand_ = XWLeptonic(bCand2_, tCand_);
  TBMass_ = tbCand_(nuAlgo_).mass();
  Q_ = TBMass_ - TMass_ - BMass2_;

  if( !passValidXWCut(tbCand_) ) return;
  tabulateEvent(iCut++, weight_);

  //////////////////////////////
  ///////  Analysis Cuts ///////
  //////////////////////////////



  //////////////////////////////
  ///////  Done wih Cuts ///////
  //////////////////////////////

  //AllCuts
  tabulateEvent(iCut++, weight_);

  if(wprimeUtil_->runningOnData()){
    cout<<" The following data event passed All Cuts!!!\n";
    printPassingEvent(event);
    cout<<" ------------------\n";
  }
  
}

inline void
TBAnalyzer::clearEvtVariables(){
  looseJets_.clear();
  looseBJets_.clear();
  tightBJets_.clear();
  allElectrons_.clear();
  looseElectrons_.clear();
  allMuons_.clear();
  looseMuons_.clear();
  met_ = pat::MET();
  bCand1_ = pat::Jet();
  bCand2_ = pat::Jet();
  wCand_ = WCandidate();
  tCand_ = XWLeptonic();
  tbCand_ = XWLeptonic();
  runNumber_ = 0;
  lumiNumber_ = 0;
  evtNumber_ = 0;
  TBMass_ = -999;
  evtType_ = -999;
  TMass_ = Tpt_ = -999;
  BMass1_ = BDisc1_ = Bpt1_ = -999;
  BMass2_ = BDisc2_ = Bpt2_ = -999;
  Wpt_ = WTransMass_ = -999;
  Q_ = -999;
  MET_ = METSig_ = -999;
  NVtxs_ = -999;
  Mbl_ = -999;
  weight_ = 0;
}

/////////////////Cuts///////////////////////
inline bool TBAnalyzer::passValidBCut(pat::Jet& b){
  //calcBVariables();
  return b.mass()>0.;
}

inline bool TBAnalyzer::passValidTCut(XWLeptonic & t){
  //calcTVariables();
  return AnalyzerBase::passValidXWCut(t);
}

/////////Check TB Properties/////
inline bool TBAnalyzer::passValidTBCut(XWLeptonic & tb){
  //calcTBVariables();
  return AnalyzerBase::passValidXWCut(tb);
}
