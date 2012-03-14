#include "UserCode/CMGWPrimeGroup/interface/TTbarAnalyzer.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember)) 

using namespace std;
TTbarAnalyzer::TTbarAnalyzer(){}
TTbarAnalyzer::TTbarAnalyzer(const edm::ParameterSet & cfg, int fileToRun) :
  AnalyzerBase(cfg, fileToRun){
  //setupCutOrder();
  if(debug_) printf("Using %i cuts\n",NCuts_);

// +++++++++++++++++++Event characteristics
 
// +++++++++++++++++++General Cut values
  minNLeptons_ = cfg.getUntrackedParameter<uint>("minNLeptons", 0);
  minNJets_ = cfg.getUntrackedParameter<uint>("minNJets", 0);
  
// +++++++++++++++++++Z Cuts
  minZmass_ = cfg.getUntrackedParameter<double>("minZmass", 0.);
  maxZmass_ = cfg.getUntrackedParameter<double>("maxZmass", 9e9);
  minZpt_ = cfg.getUntrackedParameter<double>("minZpt", 0.);

// +++++++++++++++++++V Cuts
  minVmass_ = cfg.getUntrackedParameter<double>("minVmass", 0.);
  maxVmass_ = cfg.getUntrackedParameter<double>("maxVmass", 9e9);
  minVpt_ = cfg.getUntrackedParameter<double>("minVpt", 0.);

  maxAngleBetweenJets = cfg.getParameter<double>("maxAngleBetweenJets");


// +++++++++++++++++++Hadronic Boson Cuts


  genLabel_ = cfg.getParameter<edm::InputTag>("genParticles" );

  //Selectors
  Pset eSelectorPset = cfg.getParameter<Pset>("electronSelectors");
  string looseElectronType = cfg.getUntrackedParameter<string>("LooseElectronType", "wp95");
  string tightElectronType = cfg.getUntrackedParameter<string>("TightElectronType", "wp95");
  looseElectron_ = ElectronSelector(eSelectorPset, looseElectronType);
  tightElectron_ = ElectronSelector(eSelectorPset, tightElectronType);
  if(debug_) cout<<"Using "<<looseElectronType<<" for loose electrons and "
                  <<tightElectronType<<" for tight electrons\n";

  Pset mSelectorPset = cfg.getParameter<Pset>("muonSelectors");
  string looseMuonType = cfg.getUntrackedParameter<string>("LooseMuonType", "exotica");
  string tightMuonType = cfg.getUntrackedParameter<string>("TightMuonType", "exotica");
  looseMuon_ = MuonSelector(mSelectorPset, looseMuonType);
  tightMuon_ = MuonSelector(mSelectorPset, tightMuonType);
  if(debug_) cout<<"Using "<<looseMuonType<<" for loose muons and "
                  <<tightMuonType<<" for tight muons\n";

  Pset jSelectorPset = cfg.getParameter<Pset>("jetSelectors");
  string looseJetType = cfg.getUntrackedParameter<string>("LooseJetType", "Base");
  looseJet_ = JetSelector(jSelectorPset, looseJetType);
  if(debug_) cout<<"Using "<<looseJetType<<" for jets\n";

}

TTbarAnalyzer::~TTbarAnalyzer(){
}


/// Declare Histograms
void TTbarAnalyzer::defineHistos(const TFileDirectory & dir){
  printf("Declare histos\n");

  /////////////////
  /////////////TTbar
//T Mass Histos
  defineHistoSet("hTMass" , "Reconstructed Mass of T",
                 "M_{t} (GeV)", 80, 0, 400, "GeV", hTMass,dir);
  defineHistoSet("hTenuMass","Reconstructed Mass of Te#nu",
                 "M_{t}^{e#nu} (GeV)", 80, 0, 400, "GeV", hTenuMass,dir);
  defineHistoSet("hTmnuMass","Reconstructed Mass of T#mu#nu",
                 "M_{t}^{#mu#nu} (GeV)", 80, 0, 400, "GeV", hTmnuMass,dir);
  defineHistoSet("hEvtType", "Event Type",
                 "N_{#mu}", 2, 0, 2, "NONE", hEvtType,dir);
//T2 Mass Histos
  defineHistoSet("hT2Mass" , "Reconstructed Mass of T2",
                 "M_{t}^{Had} (GeV)", 80, 0, 400, "GeV", hT2Mass,dir);
    
//Tpt Histos
  defineHistoSet("hTpt", "p_{T}^{t}", 
                 "p_{T}^{t} (GeV)", 50, 0, 1000, "GeV", hTpt,dir);
  defineHistoSet("hTenupt", "p_{T}^{t#rightarrowe#nu}", 
                 "p_{T}^{t#rightarrowbe#nu} (GeV)", 50, 0, 1000, "GeV", hTenupt,dir);
  defineHistoSet("hTmnupt", "p_{T}^{t#rightarrow#mu#nu}", 
                 "p_{T}^{t#rightarrowb#mu#nu} (GeV)", 50, 0, 1000, "GeV", hTmnupt,dir);
//T2pt Histos
  defineHistoSet("hT2pt", "p_{T}^{t had}", 
                 "p_{T}^{t had} (GeV)", 50, 0, 1000, "GeV", hT2pt,dir);
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
//W2 Trans Mass Histos
  defineHistoSet("hW2TransMass", "Reconstructed Transverse Mass of W had",
                 "M_{T}^{W had} (GeV)", 20, 0, 100, "GeV", hW2TransMass,dir);
//W2 Mass Histos
  defineHistoSet("hW2Mass", "Reconstructed Mass of W had",
                 "M^{W had} (GeV)", 30, 0, 150, "GeV", hW2Mass,dir);
    
//Wpt Histos
  defineHistoSet("hWpt", "p_{T}^{W}", 
                 "p_{T}^{W} (GeV)", 50, 0, 1000, "GeV", hWpt,dir);
//W2pt Histos
  defineHistoSet("hW2pt", "p_{T}^{W had}", 
                 "p_{T}^{W had} (GeV)", 50, 0, 1000, "GeV", hW2pt,dir);
    
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
    
  defineHistoSet("hNVtxs", "Number of Vertexs in Event",
                 "N_{Vtx}", 50, 0, 50, "NONE", hNVtxs,dir);

  defineHistoSet("hWeight", "PU Weight",
                 "Weight", 40, 0, 2, "NONE", hWeight,dir);

  ///////////////
  hW2PtVsW2Mass = dir.make<TH2F>("hW2PtVsW2Mass", "hW2PtVsW2Mass; M_{W} (GeV); p_{T}^{W} (GeV)", 20, 0., 100., 50, 0., 1000.);

  cout << "Histos declared" << endl;



}//defineHistos


void 
TTbarAnalyzer::eventLoop(edm::EventBase const & event){
  clearEvtVariables();
  runNumber_ = event.id().run();
  lumiNumber_ = event.id().luminosityBlock();
  evtNumber_ = event.id().event();
  if(debug_){
    cout<<" New event: ";
    WPrimeUtil::printEvent(event);
  }
  weight_ = wprimeUtil_->getWeight();


  // Preselection - skip events that don't look promising
  if (doPreselect_){
    if(debug_) cout<<"Testing Preselection...\n";
    // We could setup some preselection here. To be implemented.
  }
  

  //////////////////////
  //Deal With Leptons///
  //////////////////////
  event.getByLabel(electronsLabel_,patElectronsH_);
  event.getByLabel(muonsLabel_,patMuonsH_);
  if(debug_){
    printf("    Contains: %i pat electron(s), %i pat muon(s)\n",
           (int)patElectronsH_->size(), (int)patMuonsH_->size());
  }
  if(patElectronsH_->size() + patMuonsH_->size() == 0){
    cout << "Not enough leptons. Bad bad event, returning now..." << endl;
    return;
  }
  
  // Make vectors of leptons passing various criteria
  // Loop over muons, and see if they pass the TeVMuon criteria  
  for (size_t i = 0; i < patMuonsH_->size(); i++) {
    allMuons_.push_back(TeVMuon((*patMuonsH_)[i],muReconstructor_));   
    
    if (looseMuon_(allMuons_[i]) )
      looseMuons_.push_back(allMuons_[i]);
    
    if (tightMuon_(allMuons_[i]) )
      tightMuons_.push_back(allMuons_[i]);
  }
  // Loop over electrons, and see if they pass the criteria
  for (size_t i = 0; i < patElectronsH_->size(); i++) {
    allElectrons_.push_back(heep::Ele((*patElectronsH_)[i]));   
    
    if(Overlap(allElectrons_[i].patEle(), *patMuonsH_.product(), 0.01)) continue;
    
    if (looseElectron_(allElectrons_[i].patEle()))
      looseElectrons_.push_back(allElectrons_[i]);

    if (tightElectron_(allElectrons_[i].patEle()))
      tightElectrons_.push_back(allElectrons_[i]);
  }


  if(debug_){
    printf("    Contains: %i electron(s), %i muon(s)\n",
           (int)allElectrons_.size(), (int)allMuons_.size());
    printf("    Contains: %i loose electron(s), %i loose muon(s)\n",
           (int)looseElectrons_.size(), (int)looseMuons_.size());
    printf("    Contains: %i tight electron(s), %i tightmuon(s)\n",
           (int)tightElectrons_.size(), (int)tightMuons_.size());
    printLeptons();
  }


  // if(useAdjustedMET_) event.getByLabel(pfCandsLabel_, pfCandidatesH_);
  //Deal with MET
  event.getByLabel(metLabel_, metH_);
  WPrimeUtil::getLeptonsMET(patElectronsH_, allElectrons_,
                            patMuonsH_, muReconstructor_, allMuons_,
                            metH_, useAdjustedMET_, met_,
                            pfCandidatesH_);

  if (looseMuons_.size() > 0){
    sort(looseMuons_.begin(), looseMuons_.end(), highestMuonPt());
    if(debug_){
      for (uint k=0; k<looseMuons_.size(); k++){
	      cout << "looseMuons_ " << k << " has pt of " << looseMuons_.at(k).pt() <<  endl;
	    }
    }
  }

  if (tightMuons_.size() > 0){
    sort(tightMuons_.begin(), tightMuons_.end(), highestMuonPt());
  }

  if (looseElectrons_.size() > 0){
    sort(looseElectrons_.begin(), looseElectrons_.end(), highestElectronPt());
    if(debug_){
      for (uint k=0; k<looseElectrons_.size(); k++){
	cout << "looseElectrons_ " << k << " has pt of " << looseElectrons_.at(k).patEle().pt() <<  endl;
	    }
    }
  }
  if (tightElectrons_.size() > 0){
    sort(tightElectrons_.begin(), tightElectrons_.end(), highestElectronPt());
  }

  //////////////////////
  ////Deal With Jets////
  //////////////////////

  allJets_      = getProduct<vector<pat::Jet     > >(event, jetsLabel_);
  if(allJets_.size() < 1){
    if (debug_) 
      cout << "Not enough jets. Bad bad event, returning now..." << endl;
    return;
  }
  if(debug_)
    printf("    Contains: %i pat jet(s)\n",
           (int)allJets_.size());


  // Loop over jets, and see if they pass the jet criteria
  for (size_t i = 0; i < allJets_.size(); ++i) {
    if (looseJet_(allJets_[i]) && !Overlap(allJets_[i], looseMuons_, 1.0, 2) && !Overlap(allJets_[i], looseElectrons_, 1.0, 2)){
      looseJets_.push_back(allJets_[i]);
      if(allJets_[i].bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 0.679){
        looseBJets_.push_back(allJets_[i]);
      }
    }
  }
  
  if (debug_){
    printf("    Contains: %i loose jets(s)\n",
           (int)looseJets_.size());
    printf("    Contains: %i loose B jets(s)\n",
           (int)looseBJets_.size());
  }

  //get Vertex
  event.getByLabel(vertexLabel_, verticesH_);


  //////////////////////////////////////////////
  /// Start Applying Cuts///////////////////////
  //////////////////////////////////////////////

  vCand_ = ZCandidate();
  wCand_ = WCandidate();
  hadVZ_ = VZCandidate();
  int iCut=0;
  if( !passNoCut() ) return;
  tabulateEvent(iCut, weight_); ++iCut;

  if( !passMinNLeptonsCut(looseElectrons_, looseMuons_, 1) ) return;
  tabulateEvent(iCut, weight_); ++iCut;

  if( !passMinNJetsCut(looseJets_, 3) ) return;
  tabulateEvent(iCut, weight_); ++iCut;

  if( !passMinNJetsCut(looseBJets_, 1) ) return;
  tabulateEvent(iCut++, weight_);

  if( !passMinMETCut(met_, 30.) ) return;
  tabulateEvent(iCut++, weight_);

  //////////////////////////////
  ///////  Make Ws /////////////
  //////////////////////////////

  WCandV wCands = getWCandidates(looseElectrons_, looseMuons_, met_);

  //Clean up w's
  for (WCandV::iterator w = wCands.begin(); w != wCands.end(); ++w){
    if(!passWtransMassCut(*w, 30.) ||
       !passWptCut(*w, 20.)){
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
      if(m_bl > MAX_MBL*1.2)  continue;

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

  if( !passValidXWCut(tCand_) ) return;
  tabulateEvent(iCut++, weight_);
  /*
  evtType_ = wCand_ && wCand_.flavor()==PDG_ID_MUON;
  TMass_ = tCand_(kMinPz).mass();
  Tpt_ = tCand_(kMinPz).pt();
  BMass1_ = bCand1_.mass();
  BDisc1_ = bCand1_.bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
  Bpt1_ = bCand1_.pt();
  WTransMass_ = wCand_.mt();
  Wpt_  = wCand_.pt();
  */
  Mbl_ = (bCand1_.p4() + wCand_.daughter(0)->p4()).M();
  
  if( !passXWMassCut(tCand_, kMinPz, 150, 200) )return;
  tabulateEvent(iCut++, weight_);

  //////////////////////////////
  ///////  Make 2nd B  /////////
  //////////////////////////////

  //Take highest bjet that's not used.
  uint bestbjet2 = 0;
  for(uint ib=0; ib<looseJets_.size(); ++ib){
    if(areIdentical(looseJets_[ib], bCand1_)) continue;
    if(looseJets_[ib].bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > looseJets_[bestbjet2].bDiscriminator("simpleSecondaryVertexHighEffBJetTags"))
      bestbjet2 = ib;
  }
  bCand2_ = looseJets_[bestbjet2];
  /*
  BMass2_ = bCand2_.mass();
  BDisc2_ = bCand2_.bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
  Bpt2_ = bCand2_.pt();
  */
  if( 0 ) return;
  tabulateEvent(iCut++, weight_);

  if( !passMinPtCut(bCand2_, 30.) ) return;
  tabulateEvent(iCut++, weight_);

  //////////////////////////////
  ///////  Make 2nd Top ////////
  //////////////////////////////
  /*
  //ignore the class name, this is really a Hadronic Top
  vector<VZCandidate> hadTops;
  bestTMass = 9e9;
  iBestT = -1;
  int iBestJ(0);
  for(uint ij=0; ij<looseJets_.size(); ++ij){
    if(areIdentical(looseJets_[ij], bCand1_)) continue;
    if(areIdentical(looseJets_[ij], bCand2_)) continue;

    hadTops.push_back(VZCandidate(bCand2_, looseJets_[ij]));
    int it = hadTops.size() - 1;
    if(debug_) 
      cout<<"ij:Best:curr = "<<ij<<"\t"
          <<bestTMass<<"\t"<<hadTops[it].mass()<<endl;
    
    if(fabs(hadTops[it].mass() - TMASS) < fabs(bestTMass - TMASS)){
      if(debug_) cout<<"Found a new best mass"<<endl;
        iBestT = it; iBestJ = ij; 
        bestTMass = hadTops[it].mass();
      }

  }
  wJet_ = looseJets_[iBestJ];
  hadTop_ = iBestT >= 0 ? hadTops[iBestT] : VZCandidate();
  */

  ////////
  size_t bestPos = 0;
  double maxPt = -1.0;
  for(uint ij=0; ij<looseJets_.size(); ++ij){
    if(areIdentical(looseJets_[ij], bCand1_)) continue;
    if(areIdentical(looseJets_[ij], bCand2_)) continue;
    if(looseJets_[ij].pt() > maxPt) {
      maxPt = looseJets_[ij].pt();
      bestPos = ij;
    }
  }
  wJet_ = looseJets_[bestPos];
  hadTop_ = VZCandidate(bCand2_, looseJets_[bestPos]);
  
  hW2PtVsW2Mass->Fill(wJet_.mass(), wJet_.pt(), weight_);

  if( 0 ) return;
  tabulateEvent(iCut++, weight_);

  if( 150 > wJet_.pt() ) return;
  tabulateEvent(iCut++, weight_);

  if( 150 > hadTop_.mass() || hadTop_.mass() > 200. ) return;
  tabulateEvent(iCut++, weight_);

  if( 70 > wJet_.mass() || wJet_.mass() > 120. ) return;
  tabulateEvent(iCut++, weight_);

  if( 200 > wJet_.pt() ) return;
  tabulateEvent(iCut++, weight_);

  if( 250 > wJet_.pt() ) return;
  tabulateEvent(iCut++, weight_);

  //////////////////////////////
  ///////  Wrapping Up  ////////
  //////////////////////////////

  //AllCuts
  tabulateEvent(iCut, weight_); ++iCut;

  if(wprimeUtil_->runningOnData()){
    cout<<" The following data event passed All Cuts!!!\n";
    printPassingEvent(event);
    if(1 || debug_){ 
      //printEventLeptons();
      printElectrons();
      printMuons();
      printJets();
    }
    cout<<" ------------------\n";
  }
  
}//End of Event Loop

/////////////////Accessors///////////////////////

/////////////////Modifies///////////////////////
  
/////////////////Cuts///////////////////////

bool
TTbarAnalyzer::passValidVZCandCut(){
  return hadVZ_ && hadVZ_.mass()>0.;
}


////////////////////////////////
/////////Check Z Properties/////
////////////////////////////////

////////////////////////////////
///////Check Had. V Properties//
////////////////////////////////

////////////////////////////////
//////Check Muon Properties/////
////////////////////////////////

///////////////Utilities//////////////////
void TTbarAnalyzer::printEventDetails() const{
  if(zCand_){
    cout<<" Z Flavor: "<<zCand_.flavor()
        <<" Z Mass: "<<zCand_.mass()
        <<" Z Eta: "<<zCand_.eta()
        <<" Z Phi: "<<zCand_.phi()
        <<endl;
  }
  if(vCand_){
    cout<<" V Mass: "<<vCand_.mass()
        <<" V Eta: "<<vCand_.eta()
        <<" V Phi: "<<vCand_.phi()
        <<endl;
  }
  if(zCand_ && vCand_ && hadVZ_.mass()>0.){
    cout<<" VZ Mass: "<<hadVZ_.mass()
      <<" VZ Eta: "<<hadVZ_.eta()
      <<" VZ Phi: "<<hadVZ_.phi()
        <<" VZpt: "<<hadVZ_.pt()
        <<" Zpt: "<<zCand_.pt()
        <<" Vpt: "<<vCand_.pt()
        <<endl;
  }
}

void
TTbarAnalyzer::clearEvtVariables(){
  AnalyzerBase::clearEvtVariables();
  zCand_ = ZCandidate();
  vCand_ = ZCandidate();
  hadVZ_ = VZCandidate();
  VZMass_ = -999;
  evtType_ = -999;
  ZMass_ = -999;
  Zpt_ = -999;
  VMass_=-999;
  Vpt_ = -999;
  Q_ = -999;
  DeltaRll_ = DeltaRVZ_ = -999;
  weight_ = 0;
  genMuons.clear();
  ak7GenJet.clear();
  looseBJets_.clear();
  wJet_ = pat::Jet();
  hadTop_ = VZCandidate();

}

//fill Histograms
void TTbarAnalyzer::fillHistos(const int& index, const float& weight){
  if(debug_) printf("filling Histos\n");
  if(tCand_){
    hTMass[index]->Fill(tCand_().mass(), weight);
    hTpt[index]->Fill(tCand_().pt(), weight);
    //hEvtType[index]->Fill(evtType_, weight);
    //hQ[index]->Fill(Q_, weight);
    hMbl[index]->Fill(Mbl_, weight);
    if      (wCand_.flavor() == PDG_ID_ELEC){
      hTenuMass[index]->Fill(tCand_().mass(), weight);
      hTenupt[index]->Fill(tCand_().pt(), weight);
    }else if (wCand_.flavor() == PDG_ID_MUON){
      hTmnuMass[index]->Fill(tCand_().mass(), weight);
      hTmnupt[index]->Fill(tCand_().pt(), weight);
    }
  }
  if(hadTop_){
    hT2Mass[index]->Fill(hadTop_.mass(), weight);
    hT2pt[index]->Fill(hadTop_.pt(), weight);
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
  if(wJet_.mass()){
    hW2Mass[index]->Fill(wJet_.mass(), weight);
    hW2TransMass[index]->Fill(wJet_.mt(), weight);
    hW2pt[index]->Fill(wJet_.pt(), weight);
  }  
  if(bCand1_.mass()){
    hB1Mass[index]->Fill(bCand1_.mass(), weight);
    hB1Disc[index]->Fill(bCand1_.bDiscriminator("simpleSecondaryVertexHighEffBJetTags"), weight);
    hB1pt[index]->Fill(bCand1_.pt(), weight);
  }
  if(bCand2_.mass()){
    hB2Mass[index]->Fill(bCand2_.mass(), weight);
    hB2Disc[index]->Fill(bCand2_.bDiscriminator("simpleSecondaryVertexHighEffBJetTags"), weight);
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

  hNVtxs[index]->Fill((*verticesH_).size(), weight);
  hWeight[index]->Fill(weight_/wprimeUtil_->getSampleWeight(), 1.);//Don't weight
  
}//fillHistos

void TTbarAnalyzer::fillJetMultiplicityHists(){
  h_jet_mult->Fill(looseJets_.size(), weight_);
  uint max = min(10, (int)looseJets_.size());
  for(uint i=0; i<=max; ++i)
    h_jet_mult_inc->Fill(i, weight_);
}


void TTbarAnalyzer::fillGoodZHistos(){
  if     (zCand_.flavor() == PDG_ID_ELEC){
    const heep::Ele & e1 = WPrimeUtil::Find(*zCand_.daughter(0), allElectrons_);
    const heep::Ele & e2 = WPrimeUtil::Find(*zCand_.daughter(1), allElectrons_);
    if (debug_)
      cout << "Found my electrons from loose Z" << endl;
    h_Zelec1_pt->Fill(e1.patEle().pt(), weight_);
    h_Zelec1_eta->Fill(e1.eta(), weight_);
    h_Zelec1_phi->Fill(e1.phi(), weight_);
    h_Zelec2_pt->Fill(e2.patEle().pt(), weight_);
    h_Zelec2_eta->Fill(e2.eta(), weight_);
    h_Zelec2_phi->Fill(e2.phi(), weight_);	  
    h_deltaR_elec1elec2->Fill(reco::deltaR(e1, e2), weight_);
  }else if(zCand_.flavor() == PDG_ID_MUON){
    const TeVMuon & m1 = WPrimeUtil::Find(*zCand_.daughter(0), allMuons_);
    const TeVMuon & m2 = WPrimeUtil::Find(*zCand_.daughter(1), allMuons_);
    if (debug_)
      cout << "Found my muons from loose Z" << endl;
    h_Zmuon1_pt->Fill(m1.pt(), weight_);
    h_Zmuon1_eta->Fill(m1.eta(), weight_);
    h_Zmuon1_phi->Fill(m1.phi(), weight_);
    h_Zmuon2_pt->Fill(m2.pt(), weight_);
    h_Zmuon2_eta->Fill(m2.eta(), weight_);
    h_Zmuon2_phi->Fill(m2.phi(), weight_);	  
    h_deltaR_muon1muon2->Fill(reco::deltaR(m1, m2), weight_);
  }
}

void TTbarAnalyzer::fillGoodHadVHistos(){
  h_jet_HadV_pt->Fill(vCand_.pt(), weight_);
  h_jet_HadV_eta->Fill(vCand_.eta(), weight_);
  h_jet_HadV_phi->Fill(vCand_.phi(), weight_);
  if (debug_)
    cout << "filled my HadV histos" << endl;
}

void TTbarAnalyzer::fillValidVZHistos(){
  h_HadVZMass->Fill(hadVZ_.mass()/1000, weight_);
  h_VZMass_nJets->Fill(looseJets_.size(), hadVZ_.mass());

  h_HadVZpt->Fill(hadVZ_.pt(), weight_);
  h_HadVZeta->Fill(hadVZ_.eta(), weight_);
  h_HadVZphi->Fill(hadVZ_.phi(), weight_);
  if (gravMass_>0)
    h_HadVZ_res->Fill((hadVZ_.mass()-gravMass_), weight_);
  if (debug_)
    cout << "filled my histos from HadVZ" << endl;
  if     (zCand_.flavor() == PDG_ID_ELEC){
    const heep::Ele & VZe1 = WPrimeUtil::Find(*zCand_.daughter(0), allElectrons_);
    const heep::Ele & VZe2 = WPrimeUtil::Find(*zCand_.daughter(1), allElectrons_);
    //cout << "Electron from loose zCand" << endl;
    h_deltaR_HadVelec1->Fill(reco::deltaR(vCand_, VZe1), weight_);
    h_deltaR_HadVelec2->Fill(reco::deltaR(vCand_, VZe2), weight_);
  }else if(zCand_.flavor() == PDG_ID_MUON){
    const TeVMuon & VZm1 = WPrimeUtil::Find(*zCand_.daughter(0), allMuons_);
    const TeVMuon & VZm2 = WPrimeUtil::Find(*zCand_.daughter(1), allMuons_);
    //cout << "Muon from loose zCand" << endl;
    h_deltaR_HadVmuon1->Fill(reco::deltaR(vCand_, VZm1), weight_);
    h_deltaR_HadVmuon2->Fill(reco::deltaR(vCand_, VZm2), weight_);
  }
}

void TTbarAnalyzer::fillJetMergingHistos(){

  double jet1jet2mass=-99.0;


  if (looseJets_.size()>1){
    h_deltaR_jet1jet2->Fill(reco::deltaR(looseJets_.at(0), looseJets_.at(1)), weight_);
    
    //Region histos
    if (looseJets_.at(0).mass() < 40){
      h_deltaR_jet1jet2_R1_cut40->Fill(reco::deltaR(looseJets_.at(0), looseJets_.at(1)), weight_);
    }
    if (looseJets_.at(0).mass() > 40){
      h_deltaR_jet1jet2_R2_cut40->Fill(reco::deltaR(looseJets_.at(0), looseJets_.at(1)), weight_);  
    }

    if (looseJets_.at(0).mass() < 50){
      h_deltaR_jet1jet2_R1_cut50->Fill(reco::deltaR(looseJets_.at(0), looseJets_.at(1)), weight_);
    }
    if (looseJets_.at(0).mass() > 50){
      h_deltaR_jet1jet2_R2_cut50->Fill(reco::deltaR(looseJets_.at(0), looseJets_.at(1)), weight_);  
    }

    if (looseJets_.at(0).mass() < 60){
      h_deltaR_jet1jet2_R1_cut60->Fill(reco::deltaR(looseJets_.at(0), looseJets_.at(1)), weight_);
    }
    if (looseJets_.at(0).mass() > 60){
      h_deltaR_jet1jet2_R2_cut60->Fill(reco::deltaR(looseJets_.at(0), looseJets_.at(1)), weight_);  
    }

    if (looseJets_.at(0).mass() < 65){
      h_deltaR_jet1jet2_R1_cut65->Fill(reco::deltaR(looseJets_.at(0), looseJets_.at(1)), weight_);
    }
    if (looseJets_.at(0).mass() > 65){
      h_deltaR_jet1jet2_R2_cut65->Fill(reco::deltaR(looseJets_.at(0), looseJets_.at(1)), weight_);  
    }

    if (looseJets_.at(0).mass() < 70){
      h_deltaR_jet1jet2_R1_cut70->Fill(reco::deltaR(looseJets_.at(0), looseJets_.at(1)), weight_);
    }
    if (looseJets_.at(0).mass() > 70){
      h_deltaR_jet1jet2_R2_cut70->Fill(reco::deltaR(looseJets_.at(0), looseJets_.at(1)), weight_);  
    }

    if (looseJets_.at(0).mass() < 80){
      h_deltaR_jet1jet2_R1_cut80->Fill(reco::deltaR(looseJets_.at(0), looseJets_.at(1)), weight_);
    }
    if (looseJets_.at(0).mass() > 80){
      h_deltaR_jet1jet2_R2_cut80->Fill(reco::deltaR(looseJets_.at(0), looseJets_.at(1)), weight_);  
    }

    if (looseJets_.at(0).mass() < 90){
      h_deltaR_jet1jet2_R1_cut90->Fill(reco::deltaR(looseJets_.at(0), looseJets_.at(1)), weight_);
    }
    if (looseJets_.at(0).mass() > 90){
      h_deltaR_jet1jet2_R2_cut90->Fill(reco::deltaR(looseJets_.at(0), looseJets_.at(1)), weight_);  
    }
    

    if (looseJets_.at(0).mass() < 65){
      h_deltaR_jet2Z_R1->Fill(reco::deltaR(looseJets_.at(1), zCand_), weight_);
      if (looseJets_.size()>2)
	h_deltaR_jet3Z_R1->Fill(reco::deltaR(looseJets_.at(2), zCand_), weight_);
    }

    if (looseJets_.at(0).mass() > 65){
      h_deltaR_jet2Z_R2->Fill(reco::deltaR(looseJets_.at(1), zCand_), weight_);
      if (looseJets_.size()>2)
	h_deltaR_jet3Z_R2->Fill(reco::deltaR(looseJets_.at(2), zCand_), weight_);
    }


    reco::CompositeCandidate j1j2;
    j1j2.addDaughter(looseJets_.at(0));
    j1j2.addDaughter(looseJets_.at(1));
    AddFourMomenta addFM;
    addFM.set(j1j2);
    jet1jet2mass = j1j2.mass();
    
    if (looseJets_.at(0).mass()>110.0 || looseJets_.at(0).mass()<60)
      {
	h_jet1jet2_mass_Restricted->Fill(j1j2.mass(), weight_);
      }

    h_jet1mass_jet2mass->Fill(looseJets_.at(0).mass(), looseJets_.at(1).mass(), weight_);


  }//# jets > 1 loop


  if (looseJets_.at(0).mass() < 65){
    if (looseJets_.size()>0)
      h_deltaR_jet1Z_R1->Fill(reco::deltaR(looseJets_.at(0), zCand_), weight_);
  }

  if (looseJets_.at(0).mass() > 65){
    if (looseJets_.size()>0)
      h_deltaR_jet1Z_R2->Fill(reco::deltaR(looseJets_.at(0), zCand_), weight_);
  }


  //Out of jets>1 loop
  int k = 0; //if zero, only 1 jet to HadV candidate
  double jet1mass=-99.0;
  double deltaR_jet1jet2=-99.0;
  reco::CompositeCandidate hadronicVZ;
  reco::CompositeCandidate hadronicVZF;
 

  if (looseJets_.at(0).mass() > 60 && looseJets_.at(0).mass() < 110)
    {
      jet1mass = looseJets_.at(0).mass(); 
 
    }
  
  if (looseJets_.size()>1)
    {
      deltaR_jet1jet2=reco::deltaR(looseJets_.at(0), looseJets_.at(1));
      if (fabs(jet1mass - 85.0) > fabs(jet1jet2mass - 85.0))
	{
	  k = 1;
	}

    }
  
  if (jet1mass > 0 && k==0)
    {
      h_HadV_mass_Cory->Fill(jet1mass,weight_);
      
      hadronicVZ.addDaughter(looseJets_.at(0));
      hadronicVZ.addDaughter(zCand_);
      AddFourMomenta addFRM;
      addFRM.set(hadronicVZ);
      
      h_HadVZmass_Cory->Fill(hadronicVZ.mass(),weight_);
    }

  else
    if (jet1jet2mass > 0 && k==1)
      {
	h_HadV_mass_Cory->Fill(jet1jet2mass, weight_);

	reco::CompositeCandidate j1j2;
	j1j2.addDaughter(looseJets_.at(0));
	j1j2.addDaughter(looseJets_.at(1));
	AddFourMomenta addFM;
	addFM.set(j1j2);
	
	hadronicVZ.addDaughter(j1j2);
	hadronicVZ.addDaughter(zCand_);
	AddFourMomenta addFRM;
	addFRM.set(hadronicVZ);

	h_HadVZmass_Cory->Fill(hadronicVZ.mass(),weight_);
      }

 
  if (jet1mass > 0 && k==0)
    {
      h_HadV_mass_Flavia->Fill(jet1mass, weight_);

      hadronicVZF.addDaughter(looseJets_.at(0));
      hadronicVZF.addDaughter(zCand_);
      AddFourMomenta addFRM;
      addFRM.set(hadronicVZF);

      h_HadVZmass_Flavia->Fill(hadronicVZF.mass(), weight_);
    }
  else
    if (jet1jet2mass > 0 && k==1 && deltaR_jet1jet2<2)
      {
  	h_HadV_mass_Flavia->Fill(jet1jet2mass, weight_);

	reco::CompositeCandidate j1j2;
	j1j2.addDaughter(looseJets_.at(0));
	j1j2.addDaughter(looseJets_.at(1));
	AddFourMomenta addFM;
	addFM.set(j1j2);
	
	hadronicVZF.addDaughter(j1j2);
	hadronicVZF.addDaughter(zCand_);
	AddFourMomenta addFRM;
	addFRM.set(hadronicVZF);

	h_HadVZmass_Flavia->Fill(hadronicVZF.mass(),weight_);
      }

}



void TTbarAnalyzer::fillPOGMuonHists(){
  for (size_t nM=0; nM<looseMuons_.size(); nM++)
    {
      //Muon POG work
      double trackpT = looseMuons_.at(nM).track()->pt();
      double trackpT2 = trackpT*trackpT;
      double invpt = 1/trackpT;
      double trackpTError = looseMuons_.at(nM).track()->ptError();
      double dptpt = trackpTError/trackpT;
      double dptpt2 =trackpTError/trackpT2;

      for (size_t i=0; i<looseMuons_.size(); i++)
	{
	  for (size_t j=0; j<genMuons.size(); j++)
	    {
	      if (deltaR(looseMuons_.at(i),genMuons.at(j)) < 0.01)
		{

		  double genpt = genMuons.at(j).pt();
		  double invgenpt = 1/genpt;
		  h_dptpt_vs_genpt->Fill(dptpt, genpt);
		  h_dptpt2_vs_genpt->Fill(dptpt2, genpt);

		  h_dptpt_vs_invgenpt->Fill(dptpt, invgenpt);
		  h_dptpt2_vs_invgenpt->Fill(dptpt2, invgenpt);

		}
	    }
	}

      double slope = 0.0115695;
      //Slope correction
      dptpt2 = dptpt2-(slope*invpt);

      h_dptpt2->Fill(dptpt2);
      h_dptpt_vs_pt->Fill(dptpt, trackpT);
      h_dptpt2_vs_pt->Fill(dptpt2, trackpT);

      h_dptpt_vs_invpt->Fill(dptpt, invpt);
      h_dptpt2_vs_invpt->Fill(dptpt2, invpt);


    }
}
