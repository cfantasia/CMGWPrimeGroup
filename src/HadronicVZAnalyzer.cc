#include "UserCode/CMGWPrimeGroup/interface/HadronicVZAnalyzer.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember)) 

using namespace std;
HadronicVZAnalyzer::HadronicVZAnalyzer(){}
HadronicVZAnalyzer::HadronicVZAnalyzer(const edm::ParameterSet & cfg, WPrimeUtil * wprimeUtil){

  wprimeUtil_ = wprimeUtil;
  assert(wprimeUtil_);

   Cuts_          = cfg.getParameter<vstring>("Cuts");
  NCuts_          = Cuts_.size();
  //FillCutFns();

  SetCandEvtFile(cfg.getParameter<string>("candEvtFile"));
  results_.assign(NCuts_,wprime::FilterEff());

  doPreselect_ = cfg.getParameter<bool>("preselect");

  debugme = cfg.getParameter<bool>("debugme");

  muonsLabel_ = cfg.getParameter<string>("muons");
  jetsLabel_ = cfg.getParameter<string>("jets");

  muonAlgo_ = cfg.getParameter<uint>("muonAlgo");
  if(debugme) cout<<"Using muon algo "<<muonAlgo_<<endl;

  
  hltEventLabel_ = cfg.getParameter<string>("hltEventTag");
  pileupLabel_ = cfg.getParameter<string>("pileupTag");

  triggersToUse_          = cfg.getParameter<vstring>("triggersToUse");


// +++++++++++++++++++Event characteristics
  

// +++++++++++++++++++General Cut values
  maxNumZs = cfg.getParameter<uint>("maxNumZs");
  minNLeptons = cfg.getParameter<uint>("minNLeptons");
  minNJets = cfg.getParameter<uint>("minNJets");
  maxNJets = cfg.getParameter<uint>("maxNJets");
  minLeadPt = cfg.getParameter<double>("minLeadPt");
  maxAngleBetweenJets = cfg.getParameter<double>("maxAngleBetweenJets");

// +++++++++++++++++++Z Cuts
  minZpt = cfg.getParameter<double>("minZpt");
  minZmass = cfg.getParameter<double>("minZmass");
  maxZmass = cfg.getParameter<double>("maxZmass");


// +++++++++++++++++++Hadronic Boson Cuts
  minHadVpt = cfg.getParameter<double>("minHadVPt");
  minHadVmass = cfg.getParameter<double>("minHadVmass");
  maxHadVmass = cfg.getParameter<double>("maxHadVmass");

// +++++++++++++++++++Muon General Cuts
  maxMuonEta = cfg.getParameter<double>("maxMuonEta");
  minMuonLoosePt = cfg.getParameter<double>("minMuonLoosePt");
  minMuonTightPt = cfg.getParameter<double>("minMuonTightPt");
//VBTF Recommended Cuts
  maxMuonDxy = cfg.getParameter<double>("maxMuonDxy");
  maxMuonNormChi2 = cfg.getParameter<double>("maxMuonNormChi2");
  minMuonNPixHit = cfg.getParameter<uint>("minMuonNPixHit");
  minMuonNTrkHit = cfg.getParameter<uint>("minMuonNTrkHit");
  minMuonStations = cfg.getParameter<uint>("minMuonStations");
  minMuonHitsUsed = cfg.getParameter<uint>("minMuonHitsUsed");

// +++++++++++++++++++Jet General Cuts
  minJetPt = cfg.getParameter<double>("minJetPt");
  maxJetEta = cfg.getParameter<double>("maxJetEta");
  maxJetNHF = cfg.getParameter<double>("maxJetNHF");
  maxJetNEF = cfg.getParameter<double>("maxJetNEF");
  minJetCHF = cfg.getParameter<double>("minJetCHF");
  maxJetCEF = cfg.getParameter<double>("maxJetCEF");
  minJetnumConst = cfg.getParameter<uint>("minJetnumConst");
  minJetcMult = cfg.getParameter<uint>("minJetcMult");
}

HadronicVZAnalyzer::~HadronicVZAnalyzer(){
  outCandEvt.close(); 
}


/// Declare Histograms
//--------------------------------------------------------------
void HadronicVZAnalyzer::Declare_Histos(TFileDirectory & dir)
{
  // Extend later.
  printf("Declare histos\n");
  //Loose histos
  h_HadVZMass = dir.make<TH1F>("h_HadVZMass","h_HadVZMass",250,0.0,2500.0);
  h_Zmuon1_pt = dir.make<TH1F>("h_Zmuon1_pt", "h_Zmuon1_pt", 100, 0.0, 1000.0);
  h_Zmuon1_eta = dir.make<TH1F>("h_Zmuon1_eta", "h_Zmuon1_eta", 40, -5.0, 5.0);
  h_Zmuon1_phi = dir.make<TH1F>("h_Zmuon1_phi", "h_Zmuon1_phi", 20, -4.0, 4.0);
  h_Zmuon2_pt = dir.make<TH1F>("h_Zmuon2_pt", "h_Zmuon2_pt", 100, 0.0, 1000.0);
  h_Zmuon2_eta = dir.make<TH1F>("h_Zmuon2_eta", "h_Zmuon2_eta", 40, -5.0, 5.0);
  h_Zmuon2_phi = dir.make<TH1F>("h_Zmuon2_phi", "h_Zmuon2_phi", 20, -4.0, 4.0);
  h_Zmuon1_VZCut_pt = dir.make<TH1F>("h_Zmuon1_VZCut_pt", "h_Zmuon1_VZCut_pt", 100, 0.0, 1000.0);
  h_Zmuon1_VZCut_eta = dir.make<TH1F>("h_Zmuon1_VZCut_eta", "h_Zmuon1_VZCut_eta", 40, -5.0, 5.0);
  h_Zmuon1_VZCut_phi = dir.make<TH1F>("h_Zmuon1_VZCut_phi", "h_Zmuon1_VZCut_phi", 20, -4.0, 4.0);
  h_Zmuon2_VZCut_pt = dir.make<TH1F>("h_Zmuon2_VZCut_pt", "h_Zmuon2_VZCut_pt", 100, 0.0, 1000.0);
  h_Zmuon2_VZCut_eta = dir.make<TH1F>("h_Zmuon2_VZCut_eta", "h_Zmuon2_VZCut_eta", 40, -5.0, 5.0);
  h_Zmuon2_VZCut_phi = dir.make<TH1F>("h_Zmuon2_VZCut_phi", "h_Zmuon2_VZCut_phi", 20, -4.0, 4.0);
  h_deltaR_muon1muon2 = dir.make<TH1F>("h_deltaR_muon1muon2", "h_deltaR_muon1muon2", 50, 0., 5.);
  h_deltaR_HadVmuon1 = dir.make<TH1F>("h_deltaR_HadVmuon1", "h_deltaR_HadVmuon1", 50, 0., 5.);
  h_deltaR_HadVmuon2 = dir.make<TH1F>("h_deltaR_HadVmuon2", "h_deltaR_HadVmuon2", 50, 0., 5.);
  h_jet1_pt = dir.make<TH1F>("h_jet1_pt", "h_jet1_pt", 100, 0.0, 1000.0);
  h_jet1_eta = dir.make<TH1F>("h_jet1_eta", "h_jet1_eta", 40, -5.0, 5.0);
  h_jet1_phi = dir.make<TH1F>("h_jet1_phi", "h_jet1_phi", 20, -4.0, 4.0);
  h_jet2_pt = dir.make<TH1F>("h_jet2_pt", "h_jet2_pt", 100, 0.0, 1000.0);
  h_jet2_eta = dir.make<TH1F>("h_jet2_eta", "h_jet2_eta", 40, -5.0, 5.0);
  h_jet2_phi = dir.make<TH1F>("h_jet2_phi", "h_jet2_phi", 20, -4.0, 4.0);
  h_jet1_mass = dir.make<TH1F>("h_jet1_mass", "h_jet1_mass", 100, 0.0, 500.);
  h_jet2_mass = dir.make<TH1F>("h_jet2_mass", "h_jet2_mass", 100, 0.0, 500.);
  h_deltaR_jet1muon1 = dir.make<TH1F>("h_deltaR_jet1muon1", "h_deltaR_jet1muon1", 50, 0., 5.);
  h_deltaR_jet1muon2 = dir.make<TH1F>("h_deltaR_jet1muon2", "h_deltaR_jet1muon2", 50, 0., 5.);
  h_deltaR_jet2muon1 = dir.make<TH1F>("h_deltaR_jet2muon1", "h_deltaR_jet2muon1", 50, 0., 5.);
  h_deltaR_jet2muon2 = dir.make<TH1F>("h_deltaR_jet2muon2", "h_deltaR_jet2muon2", 50, 0., 5.);
  h_HadVZpt = dir.make<TH1F>("h_HadVZpt", "h_HadVZpt", 30, 0.0, 300.0);
  h_HadVZeta = dir.make<TH1F>("h_HadVZeta", "h_HadVZeta", 40, -5.0, 5.0);
  h_HadVZphi = dir.make<TH1F>("h_HadVZphi", "h_HadVZphi", 20, -4.0, 4.0);
  h_muons_pt = dir.make<TH1F>("h_muons_pt", "h_muons_pt", 100, 0.0, 1000.);
  h_muons_eta = dir.make<TH1F>("h_muons_eta", "h_muons_eta", 40, -5.0, 5.);  
  h_muons_phi = dir.make<TH1F>("h_muons_phi", "h_muons_phi", 20, -4.0, 4.);
  h_jets_pt = dir.make<TH1F>("h_jets_pt", "h_jets_pt", 100, 0.0, 1000.);
  h_jets_eta = dir.make<TH1F>("h_jets_eta", "h_jets_eta", 40, -5.0, 5.);  
  h_jets_phi = dir.make<TH1F>("h_jets_phi", "h_jets_phi", 20, -4.0, 4.);
  h_jet_HadV_pt = dir.make<TH1F>("h_jet_HadV_pt", "h_jet_HadV_pt", 100, 0.0, 1000.);
  h_jet_HadV_eta = dir.make<TH1F>("h_jet_HadV_eta", "h_jet_HadV_eta", 40, -5.0, 5.);  
  h_jet_HadV_phi = dir.make<TH1F>("h_jet_HadV_phi", "h_jet_HadV_phi", 20, -4.0, 4.);
  h_jet_VZCut_pt = dir.make<TH1F>("h_jet_VZCut_pt", "h_jet_VZCut_pt", 100, 0.0, 1000.);
  h_jet_VZCut_eta = dir.make<TH1F>("h_jet_VZCut_eta", "h_jet_VZCut_eta", 40, -5.0, 5.);  
  h_jet_VZCut_phi = dir.make<TH1F>("h_jet_VZCut_phi", "h_jet_VZCut_phi", 20, -4.0, 4.);
  h_jet_mult = dir.make<TH1F>("h_jet_mult", "h_jet_mult", 10, -0.5, 9.5); 
  h_jet_mult_inc = dir.make<TH1F>("h_jet_mult_incl", "h_jet_mult_incl", 10, -0.5, 9.5); 
  h_jet1jet2_mass = dir.make<TH1F>("h_jet1jet2_mass", "h_jet1jet2_mass", 100, 0.0, 500.);
  //tight histos
  h_tight_HadVZMass = dir.make<TH1F>("h_tight_HadVZMass","h_tight_HadVZMass",250,0.0,2500.0);
  h_tight_Zmuon1_pt = dir.make<TH1F>("h_tight_Zmuon1_pt", "h_tight_Zmuon1_pt", 100, 0.0, 1000.0);
  h_tight_Zmuon1_eta = dir.make<TH1F>("h_tight_Zmuon1_eta", "h_tight_Zmuon1_eta", 40, -5.0, 5.0);
  h_tight_Zmuon1_phi = dir.make<TH1F>("h_tight_Zmuon1_phi", "h_tight_Zmuon1_phi", 20, -4.0, 4.0);
  h_tight_Zmuon2_pt = dir.make<TH1F>("h_tight_Zmuon2_pt", "h_tight_Zmuon2_pt", 100, 0.0, 1000.0);
  h_tight_Zmuon2_eta = dir.make<TH1F>("h_tight_Zmuon2_eta", "h_tight_Zmuon2_eta", 40, -5.0, 5.0);
  h_tight_Zmuon2_phi = dir.make<TH1F>("h_tight_Zmuon2_phi", "h_tight_Zmuon2_phi", 20, -4.0, 4.0);
  h_tight_Zmuon1_VZCut_pt = dir.make<TH1F>("h_tight_Zmuon1_VZCut_pt", "h_tight_Zmuon1_VZCut_pt", 100, 0.0, 1000.0);
  h_tight_Zmuon1_VZCut_eta = dir.make<TH1F>("h_tight_Zmuon1_VZCut_eta", "h_tight_Zmuon1_VZCut_eta", 40, -5.0, 5.0);
  h_tight_Zmuon1_VZCut_phi = dir.make<TH1F>("h_tight_Zmuon1_VZCut_phi", "h_tight_Zmuon1_VZCut_phi", 20, -4.0, 4.0);
  h_tight_Zmuon2_VZCut_pt = dir.make<TH1F>("h_tight_Zmuon2_VZCut_pt", "h_tight_Zmuon2_VZCut_pt", 100, 0.0, 1000.0);
  h_tight_Zmuon2_VZCut_eta = dir.make<TH1F>("h_tight_Zmuon2_VZCut_eta", "h_tight_Zmuon2_VZCut_eta", 40, -5.0, 5.0);
  h_tight_Zmuon2_VZCut_phi = dir.make<TH1F>("h_tight_Zmuon2_VZCut_phi", "h_tight_Zmuon2_VZCut_phi", 20, -4.0, 4.0);
  h_tight_deltaR_muon1muon2 = dir.make<TH1F>("h_tight_deltaR_muon1muon2", "h_tight_deltaR_muon1muon2", 50, 0., 5.);
  h_tight_deltaR_HadVmuon1 = dir.make<TH1F>("h_tight_deltaR_HadVmuon1", "h_tight_deltaR_HadVmuon1", 50, 0., 5.);
  h_tight_deltaR_HadVmuon2 = dir.make<TH1F>("h_tight_deltaR_HadVmuon2", "h_tight_deltaR_HadVmuon2", 50, 0., 5.);
  /*  h_tight_jet1_pt = dir.make<TH1F>("h_tight_jet1_pt", "h_tight_jet1_pt", 100, 0.0, 1000.0);
      h_tight_jet1_eta = dir.make<TH1F>("h_tight_jet1_eta", "h_tight_jet1_eta", 40, -5.0, 5.0);
      h_tight_jet1_phi = dir.make<TH1F>("h_tight_jet1_phi", "h_tight_jet1_phi", 20, -4.0, 4.0);
      h_tight_jet2_pt = dir.make<TH1F>("h_tight_jet2_pt", "h_tight_jet2_pt", 100, 0.0, 1000.0);
      h_tight_jet2_eta = dir.make<TH1F>("h_tight_jet2_eta", "h_tight_jet2_eta", 40, -5.0, 5.0);
      h_tight_jet2_phi = dir.make<TH1F>("h_tight_jet2_phi", "h_tight_jet2_phi", 20, -4.0, 4.0);
      h_tight_jet1_mass = dir.make<TH1F>("h_tight_jet1_mass", "h_tight_jet1_mass", 100, 0.0, 500.);
      h_tight_jet2_mass = dir.make<TH1F>("h_tight_jet2_mass", "h_tight_jet2_mass", 100, 0.0, 500.);
      h_tight_deltaR_jet1muon1 = dir.make<TH1F>("h_tight_deltaR_jet1muon1", "h_tight_deltaR_jet1muon1", 50, 0., 5.);
      h_tight_deltaR_jet1muon2 = dir.make<TH1F>("h_tight_deltaR_jet1muon2", "h_tight_deltaR_jet1muon2", 50, 0., 5.);
      h_tight_deltaR_jet2muon1 = dir.make<TH1F>("h_tight_deltaR_jet2muon1", "h_tight_deltaR_jet2muon1", 50, 0., 5.);
      h_tight_deltaR_jet2muon2 = dir.make<TH1F>("h_tight_deltaR_jet2muon2", "h_tight_deltaR_jet2muon2", 50, 0., 5.);*/
  h_tight_HadVZpt = dir.make<TH1F>("h_tight_HadVZpt", "h_tight_HadVZpt", 30, 0.0, 300.0);
  h_tight_HadVZeta = dir.make<TH1F>("h_tight_HadVZeta", "h_tight_HadVZeta", 40, -5.0, 5.0);
  h_tight_HadVZphi = dir.make<TH1F>("h_tight_HadVZphi", "h_tight_HadVZphi", 20, -4.0, 4.0);
  h_tight_muons_pt = dir.make<TH1F>("h_tight_muons_pt", "h_tight_muons_pt", 100, 0.0, 1000.);
  h_tight_muons_eta = dir.make<TH1F>("h_tight_muons_eta", "h_tight_muons_eta", 40, -5.0, 5.);  
  h_tight_muons_phi = dir.make<TH1F>("h_tight_muons_phi", "h_tight_muons_phi", 20, -4.0, 4.);
  /*  h_tight_jets_pt = dir.make<TH1F>("h_tight_jets_pt", "h_tight_jets_pt", 100, 0.0, 1000.);
      h_tight_jets_eta = dir.make<TH1F>("h_tight_jets_eta", "h_tight_jets_eta", 40, -5.0, 5.);  
      h_tight_jets_phi = dir.make<TH1F>("h_tight_jets_phi", "h_tight_jets_phi", 20, -4.0, 4.);
      h_tight_jet_HadV_pt = dir.make<TH1F>("h_tight_jet_HadV_pt", "h_tight_jet_HadV_pt", 100, 0.0, 1000.);
      h_tight_jet_HadV_eta = dir.make<TH1F>("h_tight_jet_HadV_eta", "h_tight_jet_HadV_eta", 40, -5.0, 5.);  
      h_tight_jet_HadV_phi = dir.make<TH1F>("h_tight_jet_HadV_phi", "h_tight_jet_HadV_phi", 20, -4.0, 4.);*/
  h_tight_jet_VZCut_pt = dir.make<TH1F>("h_tight_jet_VZCut_pt", "h_tight_jet_VZCut_pt", 100, 0.0, 1000.);
  h_tight_jet_VZCut_eta = dir.make<TH1F>("h_tight_jet_VZCut_eta", "h_tight_jet_VZCut_eta", 40, -5.0, 5.);  
  h_tight_jet_VZCut_phi = dir.make<TH1F>("h_tight_jet_VZCut_phi", "h_tight_jet_VZCut_phi", 20, -4.0, 4.);
  //  h_tight_jet_mult = dir.make<TH1F>("h_tight_jet_mult", "h_tight_jet_mult", 10, -0.5, 9.5); 
  // h_tight_jet_mult_inc = dir.make<TH1F>("h_tight_jet_mult_incl", "h_tight_jet_mult_incl", 10, -0.5, 9.5); 
  // h_tight_jet1jet2_mass = dir.make<TH1F>("h_tight_jet1jet2_mass", "h_tight_jet1jet2_mass", 100, 0.0, 500.);

  // Muon POG histos
  h_dptpt2 = dir.make<TH1F>("h_dptpt2", "h_dptpt2", 500, 0., 0.1);
  h_dptpt_vs_pt = dir.make<TH2F>("h_dptpt_vs_pt", "h_dptpt_vs_pt", 50, 0., 0.5, 50, 0.0, 500.);
  h_dptpt2_vs_pt = dir.make<TH2F>("h_dptpt2_vs_pt", "h_dptpt2_vs_pt", 500, 0., 0.1, 150, 0.0, 500.);
  h_dptpt_vs_invpt = dir.make<TH2F>("h_dptpt_vs_invpt","h_dptpt_vs_invpt", 50, 0., 0.5, 50, 0., 0.1);
  h_dptpt2_vs_invpt = dir.make<TH2F>("h_dptpt2_vs_invpt","h_dptpt2_vs_invpt", 500, 0., 0.1, 150, 0., 0.1);
  h_dptpt_vs_eta = dir.make<TH2F>("h_dptpt_vs_eta","h_dptpt_vs_eta", 50, 0., 0.5, 50, -5., 5.);
  h_dptpt2_vs_eta = dir.make<TH2F>("h_dptpt2_vs_eta","h_dptpt2_vs_eta", 500, 0., 0.1, 150, -5., 5.);
  h_dptpt_vs_pt_meta09 = dir.make<TH2F>("h_dptpt_vs_pt_meta09", "h_dptpt_vs_pt_meta09", 50, 0., 0.5, 50, 0.0, 500.);
  h_dptpt2_vs_pt_meta09 = dir.make<TH2F>("h_dptpt2_vs_pt_meta09", "h_dptpt2_vs_pt_meta09", 500, 0., 0.1, 150, 0.0, 500.);
  h_dptpt_vs_invpt_meta09 = dir.make<TH2F>("h_dptpt_vs_invpt_meta09","h_dptpt_vs_invpt_meta09", 50, 0., 0.5, 50, 0., 0.1);
  h_dptpt2_vs_invpt_meta09 = dir.make<TH2F>("h_dptpt2_vs_invpt_meta09","h_dptpt2_vs_invpt_meta09", 500, 0., 0.1, 150, 0., 0.1);
  h_dptpt_vs_pt_meta0912 = dir.make<TH2F>("h_dptpt_vs_pt_meta0912", "h_dptpt_vs_pt_meta0912", 50, 0., 0.5, 50, 0.0, 500.);
  h_dptpt2_vs_pt_meta0912 = dir.make<TH2F>("h_dptpt2_vs_pt_meta0912", "h_dptpt2_vs_pt_meta0912", 500, 0., 0.1, 150, 0.0, 500.);
  h_dptpt_vs_invpt_meta0912 = dir.make<TH2F>("h_dptpt_vs_invpt_meta0912","h_dptpt_vs_invpt_meta0912", 50, 0., 0.5, 50, 0., 0.1);
  h_dptpt2_vs_invpt_meta0912 = dir.make<TH2F>("h_dptpt2_vs_invpt_meta0912","h_dptpt2_vs_invpt_meta0912", 500, 0., 0.1, 150, 0., 0.1);
  h_dptpt_vs_pt_meta1225 = dir.make<TH2F>("h_dptpt_vs_pt_meta1225", "h_dptpt_vs_pt_meta1225", 50, 0., 0.5, 50, 0.0, 500.);
  h_dptpt2_vs_pt_meta1225 = dir.make<TH2F>("h_dptpt2_vs_pt_meta1225", "h_dptpt2_vs_pt_meta1225", 500, 0., 0.1, 150, 0.0, 500.);
  h_dptpt_vs_invpt_meta1225 = dir.make<TH2F>("h_dptpt_vs_invpt_meta1225","h_dptpt_vs_invpt_meta1225", 50, 0., 0.5, 50, 0., 0.1);
  h_dptpt2_vs_invpt_meta1225 = dir.make<TH2F>("h_dptpt2_vs_invpt_meta1225","h_dptpt2_vs_invpt_meta1225", 500, 0., 0.1, 150, 0., 0.1);
  h_dptpt_vs_pt_teta09 = dir.make<TH2F>("h_dptpt_vs_pt_teta09", "h_dptpt_vs_pt_teta09", 50, 0., 0.5, 50, 0.0, 500.);
  h_dptpt2_vs_pt_teta09 = dir.make<TH2F>("h_dptpt2_vs_pt_teta09", "h_dptpt2_vs_pt_teta09", 500, 0., 0.1, 150, 0.0, 500.);
  h_dptpt_vs_invpt_teta09 = dir.make<TH2F>("h_dptpt_vs_invpt_teta09","h_dptpt_vs_invpt_teta09", 50, 0., 0.5, 50, 0., 0.1);
  h_dptpt2_vs_invpt_teta09 = dir.make<TH2F>("h_dptpt2_vs_invpt_teta09","h_dptpt2_vs_invpt_teta09", 500, 0., 0.1, 150, 0., 0.1);
  h_dptpt_vs_pt_teta0915 = dir.make<TH2F>("h_dptpt_vs_pt_teta0915", "h_dptpt_vs_pt_teta0915", 50, 0., 0.5, 50, 0.0, 500.);
  h_dptpt2_vs_pt_teta0915 = dir.make<TH2F>("h_dptpt2_vs_pt_teta0915", "h_dptpt2_vs_pt_teta0915", 500, 0., 0.1, 150, 0.0, 500.);
  h_dptpt_vs_invpt_teta0915 = dir.make<TH2F>("h_dptpt_vs_invpt_teta0915","h_dptpt_vs_invpt_teta0915", 50, 0., 0.5, 50, 0., 0.1);
  h_dptpt2_vs_invpt_teta0915 = dir.make<TH2F>("h_dptpt2_vs_invpt_teta0915","h_dptpt2_vs_invpt_teta0915", 500, 0., 0.1, 150, 0., 0.1);
  h_dptpt_vs_pt_teta1524 = dir.make<TH2F>("h_dptpt_vs_pt_teta1524", "h_dptpt_vs_pt_teta1524", 50, 0., 0.5, 50, 0.0, 500.);
  h_dptpt2_vs_pt_teta1524 = dir.make<TH2F>("h_dptpt2_vs_pt_teta1524", "h_dptpt2_vs_pt_teta1524", 500, 0., 0.1, 150, 0.0, 500.);
  h_dptpt_vs_invpt_teta1524 = dir.make<TH2F>("h_dptpt_vs_invpt_teta1524","h_dptpt_vs_invpt_teta1524", 50, 0., 0.5, 50, 0., 0.1);
  h_dptpt2_vs_invpt_teta1524 = dir.make<TH2F>("h_dptpt2_vs_invpt_teta1524","h_dptpt2_vs_invpt_teta1524", 500, 0., 0.1, 150, 0., 0.1);
  //

  cout << "Histos declared" << endl;

}//Declare_Histos

//Fill Histograms
//-----------------------------------------------------------
void HadronicVZAnalyzer::Fill_Histos(int index, float weight)
{
  if(debugme) printf("Filling Histos\n");
}//Fill_Histos

void HadronicVZAnalyzer::FillJetMultiplicityHists(){
  h_jet_mult->Fill(jets_.size(), weight_);

  double zero = 0.0;
  h_jet_mult_inc->Fill(zero, weight_);

  if (jets_.size()>=1)
  {
    h_jet_mult_inc->Fill(1, weight_);
  }
  if (jets_.size()>=2)
  {
    h_jet_mult_inc->Fill(2, weight_);
  }
  if (jets_.size()>=3)
  {
    h_jet_mult_inc->Fill(3, weight_);
  }
  if (jets_.size()>=4)
  {
    h_jet_mult_inc->Fill(4, weight_);
  }
  if (jets_.size()>=5)
  {
    h_jet_mult_inc->Fill(5, weight_);
  }
  if (jets_.size()>=6)
  {
    h_jet_mult_inc->Fill(6, weight_);
  }
  if (jets_.size()>=7)
  {
    h_jet_mult_inc->Fill(7, weight_);
  }
  if (jets_.size()>=8)
  {
    h_jet_mult_inc->Fill(8, weight_);
  }
  if (jets_.size()>=9)
  {
    h_jet_mult_inc->Fill(9, weight_);
  }
  if (jets_.size()>=10)
  {
    h_jet_mult_inc->Fill(10, weight_);
  }

}

void HadronicVZAnalyzer::FillLooseMuonHists(){
  for (size_t nM=0; nM<looseMuons_.size(); nM++)
  {
    h_muons_pt->Fill(looseMuons_.at(nM).pt(), weight_);
    h_muons_eta->Fill(looseMuons_.at(nM).eta(), weight_);
    h_muons_phi->Fill(looseMuons_.at(nM).phi(), weight_);

    //Muon POG work
    double trackpT = looseMuons_.at(nM).track()->pt();
    double trackpT2 = trackpT*trackpT;
    double invpt = 1/trackpT;
    double trackpTError = looseMuons_.at(nM).track()->ptError();
    double dptpt = trackpTError/trackpT;
    double dptpt2 =trackpTError/trackpT2;

    h_dptpt2->Fill(dptpt2);
    h_dptpt_vs_pt->Fill(dptpt, trackpT);
    h_dptpt2_vs_pt->Fill(dptpt2, trackpT);
    h_dptpt_vs_invpt->Fill(dptpt, invpt);
    h_dptpt2_vs_invpt->Fill(dptpt2, invpt);
    h_dptpt_vs_eta->Fill(dptpt, looseMuons_.at(nM).track()->eta());
    h_dptpt2_vs_eta->Fill(dptpt2, looseMuons_.at(nM).track()->eta());

    if (fabs(looseMuons_.at(nM).track()->eta()) < 0.9)
    {
      h_dptpt_vs_pt_meta09->Fill(dptpt, trackpT);
      h_dptpt2_vs_pt_meta09->Fill(dptpt2, trackpT);
      h_dptpt_vs_invpt_meta09->Fill(dptpt, invpt);
      h_dptpt2_vs_invpt_meta09->Fill(dptpt2, invpt);
    }
    if (fabs(looseMuons_.at(nM).track()->eta()) >= 0.9 && fabs(looseMuons_.at(nM).track()->eta()) < 1.2)
    {
      h_dptpt_vs_pt_meta0912->Fill(dptpt, trackpT);
      h_dptpt2_vs_pt_meta0912->Fill(dptpt2, trackpT);
      h_dptpt_vs_invpt_meta0912->Fill(dptpt, invpt);
      h_dptpt2_vs_invpt_meta0912->Fill(dptpt2, invpt);
    }
    if (fabs(looseMuons_.at(nM).track()->eta()) >= 1.2 && fabs(looseMuons_.at(nM).track()->eta()) < 2.5)
    {
      h_dptpt_vs_pt_meta1225->Fill(dptpt, trackpT);
      h_dptpt2_vs_pt_meta1225->Fill(dptpt2, trackpT);
      h_dptpt_vs_invpt_meta1225->Fill(dptpt, invpt);
      h_dptpt2_vs_invpt_meta1225->Fill(dptpt2, invpt);
    }
    if (fabs(looseMuons_.at(nM).eta()) < 0.9)
    {
      h_dptpt_vs_pt_teta09->Fill(dptpt, trackpT);
      h_dptpt2_vs_pt_teta09->Fill(dptpt2, trackpT);
      h_dptpt_vs_invpt_teta09->Fill(dptpt, invpt);
      h_dptpt2_vs_invpt_teta09->Fill(dptpt2, invpt);
    }
    if (fabs(looseMuons_.at(nM).eta()) >= 0.9 && fabs(looseMuons_.at(nM).eta()) < 1.5)
    {
      h_dptpt_vs_pt_teta0915->Fill(dptpt, trackpT);
      h_dptpt2_vs_pt_teta0915->Fill(dptpt2, trackpT);
      h_dptpt_vs_invpt_teta0915->Fill(dptpt, invpt);
      h_dptpt2_vs_invpt_teta0915->Fill(dptpt2, invpt);
    }
    if (fabs(looseMuons_.at(nM).eta()) >= 1.5 && fabs(looseMuons_.at(nM).eta()) < 2.4)
    {
      h_dptpt_vs_pt_teta1524->Fill(dptpt, trackpT);
      h_dptpt2_vs_pt_teta1524->Fill(dptpt2, trackpT);
      h_dptpt_vs_invpt_teta1524->Fill(dptpt, invpt);
      h_dptpt2_vs_invpt_teta1524->Fill(dptpt2, invpt);
    }
  }
}

void HadronicVZAnalyzer::FillTightMuonHists(){
  for (size_t nM=0; nM<tightMuons_.size(); nM++)
  {
    h_tight_muons_pt->Fill(tightMuons_.at(nM).pt(), weight_);
    h_tight_muons_eta->Fill(tightMuons_.at(nM).eta(), weight_);
    h_tight_muons_phi->Fill(tightMuons_.at(nM).phi(), weight_);
  }
}

void HadronicVZAnalyzer::FillGoodZHistos(){
  const TeVMuon & m1 = Find(*zCand_.daughter(0), muons_);
  const TeVMuon & m2 = Find(*zCand_.daughter(1), muons_);
  if (debugme)
    cout << "Found my muons from loose Z" << endl;
  h_Zmuon1_pt->Fill(m1.pt(), weight_);
  h_Zmuon1_eta->Fill(m1.eta(), weight_);
  h_Zmuon1_phi->Fill(m1.phi(), weight_);
  h_Zmuon2_pt->Fill(m2.pt(), weight_);
  h_Zmuon2_eta->Fill(m2.eta(), weight_);
  h_Zmuon2_phi->Fill(m2.phi(), weight_);	  
  h_deltaR_muon1muon2->Fill(reco::deltaR(m1, m2), weight_);
}

void HadronicVZAnalyzer::FillGoodZTightHistos(){
  const TeVMuon & m1 = Find(*zCandTight_.daughter(0), muons_);
  const TeVMuon & m2 = Find(*zCandTight_.daughter(1), muons_);
  if (debugme)
    cout << "Found my muons from tight Z" << endl;
  h_tight_Zmuon1_pt->Fill(m1.pt(), weight_);
  h_tight_Zmuon1_eta->Fill(m1.eta(), weight_);
  h_tight_Zmuon1_phi->Fill(m1.phi(), weight_);
  h_tight_Zmuon2_pt->Fill(m2.pt(), weight_);
  h_tight_Zmuon2_eta->Fill(m2.eta(), weight_);
  h_tight_Zmuon2_phi->Fill(m2.phi(), weight_);	  
  h_tight_deltaR_muon1muon2->Fill(reco::deltaR(m1, m2), weight_);
}

void HadronicVZAnalyzer::FillGoodHadVHistos(){
  h_jet_HadV_pt->Fill(wCand_.pt(), weight_);
  h_jet_HadV_eta->Fill(wCand_.eta(), weight_);
  h_jet_HadV_phi->Fill(wCand_.phi(), weight_);
  if (debugme)
    cout << "Filled my HadV histos" << endl;
}

void HadronicVZAnalyzer::FillValidVZHistos(){
  h_HadVZMass->Fill(hadVZ_.mass(), weight_);
  h_HadVZpt->Fill(hadVZ_.pt(), weight_);
  h_HadVZeta->Fill(hadVZ_.eta(), weight_);
  h_HadVZphi->Fill(hadVZ_.phi(), weight_);
  if (debugme)
    cout << "Filled my histos from HadVZ" << endl;
  const TeVMuon & VZm1 = Find(*zCand_.daughter(0), muons_);
  const TeVMuon & VZm2 = Find(*zCand_.daughter(1), muons_);
  //cout << "Muon from loose zCand" << endl;
  h_Zmuon1_VZCut_pt->Fill(VZm1.pt(), weight_);
  h_Zmuon1_VZCut_eta->Fill(VZm1.eta(), weight_);
  h_Zmuon1_VZCut_phi->Fill(VZm1.phi(), weight_);
  h_Zmuon2_VZCut_pt->Fill(VZm2.pt(), weight_);
  h_Zmuon2_VZCut_eta->Fill(VZm2.eta(), weight_);
  h_Zmuon2_VZCut_phi->Fill(VZm2.phi(), weight_);
  
  h_jet_VZCut_pt->Fill(wCand_.pt(), weight_);
  h_jet_VZCut_eta->Fill(wCand_.eta(), weight_);
  h_jet_VZCut_phi->Fill(wCand_.phi(), weight_);
  
  h_deltaR_HadVmuon1->Fill(reco::deltaR(wCand_, VZm1), weight_);
  h_deltaR_HadVmuon2->Fill(reco::deltaR(wCand_, VZm2), weight_);
}

void HadronicVZAnalyzer::FillGoodVZHistos(){
  h_tight_HadVZMass->Fill(hadVZTight_.mass(), weight_);
  h_tight_HadVZpt->Fill(hadVZTight_.pt(), weight_);
  h_tight_HadVZeta->Fill(hadVZTight_.eta(), weight_);
  h_tight_HadVZphi->Fill(hadVZTight_.phi(), weight_);
  if (debugme)
    cout << "Filled my HadVZTight histos" << endl;
  const TeVMuon & VZm1 = Find(*zCandTight_.daughter(0), muons_);
  const TeVMuon & VZm2 = Find(*zCandTight_.daughter(1), muons_);
  if (debugme)
    cout << "Found my tight muons from HadVZ" << endl;
  h_tight_Zmuon1_VZCut_pt->Fill(VZm1.pt(), weight_);
  h_tight_Zmuon1_VZCut_eta->Fill(VZm1.eta(), weight_);
  h_tight_Zmuon1_VZCut_phi->Fill(VZm1.phi(), weight_);
  h_tight_Zmuon2_VZCut_pt->Fill(VZm2.pt(), weight_);
  h_tight_Zmuon2_VZCut_eta->Fill(VZm2.eta(), weight_);
  h_tight_Zmuon2_VZCut_phi->Fill(VZm2.phi(), weight_);

  h_tight_jet_VZCut_pt->Fill(wCand_.pt(), weight_);
  if (debugme)
    cout << "Filled more tight muons histos from HadVZ -- 1" << endl;
  h_tight_jet_VZCut_eta->Fill(wCand_.eta(), weight_);
  if (debugme)
    cout << "Filled more tight muons histos from HadVZ -- 2" << endl;
  h_tight_jet_VZCut_phi->Fill(wCand_.phi(), weight_);
  if (debugme)
    cout << "Filled more tight muons histos from HadVZ -- 3" << endl;
	
  h_tight_deltaR_HadVmuon1->Fill(reco::deltaR(wCand_, VZm1), weight_);
  if (debugme)
    cout << "Filled more tight muons histos from HadVZ -- 4" << endl;
  h_tight_deltaR_HadVmuon2->Fill(reco::deltaR(wCand_, VZm2), weight_);
  if (debugme)
    cout << "Filled more tight muons histos from HadVZ -- 5" << endl;
}

void 
HadronicVZAnalyzer::eventLoop(edm::EventBase const & event){
  ClearEvtVariables();

  // Preselection - skip events that don't look promising
  if (doPreselect_){
    if(debugme) cout<<"Testing Preselection...\n";
    // We could setup some preselection here. To be implemented.
  }

  // Get leptons
  //////////////////////
  ////Deal With Muons///
  //////////////////////
  const vector<pat::Muon    > patMuons     = getProduct<vector<pat::Muon    > >(event, muonsLabel_);
  if(patMuons.size() < 2){
    cout << "Not enough leptons. Bad bad event, returning now..." << endl;
    return;
  }
  if(debugme){
    printf("    Contains: %i pat muon(s)\n",
           (int)patMuons.size());
  }

  //PU STUFF
  float PU_Weight = 1.;
  if(!wprimeUtil_->runningOnData()){//Don't do this for data
    PupInfo_ = getProduct<std::vector< PileupSummaryInfo > >(event, pileupLabel_);   
    PU_Weight = wprimeUtil_->getPUWeight3BX(PupInfo_);
    if(debugme) 
      cout <<" PU Weight: "<<PU_Weight
           <<endl;   
  }//MC Only If
  
  weight_ = wprimeUtil_->getWeight()*PU_Weight;

  // Make vectors of leptons passing various criteria
  // Loop over muons, and see if they pass the TeVMuon criteria  
  for (size_t i = 0; i < patMuons.size(); i++) {
    muons_.push_back(TeVMuon(patMuons[i],muonAlgo_));   

    //  bool LooseMuonCuts = PassMuonLoosePtCut(&muons_[i]) && PassMuonGlobalCut(&muons_[i]) && PassMuonNpixhitCut(&muons_[i]) && PassMuonNtrkhitCut(&muons_[i]) && PassMuonHitsUsedCut(&muons_[i]) && PassMuonStationsCut(&muons_[i]) && PassMuonEtaCut(&muons_[i]) && PassMuonDxyCut(&muons_[i]);
    //    bool LooseMuonCuts = PassMuonLoosePtCut(&muons_[i]) && PassMuonGlobalCut(&muons_[i]) && PassMuonNpixhitCut(&muons_[i]) && PassMuonNtrkhitCut(&muons_[i]) && PassMuonHitsUsedCut(&muons_[i]) && PassMuonStationsCut(&muons_[i]) && PassMuonEtaCut(&muons_[i]) && PassMuonDxyCut(&muons_[i]);
    
    bool LooseMuonCuts = PassMuonLoosePtCut(&muons_[i]) && PassMuonGlobalCut(&muons_[i]) && PassMuonNpixhitCut(&muons_[i]) && PassMuonHitsUsedCut(&muons_[i]) && PassMuonStationsCut(&muons_[i]) && PassMuonEtaCut(&muons_[i]) && PassMuonDxyCut(&muons_[i]);
    if (LooseMuonCuts)
      looseMuons_.push_back(muons_[i]);
    
    bool TightMuonCuts = PassMuonTightPtCut(&muons_[i]) && PassMuonGlobalCut(&muons_[i]) && PassMuonNpixhitCut(&muons_[i]) && PassMuonNtrkhitCut(&muons_[i]) && PassMuonHitsUsedCut(&muons_[i]) && PassMuonStationsCut(&muons_[i]) && PassMuonEtaCut(&muons_[i]) && PassMuonDxyCut(&muons_[i]);
    if (TightMuonCuts)
      tightMuons_.push_back(muons_[i]);

  }

  if (debugme){
    printf("    Contains: %i looseMuons(s)\n",
           (int)looseMuons_.size());

    printf("    Contains: %i tightMuons(s)\n",
           (int)tightMuons_.size());
  }

  if (looseMuons_.size() > 0)
  {
    sort(looseMuons_.begin(), looseMuons_.end(), highestMuonPt());
    if(debugme)
    {
      for (uint k=0; k<looseMuons_.size(); k++)
	    {
	      cout << "looseMuons_ " << k << " has pt of " << looseMuons_.at(k).pt() <<  endl;
	    }
    }
  }
  if (tightMuons_.size() > 0)
  {
    sort(tightMuons_.begin(), tightMuons_.end(), highestMuonPt());
  }

  //Fill histos for all tight muons in the event
  FillTightMuonHists();

  //Fill histos for all loose muons in the event
  FillLooseMuonHists();

  //////////////////////
  ////Deal With Jets////
  //////////////////////

  const vector<pat::Jet     > patJets      = getProduct<vector<pat::Jet     > >(event, jetsLabel_);
  if(patJets.size() < 1){
    cout << "Not enough jets. Bad bad event, returning now..." << endl;
    return;
  }
  if(debugme)
    printf("    Contains: %i pat jet(s)\n",
           (int)patJets.size());

  // Loop over jets, and see if they pass the jet criteria
  for (size_t i = 0; i < patJets.size(); ++i) {
    if (PassJetCut(&patJets[i]) && !Overlap(patJets[i], looseMuons_))
      jets_.push_back(patJets[i]);
  }
  
  FillJetMultiplicityHists();

  for (size_t nJ=0; nJ<jets_.size(); nJ++)
  {
    h_jets_pt->Fill(jets_.at(nJ).pt(), weight_);
    h_jets_eta->Fill(jets_.at(nJ).eta(), weight_);
    h_jets_phi->Fill(jets_.at(nJ).phi(), weight_);
  }

  //Fill jet histos for jets who passes the criteria
  if (jets_.size() > 0)
  {
    sort(jets_.begin(), jets_.end(), highestJetPt());
    if(debugme)
    {
      for (uint k=0; k<jets_.size(); k++)
	    {
	      cout << "jets_ " << k << " has pt of " << jets_.at(k).pt() <<  endl;
	    }
    }
  }


  if (jets_.size() > 0)
  {
    h_jet1_pt->Fill(jets_.at(0).pt(), weight_);
    h_jet1_eta->Fill(jets_.at(0).eta(), weight_);    
    h_jet1_phi->Fill(jets_.at(0).phi(), weight_);
    h_jet1_mass->Fill(jets_.at(0).mass(), weight_);

    if (looseMuons_.size() >= 1)
      h_deltaR_jet1muon1->Fill(reco::deltaR(jets_.at(0), looseMuons_.at(0)), weight_);
    
    if (looseMuons_.size() >= 2)
      h_deltaR_jet1muon2->Fill(reco::deltaR(jets_.at(0), looseMuons_.at(1)), weight_);
	  
      
    if (jets_.size() > 1)
    {
      h_jet2_pt->Fill(jets_.at(1).pt(), weight_);
      h_jet2_eta->Fill(jets_.at(1).eta(), weight_);    
      h_jet2_phi->Fill(jets_.at(1).phi(), weight_);
      h_jet2_mass->Fill(jets_.at(1).mass(), weight_);

      if (looseMuons_.size() >= 1)
	      h_deltaR_jet2muon1->Fill(reco::deltaR(jets_.at(1), looseMuons_.at(0)), weight_);
	    
      if (looseMuons_.size() >= 2)
	      h_deltaR_jet2muon2->Fill(reco::deltaR(jets_.at(1), looseMuons_.at(1)), weight_);
	    
  
      reco::CompositeCandidate j1j2;
      j1j2.addDaughter(jets_.at(0));
      j1j2.addDaughter(jets_.at(1));
      AddFourMomenta addFM;
      addFM.set(j1j2);
      h_jet1jet2_mass->Fill(j1j2.mass(), weight_);
    }
  }

  if (debugme){
    printf("    Contains: %i jets(s)\n",
           (int)jets_.size());
  }

  //////////////////////////////////////////////
  /// Start Applying Cuts///////////////////////
  //////////////////////////////////////////////

  //////////////////////////////
  ///Hadronic V/////////////////
  //////////////////////////////
  int iCut=0;
  //NoCuts
  Tabulate_Me(iCut, weight_);

  //Cory: CutMinNLooseLep
  if( !PassMinNLooseLeptonsCut() ){
    return;
  }
  Tabulate_Me(iCut, weight_);

  //Cory: CutMinNJets
  if( !PassMinNJetsCut() ){
    return;
  }
  Tabulate_Me(iCut, weight_);

  // Make a W candidate out of the jets.
  //Cory: CutValidW
  wCand_ = getWCand(jets_);
  if (debugme)
    cout << "Made wCand" << endl;

  //Cory: CutValidV
  if( !PassValidHadVCut() ) return;
  Tabulate_Me(iCut, weight_);
  if (debugme) cout << "Passed vCand" << endl;

  //Cory: CutVMass
  if( !PassHadVMassCut() ) return;
  Tabulate_Me(iCut, weight_);
  if (debugme) cout << "Passed vCand Mass" << endl;

  //Cory: CutVpt
  if( !PassHadVptCut() ) return;
  Tabulate_Me(iCut, weight_);
  if (debugme) cout << "Good vCand Pt" << endl;

  if (debugme) cout << "Good V from jet" << endl;
  FillGoodHadVHistos();

  //////////////////////////////
  ///Z With Loose Muons/////////
  //////////////////////////////
  
  // Make a Z candidate out of the loose muons. 
  zCand_ = MakeZs(looseMuons_);
  if (debugme)
    cout << "Made zCand" << endl;

  if( !PassValidZCut() ) return;
  Tabulate_Me(iCut, weight_);
  if (debugme)
    cout << "Valid Z from muon" << endl;

  //Cory: CutZMass
  bool goodZ = true;
  if( !PassZMassCut() ){ goodZ = false; }
  Tabulate_Me(iCut, weight_);

  //Cory: CutZpt
  if( !PassZptCut  () ){ goodZ = false;}
  Tabulate_Me(iCut, weight_);

  if (debugme && goodZ)
    cout << "Good Z from muon" << endl;
  if (goodZ) FillGoodZHistos();

  ///////////////////////////////////////
  //// Make VZ Candidate With Loose Muons
  ///////////////////////////////////////
  if (debugme) 
    cout << "Im just before boson cands" << endl;
  
  //Cory: CutValidVZ
  if(goodZ) {
    hadVZ_ = VZCandidate(zCand_, wCand_);
    if(debugme) cout << "Made my hadVZ" << endl;  
  }
  else{
    if (debugme)
      cout << "Didnt manage to make a HadVZ Candidate" << endl;
  }

  if( !PassValidHadVZCandCut()){
    //Do Nothing
  }else
    FillValidVZHistos();
  Tabulate_Me(iCut, weight_);


  //////////////////////////////
  ///Z With Tight Muons/////////
  //////////////////////////////

  //Cory: CutMinNTightLep
  if( PassMinNTightLeptonsCut() ){ 
    if(debugme) cout << "Passed # tight leptons cut" << endl;
  }else
    if(debugme) cout << "Didnt pass # tight leptons cut" << endl;
  Tabulate_Me(iCut, weight_);

  //Make a Z cand out of tight muons
  //Cory: CutValidTightZ
  zCandTight_ = MakeZs(tightMuons_);
  if (debugme)
    cout << "Made zCandTight" << endl;

  bool validZTight = PassValidZCutTight();
  if (debugme && validZTight)
    cout << "Valid Z from tight muons" << endl;
  
  //Cory: CutTightZMass
  //Cory: CutTightZpt
  bool goodZTight = true;
  if( !validZTight || !PassZMassCutTight()){ goodZTight = false; }
  Tabulate_Me(iCut, weight_);
  if( !validZTight || !PassZptCutTight  ()){ goodZTight = false; }
  Tabulate_Me(iCut, weight_);
  
  if (debugme && goodZTight)
    cout << "Good Z from tight muons" << endl;
  if (goodZTight) FillGoodZTightHistos();
  
  ///////////////////////////////////////
  //// Make VZ Candidate With Tight Muons
  ////////////////////////////////////////

  //Cory: CutValidTightVZ
  if(goodZTight) {
    hadVZTight_ = VZCandidate(zCandTight_, wCand_);
    if(debugme) cout << "Made my hadVZ tight" << endl;
    
    if (debugme)
    {
      cout << "Filled tight muons histos from HadVZ" << endl;
      cout << "wCand pT is " << wCand_.pt() << endl;
    }
  }else{
    if (debugme)
      cout << "Didnt manage to make a HadVZ Candidate with tight muons" << endl;
  }

  if( !PassValidHadVZTightCandCut() ){
    //Do Nothing
  }else
    FillGoodVZHistos();

  Tabulate_Me(iCut, weight_);

  //AllCuts
  Tabulate_Me(iCut, weight_);

}//End of Event Loop


/////////////////Accessors///////////////////////

/////////////////Modifies///////////////////////
void HadronicVZAnalyzer::SetCandEvtFile(string s){
  outCandEvt.open(s.c_str());
  WPrimeUtil::CheckStream(outCandEvt, s);
}

ZCandidate HadronicVZAnalyzer::MakeZs(const MuonV & muons){
  //Cory: Is this opp charge needed?
  bool oppCharge=false;
  for (size_t i=0; i<muons.size(); i++){
    if (debugme) 
      cout << "muons has charge "  << muons.at(i).charge() << endl;
    for (size_t j=i+1; j<muons.size(); j++){
      if (muons[i].charge() != muons[j].charge())
        oppCharge=1;
    }
  }
  if (oppCharge== true){
    ZCandV zCands = getZCands(muons, 25.0);
    if(zCands.size()>0)
    { 
      if (debugme)
      {
        cout << "Daugh 1 pT is " << zCands[0].daughter(0)->pt() << endl;
        cout << "Daugh 2 pT is " << zCands[0].daughter(1)->pt() << endl;
      }
      return  zCands[0];
    }
  }
  return ZCandidate();
}
  
/////////////////Cuts///////////////////////

//Always true
bool HadronicVZAnalyzer::PassNoCut()
{
  return true;
}

//Trigger requirements
//-----------------------------------------------------------
bool HadronicVZAnalyzer::PassTriggersCut()
{
//-----------------------------------------------------------
  if(debugme) cout<<"Trigger requirements"<<endl;
  // To be implemented
  //return WPrimeUtil::PassTriggersCut(triggerEvent_,triggersToUse_);
  return true;
}//--- PassTriggersCut()

inline bool
HadronicVZAnalyzer::PassMinNLooseLeptonsCut(){
  return (looseMuons_.size()) > minNLeptons;
}

inline bool
HadronicVZAnalyzer::PassMinNTightLeptonsCut(){
  return tightMuons_.size() > minNLeptons;
}

bool
HadronicVZAnalyzer::PassNJetsCut(){
  return ((jets_.size() > minNJets) && (jets_.size() < maxNJets));
}

inline bool
HadronicVZAnalyzer::PassMinNJetsCut(){
  return jets_.size() > minNJets;
}

bool
HadronicVZAnalyzer::PassValidHadVCut(){
  return wCand_ && wCand_.mass()>0.;
}

bool
HadronicVZAnalyzer::PassValidZCut(){
  return zCand_ && zCand_.mass()>0.;
}

bool
HadronicVZAnalyzer::PassValidZCutTight(){
  return zCandTight_ && zCandTight_.mass()>0.;
}

bool
HadronicVZAnalyzer::PassValidHadVZCandCut(){
  return hadVZ_ && hadVZ_.mass()>0.;
}

bool
HadronicVZAnalyzer::PassValidHadVZTightCandCut(){
  return hadVZTight_ && hadVZTight_.mass()>0.;
}

bool
HadronicVZAnalyzer::PassNumberOfZsCut(){
  return numZs < maxNumZs;
}

bool
HadronicVZAnalyzer::PassLeadingLeptonPtCut(){
  return LeadPt > minLeadPt;
}

////////////////////////////////
/////////Check Z Properties/////
////////////////////////////////
bool
HadronicVZAnalyzer::PassZMassCut(){
  return (zCand_.mass() > minZmass) && (zCand_.mass() < maxZmass);  
}

bool
HadronicVZAnalyzer::PassZptCut(){
  return zCand_.pt() > minZpt;
}

bool
HadronicVZAnalyzer::PassZMassCutTight(){
  return (zCandTight_.mass() > minZmass) && (zCandTight_.mass() < maxZmass);  
}

bool
HadronicVZAnalyzer::PassZptCutTight(){
  return zCandTight_.pt() > minZpt;
}
////////////////////////////////
///////Check Had. V Properties//
////////////////////////////////

bool HadronicVZAnalyzer::PassHadVMassCut(){
  return (wCand_.mass() > minHadVmass) && (wCand_.mass() < maxHadVmass);
}

bool
HadronicVZAnalyzer::PassHadVptCut(){
  return wCand_.pt() > minHadVpt;
}

////////////////////////////////
//////Check Muon Properties/////
////////////////////////////////
/// Big fution to check all cuts.
bool HadronicVZAnalyzer::PassMuonCut(const TeVMuon* mu){
  for(uint i=0; i<MuonCutFns_.size(); ++i){
    bool result = CALL_MEMBER_FN(*this,MuonCutFns_[i])(mu);
    if(result==false) return false;
  }
  return true;
}

bool HadronicVZAnalyzer::PassMuonLoosePtCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Loose Pt Cut"<<endl;
  return (mu->pt() > minMuonLoosePt);
}

bool HadronicVZAnalyzer::PassMuonTightPtCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Tight Pt Cut"<<endl;
  return (mu->pt() > minMuonTightPt);
}//--- PassMuonPtCut

bool HadronicVZAnalyzer::PassMuonGlobalCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Global Cut"<<endl;
  return (mu->isGlobalMuon()); 
}//--- PassMuonGlobalCut

bool HadronicVZAnalyzer::PassMuonNpixhitCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon NpixhitCut"<<endl;
  return (mu->globalTrack()->hitPattern().numberOfValidPixelHits() > minMuonNPixHit);
}//--- PassMuonNpixhitCut

bool HadronicVZAnalyzer::PassMuonNtrkhitCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon NtrkhitCut"<<endl;
  return (mu->globalTrack()->hitPattern().numberOfValidTrackerHits() > minMuonNTrkHit);
}//--- PassMuonNtrkhitCut

bool HadronicVZAnalyzer::PassMuonNormChi2Cut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Chi2 Cut"<<endl;
  return (mu->globalTrack()->normalizedChi2() < maxMuonNormChi2);
}//--- PassMuonChi2Cut

bool HadronicVZAnalyzer::PassMuonHitsUsedCut(const TeVMuon* mu){
  //Num Valid Muon Hits
  if(debugme) cout<<"Check Muon Hits Used Cut"<<endl;
  return (mu->globalTrack()->hitPattern().numberOfValidMuonHits() > minMuonHitsUsed);
}//--- PassMuonHits Used Cut

bool HadronicVZAnalyzer::PassMuonStationsCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Stations Cut"<<endl;
  return (mu->numberOfMatches() > minMuonStations);
}//--- PassMuonStationsCut

bool HadronicVZAnalyzer::PassMuonEtaCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Eta Cut"<<endl;
  return (fabs(mu->eta()) < maxMuonEta);
}//--- PassMuonEta Cut

/*
  bool HadronicVZAnalyzer::PassMuonCombRelIsoCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon CombRelIso Cut"<<endl;
  return (Calc_MuonRelIso(mu) < maxWmunuCombRelIso);
  }//--- PassMuonCombRelIsoCut
*/

bool HadronicVZAnalyzer::PassMuonDxyCut(const TeVMuon* mu){
  if(debugme) cout<<"Check Muon Dxy Cut"<<endl;
  return (fabs(mu->userFloat("d0")) < maxMuonDxy);
}//--- PassMuonDxyCut

////////////////////////////////
/////Check jet Properties //////
////////////////////////////////
bool HadronicVZAnalyzer::PassJetCut(const pat::Jet* jet){
  bool jetKinStatus = PassJetPtCut(jet) && PassJetEtaCut(jet);
  bool jetIDStatus = PassJetIDCut(jet);
  return jetKinStatus && jetIDStatus;
}

bool HadronicVZAnalyzer::PassJetPtCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet Pt Cut"<<endl;
  return (jet->pt() > minJetPt);
}

bool HadronicVZAnalyzer::PassJetEtaCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet Eta Cut"<<endl;
  return (fabs(jet->eta()) < maxJetEta);
}

bool HadronicVZAnalyzer::PassJetNHFCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet NHF Cut"<<endl;
  return (jet->neutralHadronEnergyFraction() < maxJetNHF);
}

bool HadronicVZAnalyzer::PassJetNEFCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet NEF Cut"<<endl;
  return (jet->neutralEmEnergyFraction() < maxJetNEF);
}

bool HadronicVZAnalyzer::PassJetNConstCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet NConst Cut"<<endl;
  return (jet->numberOfDaughters() > minJetnumConst);
}

bool HadronicVZAnalyzer::PassJetCHFCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet CHF Cut"<<endl;
  return (jet->chargedHadronEnergyFraction() > minJetCHF);
}

bool HadronicVZAnalyzer::PassJetCMultCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet CMult Cut"<<endl;
  return ((unsigned int)jet->chargedMultiplicity() > minJetcMult);
}

bool HadronicVZAnalyzer::PassJetCEFCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Jet CEF Cut"<<endl;
  return (jet->chargedEmEnergyFraction() < maxJetCEF);
}

bool HadronicVZAnalyzer::PassJetIDCut(const pat::Jet* jet){
  if(debugme) cout<<"Check Pass JetID Cut"<<endl;
  bool neutralStatus = 
    PassJetNHFCut(jet) &&
    PassJetNEFCut(jet) &&
    PassJetNConstCut(jet);
  bool chargedStatus = (fabs(jet->eta()) > 2.4) ||
    (PassJetCHFCut(jet) &&
     PassJetCEFCut(jet) &&//Typo?, should be CEF?
     PassJetCMultCut(jet)
      );
  return (neutralStatus && chargedStatus);
}

///////////////Utilities//////////////////

//Tabulate results after the cut has been passed
//-----------------------------------------------------------
void HadronicVZAnalyzer::Tabulate_Me(int& cut_index, const float& weight)
{
//-----------------------------------------------------------
  if(debugme) cout<<"Tabulating results for cut_index = "
                  <<cut_index<<" = "<<Cuts_[cut_index]<<endl;

//increase the number of events passing the cuts
  //hNumEvts->Fill(cut_index,weight);
  
  results_[cut_index].Nsurv_evt_cut_w += weight;
  results_[cut_index].Nsurv_evt_cut++;
//fill the histograms
  Fill_Histos(cut_index,weight);
    
  ++cut_index;
}//Tabulate_Me


void
HadronicVZAnalyzer::ClearEvtVariables(){
  muons_.clear();
  jets_.clear();
  looseMuons_.clear();
  tightMuons_.clear();
  zCand_ = ZCandidate();
  zCandTight_ = ZCandidate();
  wCand_ = WCandidate();
  hadVZ_ = VZCandidate();
  hadVZTight_ = VZCandidate();
}

void HadronicVZAnalyzer::beginFile(std::vector<wprime::InputFile>::const_iterator fi){
  TFileDirectory dir = wprimeUtil_->getFileService()->mkdir(fi->samplename); 
  Declare_Histos(dir);
  //ResetCounters();
}

// operations to be done when closing input file 
// (e.g. print summary)
void HadronicVZAnalyzer::endFile(std::vector<wprime::InputFile>::const_iterator fi,
                                 ofstream & out){
  WPrimeUtil::tabulateSummary(results_);
  WPrimeUtil::printSummary(fi->samplename, fi->description, Cuts_, results_, out);  
}

void HadronicVZAnalyzer::endAnalysis(ofstream & out){
}




