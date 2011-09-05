#include "UserCode/CMGWPrimeGroup/interface/HadronicVZAnalyzer.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember)) 

using namespace std;
HadronicVZAnalyzer::HadronicVZAnalyzer(){}
HadronicVZAnalyzer::HadronicVZAnalyzer(const edm::ParameterSet & cfg, WPrimeUtil * wprimeUtil) :
  AnalyzerBase(cfg, wprimeUtil){
  //FillCutFns();
  if(debugme) printf("Using %i cuts\n",NCuts_);

// +++++++++++++++++++Event characteristics
 
// +++++++++++++++++++General Cut values
  maxNJets = cfg.getParameter<uint>("maxNJets");
  maxAngleBetweenJets = cfg.getParameter<double>("maxAngleBetweenJets");

// +++++++++++++++++++Z Cuts

// +++++++++++++++++++Hadronic Boson Cuts

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
  outCandEvt_.close(); 
}


/// Declare Histograms
//--------------------------------------------------------------
void HadronicVZAnalyzer::Declare_Histos(const TFileDirectory & dir)
{
  // Extend later.
  printf("Declare histos\n");
  //Loose histos
  h_HadVZMass = dir.make<TH1F>("h_HadVZMass","h_HadVZMass",100,0.0,2500.0);
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

  //
  h_tight_muons_pt = dir.make<TH1F>("h_tight_muons_pt", "h_tight_muons_pt", 100, 0.0, 1000.);
  h_tight_muons_eta = dir.make<TH1F>("h_tight_muons_eta", "h_tight_muons_eta", 40, -5.0, 5.);  
  h_tight_muons_phi = dir.make<TH1F>("h_tight_muons_phi", "h_tight_muons_phi", 20, -4.0, 4.);

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

  //Jet merging histos
  h_m1_vs_m12 = dir.make<TH2F>("h_m1_vs_m12", "Mass Jet1 vs Mass Jet12;M_{j1};M_{j1,j2}", 100, 0., 200., 100, 0., 200.);
  h_bestmass = dir.make<TH1F>("h_bestmass", "Best V Mass;M_{V}^{Best}", 75, 0., 150.);//Cory: Change back

  h_deltaR_jet1jet2 = dir.make<TH1F>("h_deltaR_jet1jet2", "h_deltaR_jet1jet2", 50, 0.0, 5.0);
  h_deltaR_jet1Z_R1 = dir.make<TH1F>("h_deltaR_jet1Z_R1", "h_deltaR_jet1Z_R1", 50, 0.0, 5.0);
  h_deltaR_jet2Z_R1 = dir.make<TH1F>("h_deltaR_jet2Z_R1", "h_deltaR_jet2Z_R1", 50, 0.0, 5.0);
  h_deltaR_jet3Z_R1 = dir.make<TH1F>("h_deltaR_jet3Z_R1", "h_deltaR_jet3Z_R1", 50, 0.0, 5.0);
  h_deltaR_jet1Z_R2 = dir.make<TH1F>("h_deltaR_jet1Z_R2", "h_deltaR_jet1Z_R2", 50, 0.0, 5.0);
  h_deltaR_jet2Z_R2 = dir.make<TH1F>("h_deltaR_jet2Z_R2", "h_deltaR_jet2Z_R2", 50, 0.0, 5.0);
  h_deltaR_jet3Z_R2 = dir.make<TH1F>("h_deltaR_jet3Z_R2", "h_deltaR_jet3Z_R2", 50, 0.0, 5.0);

  h_deltaR_jet1jet2_R1_cut40 = dir.make<TH1F>("h_deltaR_jet1jet2_R1_cut40", "h_deltaR_jet1jet2_R1_cut40", 50, 0.0, 5.0);
  h_deltaR_jet1jet2_R1_cut50 = dir.make<TH1F>("h_deltaR_jet1jet2_R1_cut50", "h_deltaR_jet1jet2_R1_cut50", 50, 0.0, 5.0);
  h_deltaR_jet1jet2_R1_cut60 = dir.make<TH1F>("h_deltaR_jet1jet2_R1_cut60", "h_deltaR_jet1jet2_R1_cut60", 50, 0.0, 5.0);
  h_deltaR_jet1jet2_R1_cut65 = dir.make<TH1F>("h_deltaR_jet1jet2_R1_cut65", "h_deltaR_jet1jet2_R1_cut65", 50, 0.0, 5.0);
  h_deltaR_jet1jet2_R1_cut70 = dir.make<TH1F>("h_deltaR_jet1jet2_R1_cut70", "h_deltaR_jet1jet2_R1_cut70", 50, 0.0, 5.0);
  h_deltaR_jet1jet2_R1_cut80 = dir.make<TH1F>("h_deltaR_jet1jet2_R1_cut80", "h_deltaR_jet1jet2_R1_cut80", 50, 0.0, 5.0);
  h_deltaR_jet1jet2_R1_cut90 = dir.make<TH1F>("h_deltaR_jet1jet2_R1_cut90", "h_deltaR_jet1jet2_R1_cut90", 50, 0.0, 5.0);

  h_deltaR_jet1jet2_R2_cut40 = dir.make<TH1F>("h_deltaR_jet1jet2_R2_cut40", "h_deltaR_jet1jet2_R2_cut40", 50, 0.0, 5.0);
  h_deltaR_jet1jet2_R2_cut50 = dir.make<TH1F>("h_deltaR_jet1jet2_R2_cut50", "h_deltaR_jet1jet2_R2_cut50", 50, 0.0, 5.0);
  h_deltaR_jet1jet2_R2_cut60 = dir.make<TH1F>("h_deltaR_jet1jet2_R2_cut60", "h_deltaR_jet1jet2_R2_cut60", 50, 0.0, 5.0);
  h_deltaR_jet1jet2_R2_cut65 = dir.make<TH1F>("h_deltaR_jet1jet2_R2_cut65", "h_deltaR_jet1jet2_R2_cut65", 50, 0.0, 5.0);
  h_deltaR_jet1jet2_R2_cut70 = dir.make<TH1F>("h_deltaR_jet1jet2_R2_cut70", "h_deltaR_jet1jet2_R2_cut70", 50, 0.0, 5.0);
  h_deltaR_jet1jet2_R2_cut80 = dir.make<TH1F>("h_deltaR_jet1jet2_R2_cut80", "h_deltaR_jet1jet2_R2_cut80", 50, 0.0, 5.0);
  h_deltaR_jet1jet2_R2_cut90 = dir.make<TH1F>("h_deltaR_jet1jet2_R2_cut90", "h_deltaR_jet1jet2_R2_cut90", 50, 0.0, 5.0);


  /////////////////
  DeclareHistoSet("hVZMass", "Reconstructed VZ Invariant Mass",
                  "M_{VZ} (GeV)", 100, 0, 2500, "GeV", hVZMass,dir);
  DeclareHistoSet("hVZpt", "Reconstructed VZ Transverse Momentum",
                  "p_{VZ}^{T} (GeV)", 100, 0, 1000, "GeV", hVZpt,dir);

  DeclareHistoSet("hZMass" , "Reconstructed Mass of Z",
                  "M_{Z} (GeV)", 30, 60, 120, "GeV", hZMass,dir);
  DeclareHistoSet("hZpt", "p_{T}^{Z}", 
                  "p_{T}^{Z} (GeV)", 100, 0, 1000, "GeV", hZpt,dir);
  DeclareHistoSet("hVMass" , "Reconstructed Mass of V",
                  "M_{V} (GeV)", 75, 0, 150, "GeV", hVMass,dir);//Cory: Change back
  DeclareHistoSet("hVpt", "p_{T}^{V}", 
                  "p_{T}^{V} (GeV)", 100, 0, 1000, "GeV", hVpt,dir);

  DeclareHistoSet("hNLLeps", "Number of Loose Leptons in Event",
                  "N_{l}^{Loose}", 10, 0, 10, "NONE", hNLLeps,dir);
  DeclareHistoSet("hNLJets", "Number of Loose Jets in Event",
                  "N_{Jets}^{Loose}", 10, 0, 10, "NONE", hNLJets,dir);


  string title = Form("Expected # of Events / %.0f pb^{-1}",  wprimeUtil_->getLumi_ipb());
  title = title + ";;" + title;
  hNumEvts = NULL; hNumEvts = dir.make<TH1F>("hNumEvts",title.c_str(),NCuts_,0, NCuts_);
  for(int i=0; i<NCuts_; ++i) hNumEvts->GetXaxis()->SetBinLabel(i+1,Cuts_[i].c_str());

  cout << "Histos declared" << endl;

}//Declare_Histos

//Fill Histograms
//-----------------------------------------------------------
void HadronicVZAnalyzer::Fill_Histos(const int& index, const float& weight)
{
  if(debugme) printf("Filling Histos\n");

  if(hadVZ_){
    hVZMass[index]->Fill(hadVZ_.mass(), weight);
    hVZpt[index]->Fill(hadVZ_.pt(), weight);
  }
  if(zCand_){
    hZMass[index]->Fill(zCand_.mass(), weight);
    hZpt[index]->Fill(zCand_.pt(), weight);
  }
  if(vCand_){
    hVMass[index]->Fill(vCand_.mass(), weight);
    hVpt[index]->Fill(vCand_.pt(), weight);
  }
  hNLLeps[index]->Fill(looseElectrons_.size()+looseMuons_.size(), weight);
  hNLJets[index]->Fill(looseJets_.size(), weight);

}//Fill_Histos

void HadronicVZAnalyzer::FillJetMultiplicityHists(){
  h_jet_mult->Fill(looseJets_.size(), weight_);
  uint max = min(10, (int)looseJets_.size());
  for(uint i=0; i<=max; ++i)
    h_jet_mult_inc->Fill(i, weight_);
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
  const TeVMuon & m1 = WPrimeUtil::Find(*zCand_.daughter(0), allMuons_);
  const TeVMuon & m2 = WPrimeUtil::Find(*zCand_.daughter(1), allMuons_);
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

void HadronicVZAnalyzer::FillGoodHadVHistos(){
  h_jet_HadV_pt->Fill(vCand_.pt(), weight_);
  h_jet_HadV_eta->Fill(vCand_.eta(), weight_);
  h_jet_HadV_phi->Fill(vCand_.phi(), weight_);
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
  const TeVMuon & VZm1 = WPrimeUtil::Find(*zCand_.daughter(0), allMuons_);
  const TeVMuon & VZm2 = WPrimeUtil::Find(*zCand_.daughter(1), allMuons_);
  //cout << "Muon from loose zCand" << endl;
  h_Zmuon1_VZCut_pt->Fill(VZm1.pt(), weight_);
  h_Zmuon1_VZCut_eta->Fill(VZm1.eta(), weight_);
  h_Zmuon1_VZCut_phi->Fill(VZm1.phi(), weight_);
  h_Zmuon2_VZCut_pt->Fill(VZm2.pt(), weight_);
  h_Zmuon2_VZCut_eta->Fill(VZm2.eta(), weight_);
  h_Zmuon2_VZCut_phi->Fill(VZm2.phi(), weight_);
  
  h_jet_VZCut_pt->Fill(vCand_.pt(), weight_);
  h_jet_VZCut_eta->Fill(vCand_.eta(), weight_);
  h_jet_VZCut_phi->Fill(vCand_.phi(), weight_);
  
  h_deltaR_HadVmuon1->Fill(reco::deltaR(vCand_, VZm1), weight_);
  h_deltaR_HadVmuon2->Fill(reco::deltaR(vCand_, VZm2), weight_);
}

void HadronicVZAnalyzer::FillJetMergingHistos(){

  if (allJets_.size()>1){
    h_deltaR_jet1jet2->Fill(reco::deltaR(allJets_.at(0), allJets_.at(1)), weight_);
    
    //Region histos
    if (allJets_.at(0).mass() < 40){
      h_deltaR_jet1jet2_R1_cut40->Fill(reco::deltaR(allJets_.at(0), allJets_.at(1)), weight_);
    }
    if (allJets_.at(0).mass() > 40){
      h_deltaR_jet1jet2_R2_cut40->Fill(reco::deltaR(allJets_.at(0), allJets_.at(1)), weight_);  
    }

    if (allJets_.at(0).mass() < 50){
      h_deltaR_jet1jet2_R1_cut50->Fill(reco::deltaR(allJets_.at(0), allJets_.at(1)), weight_);
    }
    if (allJets_.at(0).mass() > 50){
      h_deltaR_jet1jet2_R2_cut50->Fill(reco::deltaR(allJets_.at(0), allJets_.at(1)), weight_);  
    }

    if (allJets_.at(0).mass() < 60){
      h_deltaR_jet1jet2_R1_cut60->Fill(reco::deltaR(allJets_.at(0), allJets_.at(1)), weight_);
    }
    if (allJets_.at(0).mass() > 60){
      h_deltaR_jet1jet2_R2_cut60->Fill(reco::deltaR(allJets_.at(0), allJets_.at(1)), weight_);  
    }

    if (allJets_.at(0).mass() < 65){
      h_deltaR_jet1jet2_R1_cut65->Fill(reco::deltaR(allJets_.at(0), allJets_.at(1)), weight_);
    }
    if (allJets_.at(0).mass() > 65){
      h_deltaR_jet1jet2_R2_cut65->Fill(reco::deltaR(allJets_.at(0), allJets_.at(1)), weight_);  
    }

    if (allJets_.at(0).mass() < 70){
      h_deltaR_jet1jet2_R1_cut70->Fill(reco::deltaR(allJets_.at(0), allJets_.at(1)), weight_);
    }
    if (allJets_.at(0).mass() > 70){
      h_deltaR_jet1jet2_R2_cut70->Fill(reco::deltaR(allJets_.at(0), allJets_.at(1)), weight_);  
    }

    if (allJets_.at(0).mass() < 80){
      h_deltaR_jet1jet2_R1_cut80->Fill(reco::deltaR(allJets_.at(0), allJets_.at(1)), weight_);
    }
    if (allJets_.at(0).mass() > 80){
      h_deltaR_jet1jet2_R2_cut80->Fill(reco::deltaR(allJets_.at(0), allJets_.at(1)), weight_);  
    }

    if (allJets_.at(0).mass() < 90){
      h_deltaR_jet1jet2_R1_cut90->Fill(reco::deltaR(allJets_.at(0), allJets_.at(1)), weight_);
    }
    if (allJets_.at(0).mass() > 90){
      h_deltaR_jet1jet2_R2_cut90->Fill(reco::deltaR(allJets_.at(0), allJets_.at(1)), weight_);  
    }
    

    if (allJets_.at(0).mass() < 65){
      h_deltaR_jet2Z_R1->Fill(reco::deltaR(allJets_.at(1), zCand_), weight_);
      if (allJets_.size()>2)
	h_deltaR_jet3Z_R1->Fill(reco::deltaR(allJets_.at(2), zCand_), weight_);
    }

    if (allJets_.at(0).mass() > 65){
      h_deltaR_jet2Z_R2->Fill(reco::deltaR(allJets_.at(1), zCand_), weight_);
      if (allJets_.size()>2)
	h_deltaR_jet3Z_R2->Fill(reco::deltaR(allJets_.at(2), zCand_), weight_);
    }

  }


  if (allJets_.at(0).mass() < 65){
    if (allJets_.size()>0)
      h_deltaR_jet1Z_R1->Fill(reco::deltaR(allJets_.at(0), zCand_), weight_);
  }

  if (allJets_.at(0).mass() > 65){
    if (allJets_.size()>0)
      h_deltaR_jet1Z_R2->Fill(reco::deltaR(allJets_.at(0), zCand_), weight_);
  }



}

void 
HadronicVZAnalyzer::eventLoop(edm::EventBase const & event){
  if(debugme) WPrimeUtil::PrintEvent(event);
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
    allMuons_.push_back(TeVMuon(patMuons[i],muonAlgo_));   

    //  bool LooseMuonCuts = PassMuonLoosePtCut(&muons_[i]) && PassMuonGlobalCut(&muons_[i]) && PassMuonNpixhitCut(&muons_[i]) && PassMuonNtrkhitCut(&muons_[i]) && PassMuonHitsUsedCut(&muons_[i]) && PassMuonStationsCut(&muons_[i]) && PassMuonEtaCut(&muons_[i]) && PassMuonDxyCut(&muons_[i]);
    //    bool LooseMuonCuts = PassMuonLoosePtCut(&muons_[i]) && PassMuonGlobalCut(&muons_[i]) && PassMuonNpixhitCut(&muons_[i]) && PassMuonNtrkhitCut(&muons_[i]) && PassMuonHitsUsedCut(&muons_[i]) && PassMuonStationsCut(&muons_[i]) && PassMuonEtaCut(&muons_[i]) && PassMuonDxyCut(&muons_[i]);
    
    bool LooseMuonCuts = PassMuonLoosePtCut(&allMuons_[i]) && PassMuonGlobalCut(&allMuons_[i]) && PassMuonNpixhitCut(&allMuons_[i]) && PassMuonHitsUsedCut(&allMuons_[i]) && PassMuonStationsCut(&allMuons_[i]) && PassMuonEtaCut(&allMuons_[i]) && PassMuonDxyCut(&allMuons_[i]);
    if (LooseMuonCuts)
      looseMuons_.push_back(allMuons_[i]);
    
    bool TightMuonCuts = PassMuonTightPtCut(&allMuons_[i]) && PassMuonGlobalCut(&allMuons_[i]) && PassMuonNpixhitCut(&allMuons_[i]) && PassMuonNtrkhitCut(&allMuons_[i]) && PassMuonHitsUsedCut(&allMuons_[i]) && PassMuonStationsCut(&allMuons_[i]) && PassMuonEtaCut(&allMuons_[i]) && PassMuonDxyCut(&allMuons_[i]);
    if (TightMuonCuts)
      tightMuons_.push_back(allMuons_[i]);

  }

  if (debugme){
    printf("    Contains: %i looseMuons(s)\n",
           (int)looseMuons_.size());

    printf("    Contains: %i tightMuons(s)\n",
           (int)tightMuons_.size());
  }

  if (looseMuons_.size() > 0){
    sort(looseMuons_.begin(), looseMuons_.end(), highestMuonPt());
    if(debugme){
      for (uint k=0; k<looseMuons_.size(); k++){
	      cout << "looseMuons_ " << k << " has pt of " << looseMuons_.at(k).pt() <<  endl;
	    }
    }
  }
  if (tightMuons_.size() > 0){
    sort(tightMuons_.begin(), tightMuons_.end(), highestMuonPt());
  }

  //Fill histos for all tight muons in the event
  FillTightMuonHists();

  //Fill histos for all loose muons in the event
  FillLooseMuonHists();

  //////////////////////
  ////Deal With Jets////
  //////////////////////

  allJets_      = getProduct<vector<pat::Jet     > >(event, jetsLabel_);
  if(allJets_.size() < 1){
    cout << "Not enough jets. Bad bad event, returning now..." << endl;
    return;
  }
  if(debugme)
    printf("    Contains: %i pat jet(s)\n",
           (int)allJets_.size());

  // Loop over jets, and see if they pass the jet criteria
  for (size_t i = 0; i < allJets_.size(); ++i) {
    if (PassJetCut(&allJets_[i]) && !Overlap(allJets_[i], looseMuons_, 0.5, 2))
      looseJets_.push_back(allJets_[i]);
  }
  
  FillJetMultiplicityHists();

  for (size_t nJ=0; nJ<looseJets_.size(); nJ++)
  {
    h_jets_pt->Fill(looseJets_.at(nJ).pt(), weight_);
    h_jets_eta->Fill(looseJets_.at(nJ).eta(), weight_);
    h_jets_phi->Fill(looseJets_.at(nJ).phi(), weight_);
  }

  //Fill jet histos for jets who passes the criteria
  if (looseJets_.size() > 0){
    sort(looseJets_.begin(), looseJets_.end(), highestJetPt());
    if(debugme){
      for (uint k=0; k<looseJets_.size(); k++){
	      cout << "jets_ " << k << " has pt of " << looseJets_.at(k).pt() <<  endl;
	    }
    }
  }


  if (looseJets_.size() > 0){
    h_jet1_pt->Fill(looseJets_.at(0).pt(), weight_);
    h_jet1_eta->Fill(looseJets_.at(0).eta(), weight_);    
    h_jet1_phi->Fill(looseJets_.at(0).phi(), weight_);
    h_jet1_mass->Fill(looseJets_.at(0).mass(), weight_);

    if (looseMuons_.size() >= 1)
      h_deltaR_jet1muon1->Fill(reco::deltaR(looseJets_.at(0), looseMuons_.at(0)), weight_);
    
    if (looseMuons_.size() >= 2)
      h_deltaR_jet1muon2->Fill(reco::deltaR(looseJets_.at(0), looseMuons_.at(1)), weight_);
	  
      
    if (looseJets_.size() > 1)
    {
      h_jet2_pt->Fill(looseJets_.at(1).pt(), weight_);
      h_jet2_eta->Fill(looseJets_.at(1).eta(), weight_);    
      h_jet2_phi->Fill(looseJets_.at(1).phi(), weight_);
      h_jet2_mass->Fill(looseJets_.at(1).mass(), weight_);

      if (looseMuons_.size() >= 1)
	      h_deltaR_jet2muon1->Fill(reco::deltaR(looseJets_.at(1), looseMuons_.at(0)), weight_);
	    
      if (looseMuons_.size() >= 2)
	      h_deltaR_jet2muon2->Fill(reco::deltaR(looseJets_.at(1), looseMuons_.at(1)), weight_);
	    
  
      reco::CompositeCandidate j1j2;
      j1j2.addDaughter(looseJets_.at(0));
      j1j2.addDaughter(looseJets_.at(1));
      AddFourMomenta addFM;
      addFM.set(j1j2);
      h_jet1jet2_mass->Fill(j1j2.mass(), weight_);
    }
  }

  if (debugme){
    printf("    Contains: %i jets(s)\n",
           (int)looseJets_.size());
  }
  FillJetMergingHistos();

  //////////////////////////////////////////////
  /// Start Applying Cuts///////////////////////
  //////////////////////////////////////////////

  int iCut=0;
  if( !PassNoCut() ) return;
  Tabulate_Me(iCut, weight_); ++iCut;

  if( !PassMinNLeptonsCut() ) return;
  Tabulate_Me(iCut, weight_); ++iCut;

  if( !PassMinNJetsCut() ) return;
  Tabulate_Me(iCut, weight_); ++iCut;

  //////////////////////////////
  ///////  Z With Muons  ///////
  //////////////////////////////
  
  // Make a Z candidate out of the loose muons. 
  ZCandV zCands = getZCands(looseMuons_, 100.);
  zCand_ = zCands.size() ? zCands[0] : ZCandidate();
  if (debugme) cout << "Made zCand" << endl;

  if( !PassValidZCut() ) return;
  Tabulate_Me(iCut, weight_); ++iCut;

  if( !PassZMassCut() ) return;
  Tabulate_Me(iCut, weight_); ++iCut;

  if( !PassZptCut  () ) return;
  Tabulate_Me(iCut, weight_); ++iCut;

  FillGoodZHistos();

  ///////////////////////////////////////
  //////// Make V from Jets  ////////////
  ///////////////////////////////////////

  vCand_ = getWCand(looseJets_);
  if (debugme) cout << "Made vCand" << endl;

  if( !PassValidVCut() ) return;
  Tabulate_Me(iCut, weight_); ++iCut;

  LorentzVector j1j2 = looseJets_.at(0).p4();
  if(looseJets_.size()>1) j1j2 += looseJets_.at(1).p4();
  h_m1_vs_m12->Fill(looseJets_.at(0).mass(), j1j2.mass(), weight_);
  float bestMass = fabs(j1j2.mass() - 85.) < fabs(looseJets_.at(0).mass() - 85.) ? j1j2.mass() : looseJets_.at(0).mass();
  h_bestmass->Fill(bestMass, weight_);
  //if(bestMass < 60) cout<<" nJets: "<<looseJets_.size()<<" best: "<<bestMass<<" m1: "<<looseJets_.at(0).mass()<<" m12: "<<j1j2.mass()<<endl;

  if( !PassVMassCut() ) return;
  Tabulate_Me(iCut, weight_); ++iCut;

  if( !PassVptCut() ) return;
  Tabulate_Me(iCut, weight_); ++iCut;

  if (debugme) cout << "Good V from jet" << endl;
  FillGoodHadVHistos();

  ///////////////////////////////////////
  //////// Make VZ Candidate ////////////
  ///////////////////////////////////////
  if (debugme) cout << "Im just before boson cands" << endl;
  
  hadVZ_ = VZCandidate(zCand_, vCand_);
  if(debugme) cout << "Made my hadVZ" << endl;  

  if( !PassValidVZCandCut() ) return;
  Tabulate_Me(iCut, weight_); ++iCut;

  FillValidVZHistos();

  //AllCuts
  Tabulate_Me(iCut, weight_); ++iCut;

  if(wprimeUtil_->runningOnData()){
    cout<<" The following data event passed All Cuts!!!\n";
    PrintPassingEvent(event);
    if(1 || debugme){ 
      //PrintEventLeptons();
      PrintMuons();
      PrintJets();
    }
    cout<<" ------------------\n";
  }
  
}//End of Event Loop


/////////////////Accessors///////////////////////

/////////////////Modifies///////////////////////
  
/////////////////Cuts///////////////////////

bool
HadronicVZAnalyzer::PassValidVZCandCut(){
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
void HadronicVZAnalyzer::PrintEventDetails() const{
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
HadronicVZAnalyzer::ClearEvtVariables(){
  allJets_.clear();
  looseJets_.clear();
  allMuons_.clear();
  looseMuons_.clear();
  tightMuons_.clear();
  zCand_ = ZCandidate();
  vCand_ = WCandidate();
  hadVZ_ = VZCandidate();
}

