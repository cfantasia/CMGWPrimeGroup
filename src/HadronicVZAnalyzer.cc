#include "UserCode/CMGWPrimeGroup/interface/HadronicVZAnalyzer.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember)) 

using namespace std;
HadronicVZAnalyzer::HadronicVZAnalyzer(){}
HadronicVZAnalyzer::HadronicVZAnalyzer(const edm::ParameterSet & cfg, int fileToRun) :
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

HadronicVZAnalyzer::~HadronicVZAnalyzer(){
}


/// Declare Histograms
void HadronicVZAnalyzer::defineHistos(const TFileDirectory & dir){
  printf("Declare histos\n");

  //POG histos
  h_dptpt2 = dir.make<TH1F>("h_dptpt2", "h_dptpt2", 500, 0., 0.05);
  h_dptpt_vs_pt = dir.make<TH2F>("h_dptpt_vs_pt", "h_dptpt_vs_pt", 50, 0., 0.5, 50, 0.0, 500.);
  h_dptpt2_vs_pt = dir.make<TH2F>("h_dptpt2_vs_pt", "h_dptpt2_vs_pt", 500, 0., 0.05, 150, 0.0, 500.);
  h_dptpt_vs_invpt = dir.make<TH2F>("h_dptpt_vs_invpt","h_dptpt_vs_invpt", 50, 0., 0.5, 50, 0., 0.1);
  h_dptpt2_vs_invpt = dir.make<TH2F>("h_dptpt2_vs_invpt","h_dptpt2_vs_invpt", 500, 0., 0.05, 150, 0., 0.1);
 
  h_dptpt_vs_genpt = dir.make<TH2F>("h_dptpt_vs_genpt", "h_dptpt_vs_genpt", 50, 0., 0.5, 50, 0.0, 500.);
  h_dptpt2_vs_genpt = dir.make<TH2F>("h_dptpt2_vs_genpt", "h_dptpt2_vs_genpt", 500, 0., 0.05, 150, 0.0, 500.);
  h_dptpt_vs_invgenpt = dir.make<TH2F>("h_dptpt_vs_invgenpt","h_dptpt_vs_invgenpt", 50, 0., 0.5, 50, 0., 0.1);
  h_dptpt2_vs_invgenpt = dir.make<TH2F>("h_dptpt2_vs_invgenpt","h_dptpt2_vs_invgenpt", 500, 0., 0.05, 150, 0., 0.1);




  //Loose histos
  //Dealing with MET
  h_HadVWMass = dir.make<TH1F>("h_HadVWMass", "h_HadVWMass", 100, 0.0, 2.5);
  // h_MET_AllCuts = dir.make<TH1F>("h_MET_AllCuts", "h_MET_AllCuts", 100, 0.0, 1000.0);
  h_WMass = dir.make<TH1F>("h_WMass", "h_WMass", 100, 0.0, 300.0);
  h_genWMass = dir.make<TH1F>("h_genWMass", "h_genWMass", 100, 0.0, 300.0);
  h_genZMass = dir.make<TH1F>("h_genZMass", "h_genZMass", 100, 0.0, 300.0);
  h_HadVZgenMass = dir.make<TH1F>("h_HadVZgenMass", "h_HadVZgenMass", 100, 0.0, 2.5);
  h_HadVWgenMass = dir.make<TH1F>("h_HadVWgenMass", "h_HadVWgenMass", 100, 0.0, 2.5);

  h_VWMass_nJets = dir.make<TH2F>("h_VWMass_nJets", "h_VWMass_nJets", 10, -0.5, 9.5, 100, 0.0, 2500.0);
  h_VZMass_nJets = dir.make<TH2F>("h_VZMass_nJets", "h_VZMass_nJets", 10, -0.5, 9.5, 100, 0.0, 2500.0);
  h_ptJetCut_nJets = dir.make<TH2F>("h_ptJetCut_nJets", "h_ptJetCut_nJets", 10, -0.5,9.5, 5, 25.0, 75.0);

  //HadVZ Properties
  h_HadVZMass = dir.make<TH1F>("h_HadVZMass","h_HadVZMass; M_{VZ} (TeV)",100,0.0,2.5);
  h_HadVZpt = dir.make<TH1F>("h_HadVZpt", "h_HadVZptp_{T}^{VZ} (GeV)", 60, 0.0, 300.0);
  h_HadVZeta = dir.make<TH1F>("h_HadVZeta", "h_HadVZeta; #eta^{VZ}", 40, -5.0, 5.0);
  h_HadVZphi = dir.make<TH1F>("h_HadVZphi", "h_HadVZphi; #phi^{VZ}", 20, -4.0, 4.0);
  h_HadVZ_res = dir.make<TH1F>("h_HadVZ_res", "h_hadVZ_res; M_{VZ}^{Reco} - M_{VZ}^{Gen} (GeV)", 400, -200.0, 200.0); 
  //ZLeptons Properties
  h_Zelec1_pt = dir.make<TH1F>("h_Zelec1_pt", "h_Zelec1_pt; p_{T}^{e1} (GeV)", 100, 0.0, 1000.0);
  h_Zelec1_eta = dir.make<TH1F>("h_Zelec1_eta", "h_Zelec1_eta; #eta^{e1}", 40, -5.0, 5.0);
  h_Zelec1_phi = dir.make<TH1F>("h_Zelec1_phi", "h_Zelec1_phi; #phi^{e1}", 20, -4.0, 4.0);
  h_Zelec2_pt = dir.make<TH1F>("h_Zelec2_pt", "h_Zelec2_pt; p_{T}^{e2} (GeV)", 100, 0.0, 1000.0);
  h_Zelec2_eta = dir.make<TH1F>("h_Zelec2_eta", "h_Zelec2_eta; #eta^{e2}", 40, -5.0, 5.0);
  h_Zelec2_phi = dir.make<TH1F>("h_Zelec2_phi", "h_Zelec2_phi; #phi^{e2}", 20, -4.0, 4.0);
  h_Zmuon1_pt = dir.make<TH1F>("h_Zmuon1_pt", "h_Zmuon1_pt; p_{T}^{#mu1} (GeV)", 100, 0.0, 1000.0);
  h_Zmuon1_eta = dir.make<TH1F>("h_Zmuon1_eta", "h_Zmuon1_eta; #eta^{#mu1}", 40, -5.0, 5.0);
  h_Zmuon1_phi = dir.make<TH1F>("h_Zmuon1_phi", "h_Zmuon1_phi; #phi^{#mu1}", 20, -4.0, 4.0);
  h_Zmuon2_pt = dir.make<TH1F>("h_Zmuon2_pt", "h_Zmuon2_pt; p_{T}^{#mu2} (GeV)", 100, 0.0, 1000.0);
  h_Zmuon2_eta = dir.make<TH1F>("h_Zmuon2_eta", "h_Zmuon2_eta; #eta^{#mu2}", 40, -5.0, 5.0);
  h_Zmuon2_phi = dir.make<TH1F>("h_Zmuon2_phi", "h_Zmuon2_phi; #phi^{#mu2}", 20, -4.0, 4.0);
  h_deltaR_elec1elec2 = dir.make<TH1F>("h_deltaR_elec1elec2", "h_deltaR_elec1elec2; #Delta_{R}^{e1,e2}", 50, 0., 5.);
  h_deltaR_muon1muon2 = dir.make<TH1F>("h_deltaR_muon1muon2", "h_deltaR_muon1muon2; #Delta_{R}^{#mu1,#mu2}", 50, 0., 5.);

  //HadV vs Leptons properties
  h_deltaR_HadVelec1 = dir.make<TH1F>("h_deltaR_HadVelec1", "h_deltaR_HadVelec1; #Delta_{R}^{V,e1}", 50, 0., 5.);
  h_deltaR_HadVelec2 = dir.make<TH1F>("h_deltaR_HadVelec2", "h_deltaR_HadVelec2;\
 #Delta_{R}^{V,e2}", 50, 0., 5.);
  h_deltaR_HadVmuon1 = dir.make<TH1F>("h_deltaR_HadVmuon1", "h_deltaR_HadVmuon1;\
 #Delta_{R}^{V,#mu1}", 50, 0., 5.);
  h_deltaR_HadVmuon2 = dir.make<TH1F>("h_deltaR_HadVmuon2", "h_deltaR_HadVmuon2;\
 #Delta_{R}^{V,#mu2}", 50, 0., 5.);
  
  //Jet Properties
  h_jet1_pt = dir.make<TH1F>("h_jet1_pt", "h_jet1_pt; p_{T}^{jet1} (GeV)", 100, 0.0, 1000.0);
  h_jet1_eta = dir.make<TH1F>("h_jet1_eta", "h_jet1_eta; #eta^{jet1}", 40, -5.0, 5.0);
  h_jet1_phi = dir.make<TH1F>("h_jet1_phi", "h_jet1_phi; #phi^{jet1}", 20, -4.0, 4.0);
  h_jet2_pt = dir.make<TH1F>("h_jet2_pt", "h_jet2_pt; p_{T}^{jet2} (GeV)", 100, 0.0, 1000.0);
  h_jet2_eta = dir.make<TH1F>("h_jet2_eta", "h_jet2_eta; #eta^{jet2}", 40, -5.0, 5.0);
  h_jet2_phi = dir.make<TH1F>("h_jet2_phi", "h_jet2_phi; #phi^{jet2}", 20, -4.0, 4.0);
  h_jet1_mass = dir.make<TH1F>("h_jet1_mass", "h_jet1_mass; M_{jet1} (GeV)", 100, 0.0, 500.);
  h_jet2_mass = dir.make<TH1F>("h_jet2_mass", "h_jet2_mass; M_{jet2} (GeV)", 100, 0.0, 500.);
  // h_jet1_Vwindow = dir.make<TH1F>("h_jet1_Vwindow", "h_jet1_Vwindow", 100, 0., 500.);
  h_deltaR_jet1muon1 = dir.make<TH1F>("h_deltaR_jet1muon1", "h_deltaR_jet1muon1; #Delta_{R}^{jet1,#mu1}", 50, 0., 5.);
  h_deltaR_jet1muon2 = dir.make<TH1F>("h_deltaR_jet1muon2", "h_deltaR_jet1muon2; #Delta_{R}^{jet1,#mu2}", 50, 0., 5.);
  h_deltaR_jet2muon1 = dir.make<TH1F>("h_deltaR_jet2muon1", "h_deltaR_jet2muon1; #Delta_{R}^{jet2,#mu1}", 50, 0., 5.);
  h_deltaR_jet2muon2 = dir.make<TH1F>("h_deltaR_jet2muon2", "h_deltaR_jet2muon2; #Delta_{R}^{jet2,#mu2}", 50, 0., 5.);
  h_jet_mult = dir.make<TH1F>("h_jet_mult", "h_jet_mult; N_{Jets}", 10, -0.5, 9.5); 
  h_jet_mult_inc = dir.make<TH1F>("h_jet_mult_incl", "h_jet_mult_incl; N_{Jets}^{Inc}", 10, -0.5, 9.5); 
  h_jet1jet2_mass = dir.make<TH1F>("h_jet1jet2_mass", "h_jet1jet2_mass; M_{jet1, jet2}", 100, 0.0, 500.);
  h_deltaR_jet1jet2 = dir.make<TH1F>("h_deltaR_jet1jet2", "h_deltaR_jet1jet2; #Delta_{R}^{jet1, jet2}", 50, 0.0, 5.0);
  h_jet_HadV_pt = dir.make<TH1F>("h_jet_HadV_pt", "h_jet_HadV_pt; p_{T}^{V} (GeV)", 100, 0.0, 1000.);
  h_jet_HadV_eta = dir.make<TH1F>("h_jet_HadV_eta", "h_jet_HadV_eta; #eta_{V}", 40, -5.0, 5.);  
  h_jet_HadV_phi = dir.make<TH1F>("h_jet_HadV_phi", "h_jet_HadV_phi; #phi_{V}", 20, -4.0, 4.);


  //Jet merging histos
  h_m1_vs_m12 = dir.make<TH2F>("h_m1_vs_m12", "Mass Jet1 vs Mass Jet12;M_{j1};M_{j1,j2}", 100, 0., 200., 100, 0., 200.);
  h_bestmass = dir.make<TH1F>("h_bestmass", "Best V Mass;M_{V}^{Best}", 75, 0., 150.);//Cory: Change back
  h_jet1mass_jet2mass = dir.make<TH2F>("h_jet1mass_jet2mass", "h_jet1mass_jet2mass", 100, 0.0, 500.0, 100, 0.0, 500.0);
  h_jet1jet2_mass_Restricted = dir.make<TH1F>("h_jet1jet2_mass_Restricted", "h_jet1jet2_mass_Restricted", 100, 0.0, 500.0);
  h_HadVZmass_Cory = dir.make<TH1F>("h_HadVZmass_Cory", "h_HadVZmass_Cory", 100, 0.0, 2500.0);
  h_HadVZmass_Flavia = dir.make<TH1F>("h_HadVZmass_Flavia", "h_HadVZmass_Flavia", 100, 0.0, 2500.0);
  h_HadV_mass_Cory = dir.make<TH1F>("h_HadVmass_Cory", "h_HadVmass_Cory", 100, 0.0, 500.0);
  h_HadV_mass_Flavia = dir.make<TH1F>("h_HadVmass_Flavia", "h_HadVmass_Flavia", 100, 0.0, 500.0);
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
  defineHistoSet("hVZMass", "Reconstructed VZ Invariant Mass",
                  "M_{VZ} (GeV)", 100, 0, 2500, "GeV", hVZMass,dir);
  defineHistoSet("hVZeeMass", "Reconstructed VZee Invariant Mass",
                  "M_{VZ}^{ee} (GeV)", 100, 0, 2500, "GeV", hVZeeMass,dir);
  defineHistoSet("hVZmmMass", "Reconstructed VZmm Invariant Mass",
                  "M_{VZ}^{#mu#mu} (GeV)", 100, 0, 2500, "GeV", hVZmmMass,dir);
  defineHistoSet("hVZpt", "Reconstructed VZ Transverse Momentum",
                  "p_{VZ}^{T} (GeV)", 100, 0, 1000, "GeV", hVZpt,dir);

  defineHistoSet("hQ", "Q=M_{VZ} - M_{V} - M_{Z}",
                  "Q (GeV)", 50, 0, 2500, "GeV", hQ,dir);
  defineHistoSet("hMET", "Missing Transverse Energy", 
		 "MET (GeV)", 100, 0., 1000., "GeV", hMET,dir);
  defineHistoSet("hMETmm", "Missing Transverse Energy",
                 "MET (GeV)", 100, 0., 1000., "GeV", hMETmm,dir);
  defineHistoSet("hMETee", "Missing Transverse Energy",
                 "MET (GeV)", 100, 0., 1000., "GeV", hMETee,dir);


  defineHistoSet("hZMass" , "Reconstructed Mass of Z",
                  "M_{Z} (GeV)", 30, 60, 120, "GeV", hZMass,dir);
  defineHistoSet("hZeeMass" , "Reconstructed Mass of Zee",
                  "M_{Z}^{ee} (GeV)", 30, 60, 120, "GeV", hZeeMass,dir);
  defineHistoSet("hZmmMass" , "Reconstructed Mass of Zmm",
                  "M_{Z}^{#mu#mu} (GeV)", 30, 60, 120, "GeV", hZmmMass,dir);
  defineHistoSet("hZpt", "p_{T}^{Z}", 
                  "p_{T}^{Z} (GeV)", 100, 0, 1000, "GeV", hZpt,dir);
  defineHistoSet("hEvtType", "Event Type",
                  "N_{#mu}", 3, 0, 3, "NONE", hEvtType,dir);
  defineHistoSet("hDRee", "#Delta_{R}^{e_{1},e_{2}}",
                 "#Delta_{R}^{e_{1},e_{2}}", 50, 0, 5, "NONE", hDRee, dir);
  defineHistoSet("hDRmm", "#Delta_{R}^{#mu_{1},#mu_{2}}",
                 "#Delta_{R}^{#mu_{1},#mu_{2}}", 50, 0, 5, "NONE", hDRmm, dir);
  
  defineHistoSet("heeVMass" , "Reconstructed Mass of V (Zee channel)",
                  "M_{V} (GeV)", 75, 0, 150, "GeV", heeVMass,dir);//Cory: Change back
  defineHistoSet("heeVpt", "p_{T}^{V} (Zee channel)", 
                  "p_{T}^{V} (GeV)", 100, 0, 1000, "GeV", heeVpt,dir);


  defineHistoSet("hmmVMass" , "Reconstructed Mass of V (Zmm channel)",
		 "M_{V} (GeV)", 75, 0, 150, "GeV", hmmVMass,dir);//Cory: Change back                                                                        
  defineHistoSet("hmmVpt", "p_{T}^{V} (Zmm channel)",
		 "p_{T}^{V} (GeV)", 100, 0, 1000, "GeV", hmmVpt,dir);

  defineHistoSet("hNVtxs", "Number of Vertexs in Event",
                 "N_{Vtx}", 50, 0, 50, "NONE", hNVtxs,dir);

  defineHistoSet("hmmNVtxs", "Number of Vertexs in Event",
                 "N_{Vtx}", 50, 0, 50, "NONE", hmmNVtxs,dir);

  defineHistoSet("heeNVtxs", "Number of Vertexs in Event",
                 "N_{Vtx}", 50, 0, 50, "NONE", heeNVtxs,dir);

  defineHistoSet("hNLLeps", "Number of Loose Leptons in Event",
                  "N_{l}^{Loose}", 10, 0, 10, "NONE", hNLLeps,dir);
  defineHistoSet("hNLJets", "Number of Loose Jets in Event",
                  "N_{Jets}^{Loose}", 10, 0, 10, "NONE", hNLJets,dir);

  cout << "Histos declared" << endl;

  tVZCand = dir.make<TTree>("tVZCand", "Analysis Variables after VZCand");//Only 1 for now;
  tVZCand->Branch("Run", &runNumber_);
  tVZCand->Branch("Lumi", &lumiNumber_);
  tVZCand->Branch("Event", &evtNumber_);
  tVZCand->Branch("VZMass", &VZMass_);
  tVZCand->Branch("EvtType", &evtType_);
  tVZCand->Branch("ZMass", &ZMass_);
  tVZCand->Branch("VMass", &VMass_);
  tVZCand->Branch("Zpt", &Zpt_);
  tVZCand->Branch("Vpt", &Vpt_);
  tVZCand->Branch("Q", &Q_);
  tVZCand->Branch("DeltaRll", &DeltaRll_);
  tVZCand->Branch("DeltaRVZ", &DeltaRVZ_);
  tVZCand->Branch("weight", &weight_);


}//defineHistos


void 
HadronicVZAnalyzer::eventLoop(edm::EventBase const & event){
  clearEvtVariables();
  runNumber_ = event.id().run();
  lumiNumber_ = event.id().luminosityBlock();
  evtNumber_ = event.id().event();
  if(debug_) WPrimeUtil::printEvent(event);
  
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
  if(patElectronsH_->size() + patMuonsH_->size() == 0){
    cout << "Not enough leptons. Bad bad event, returning now..." << endl;
    return;
  }
  /*
  if(!wprimeUtil_->runningOnData())
    {

      gravMass_ = -999.0;

      // std::vector<reco::GenParticle> genMuons;                                                                                                               
      edm::Handle< vector<reco::GenJet  > > genJets;

      std::vector<reco::GenParticle> genZ;
      std::vector<reco::GenParticle> genW;

      int isW = 0;
      int isZ = 0;

      event.getByLabel(genLabel_, genParticles);
      for(size_t i = 0; i != genParticles->size(); ++i) {
	const reco::GenParticle & genP = (*genParticles)[i];
	if (genP.pdgId() == 5000039 && genP.status() == 3)
	  {
	    gravMass_ = genP.mass();
	    // cout << "Found my graviton! Mass is " << genP.mass() << endl;
	  }
	if (fabs(genP.pdgId()) == 13 && genP.status() == 3)
	  genMuons.push_back((*genParticles)[i]);

	if (fabs(genP.pdgId()) == 11 && genP.status() == 3)
          genElectrons.push_back((*genParticles)[i]);

	if (fabs(genP.pdgId()) == 24 && genP.status() == 3)
	  {
	    //cout << "Found W!" << endl;
	    h_genWMass->Fill(genP.mass(), weight_);
	    genW.push_back((*genParticles)[i]);
	    isW = 1;
	  }
	if (fabs(genP.pdgId()) == 23 && genP.status() == 3)
	  {
	    // cout << "Found Z! With mass of " << genP.mass() << endl;
	    h_genZMass->Fill(genP.mass(), weight_);
            genZ.push_back((*genParticles)[i]);
	    isZ = 1;
	  }


      } // loop over genParticles

      //Get the genJets
      event.getByLabel(edm::InputTag("ak7GenJets"), genJets);
      int jetID = 0;
      for (size_t i = 0; i < genJets->size(); i++)
	{
	  if ( ((*genJets)[i].pt() > 30.0) && (fabs((*genJets)[i].eta()) < 2.4) )
	    jetID = 1;

	  if (jetID==1 && !Overlap((*genJets)[i], genMuons, 1.0, 2) && !Overlap((*genJets)[i], genElectrons, 1.0, 2)) 
	    ak7GenJet.push_back((*genJets)[i]);
	}
      
      sort(ak7GenJet.begin(), ak7GenJet.end(), highestJetPt());

      //cout << "size of gen jet is " << ak7GenJet.size() << endl;

      reco::CompositeCandidate HadVZgen;
      reco::CompositeCandidate HadVWgen;
      
      //Make the gen combination of Z/W + extra jet only if extra jet is V-like
      if (ak7GenJet.size()>0)
	{

	  if (isW == 1){
	    if (ak7GenJet.at(0).mass()>70.0 && ak7GenJet.at(0).mass()<120.0 && ak7GenJet.at(0).pt()>290.0 && genW.at(0).pt()>150.0 && reco::deltaR(genW.at(0), ak7GenJet.at(0)) > 1.0 )
	      {
		HadVWgen.addDaughter(genW.at(0));
		HadVWgen.addDaughter(ak7GenJet.at(0));
		AddFourMomenta addFM;
		addFM.set(HadVWgen);
		h_HadVWgenMass->Fill(HadVWgen.mass()/1000, weight_);
	      }
	  }
	  
          if (isZ == 1){
            if (ak7GenJet.at(0).mass()>70.0 && ak7GenJet.at(0).mass()<120.0 && ak7GenJet.at(0).pt()>290.0 && genZ.at(0).pt()>150.0 && reco::deltaR(genZ.at(0), ak7GenJet.at(0)) > 1.0)
	      { 
		HadVZgen.addDaughter(genZ.at(0));
		HadVZgen.addDaughter(ak7GenJet.at(0));
		AddFourMomenta addFM;
		addFM.set(HadVZgen);
                h_HadVZgenMass->Fill(HadVZgen.mass()/1000, weight_);
	      }
	  }
	}//more than one jet if

    }//running on MC if
  
  */

  // for (size_t i=0; i<genMuons.size(); i++)
  //  cout << "My muon " << i << " pt is " << genMuons[i].pt() << endl;
  //  weight_ = wprimeUtil_->getWeight();

  // Make vectors of leptons passing various criteria
  // Loop over electrons, and see if they pass the criteria
  for (size_t i = 0; i < patElectronsH_->size(); i++) {
    allElectrons_.push_back(heep::Ele((*patElectronsH_)[i]));   
    
    if(Overlap(allElectrons_[i].patEle(), *patMuonsH_.product(), 0.01)) continue;
    
    if (looseElectron_(allElectrons_[i].patEle()))
      looseElectrons_.push_back(allElectrons_[i]);

    if (tightElectron_(allElectrons_[i].patEle()))
      tightElectrons_.push_back(allElectrons_[i]);
  }

  // Loop over muons, and see if they pass the TeVMuon criteria  
  for (size_t i = 0; i < patMuonsH_->size(); i++) {
    allMuons_.push_back(TeVMuon((*patMuonsH_)[i],muReconstructor_));   
    
    if (looseMuon_(allMuons_[i]) )
      looseMuons_.push_back(allMuons_[i]);
    
    if (tightMuon_(allMuons_[i]) )
      tightMuons_.push_back(allMuons_[i]);
  }

  if(debug_){
    printLeptons();
    printf("    Contains: %i electron(s), %i muon(s)\n",
           (int)allElectrons_.size(), (int)allMuons_.size());
    printf("    Contains: %i loose electron(s), %i loose muon(s)\n",
           (int)looseElectrons_.size(), (int)looseMuons_.size());
    printf("    Contains: %i tight electron(s), %i tightmuon(s)\n",
           (int)tightElectrons_.size(), (int)tightMuons_.size());
  }


  // if(useAdjustedMET_) event.getByLabel(pfCandsLabel_, pfCandidatesH_);
  //Deal with MET
  event.getByLabel(metLabel_, metH_);
  WPrimeUtil::getLeptonsMET(patElectronsH_, allElectrons_,
                            patMuonsH_, muReconstructor_, allMuons_,
                            metH_, useAdjustedMET_, met_,
                            pfCandidatesH_);

  // h_MET_AllCuts->Fill(met_.et(), weight_);

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

  fillPOGMuonHists();

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


  /*
  //Loop over jets to change their JES correction - systematics calculation
  
  JetCorrectionUncertainty *jecUnc_PF;
  jecUnc_PF =
    (JetCorrectionUncertainty*) new JetCorrectionUncertainty("GR_R_42_V19_AK7PF_Uncertainty.txt");
  

  
  for (size_t i=0; i<allJets_.size(); ++i) {

    float mass = allJets_[i].mass();
    float ptUnscaled = allJets_[i].pt();
    
    // estimate the uncertainty                                                                                                                             
    jecUnc_PF->setJetEta(allJets_[i].eta());
    jecUnc_PF->setJetPt(ptUnscaled);
    
    int scaleEnergy = 0;
    //If I want to scale up
    scaleEnergy = 1.0;
    //If I want to scale down
    //scaleEnergy = -1.0;
    
    // apply the uncertainty                                                                                                                                
    float pt = ptUnscaled + scaleEnergy*jecUnc_PF->getUncertainty(true)*ptUnscaled;
    float p = pt/fabs(sin(allJets_[i].theta()));
    float energy = sqrt(p*p+mass*mass);
    float eta = allJets_[i].eta();
    float phi =allJets_[i].phi();

    math::PtEtaPhiMLorentzVector polar = allJets_[i].polarP4();

    polar.SetPt(pt);
    allJets_[i].setP4(polar);


    // allJets_[i].setP4(LorentzVector(pt, eta, phi, energy));


    }
  */
  // Loop over jets, and see if they pass the jet criteria
  for (size_t i = 0; i < allJets_.size(); ++i) {
    if (looseJet_(allJets_[i]) && !Overlap(allJets_[i], looseMuons_, 1.0, 2) && !Overlap(allJets_[i], looseElectrons_, 1.0, 2))
      looseJets_.push_back(allJets_[i]);
  }
  
  /*  if (looseJets_.size() > 2){
    if (debug_)
      cout << "Too many jets. Returning now" << endl;
    return; 
  }
  */
  int cont40 = 0;
  int cont50 = 0;
  int cont60 = 0;
  int cont70 = 0;

  for (size_t i = 0; i < looseJets_.size(); i++)
    {

      if (looseJets_.at(i).pt()>40.0)
	cont40 = cont40+1;

      if (looseJets_.at(i).pt()>50.0)
	cont50 = cont50+1;

      if (looseJets_.at(i).pt()>60.0)
	cont60 = cont60+1;

      if (looseJets_.at(i).pt()>70.0)
	cont70 = cont70+1;
    }

  h_ptJetCut_nJets->Fill(looseJets_.size(), 30, weight_);
  h_ptJetCut_nJets->Fill(cont40, 40, weight_);
  h_ptJetCut_nJets->Fill(cont50, 50, weight_);
  h_ptJetCut_nJets->Fill(cont60, 60, weight_);
  h_ptJetCut_nJets->Fill(cont70, 70, weight_);


  fillJetMultiplicityHists();
  


  //fill jet histos for jets who passes the criteria
  if (looseJets_.size() > 0){
    sort(looseJets_.begin(), looseJets_.end(), highestJetPt());
    if(debug_){
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
    
    /* if (looseJets_.at(0).mass() > 60 && looseJets_.at(0).mass() < 110){
      h_jet1_Vwindow->Fill(looseJets_.at(0).mass(), weight_);
      }*/

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

  if (debug_){
    printf("    Contains: %i loose jets(s)\n",
           (int)looseJets_.size());
  }

  //get Vertex
  event.getByLabel(vertexLabel_, verticesH_);


  //////////////////////////////////////////////
  /// Start Applying Cuts///////////////////////
  //////////////////////////////////////////////

  //Make a WCand out of looseLeptons
  //WCandidate wCands = getWCand
  WCandidate wCand_ = getWCand(looseElectrons_, looseMuons_, met_);
  h_WMass->Fill(wCand_.mt(), weight_);
  int goodW = 0;
  if (wCand_.mt() > 70 && wCand_.mt() < 110 && wCand_.pt() > 150)
    {
      goodW = 1;
    }
  
  if (looseJets_.size()>0)
    {

      vCand_ = getVCand(looseJets_);
      if (goodW==1 && vCand_.mass() > 70 && vCand_.mass() < 120 && vCand_.pt() > 100)
	// if (wCand_.mt() > 70 && wCand_.mt() < 110 && wCand_.pt() > 150 && vCand_.mass() > 70 && vCand_.mass() < 120 && vCand_.pt() > 100)
	{
	  wzCand_ = (vCand_ && wCand_) ? XWLeptonic(vCand_, wCand_) : XWLeptonic();
	  double WZMass = wzCand_().mass();
	  if (WZMass > 0){
	    h_HadVWMass->Fill(WZMass/1000, weight_);
	    h_VWMass_nJets->Fill(looseJets_.size(), WZMass);
	    /* cout << "W lepton four vector " << endl;
	    cout << "px =  " << wCand_.daughter(0)->px() << endl; 
	    cout << "py =  " << wCand_.daughter(0)->py() << endl;
            cout << "pz =  " << wCand_.daughter(0)->pz() << endl;
            cout << "e =  " << wCand_.daughter(0)->energy() << endl;


	    cout << "MET Info" << endl;
	    cout << "MET px = " << wCand_.daughter(1)->px() << endl; 
            cout << "MET py = " << wCand_.daughter(1)->py() << endl;
            cout << "MET pz = " << wCand_.daughter(1)->pz() << endl;


	    cout << "Jet (vCand) Info" << endl;
	    cout << "Jet px = " << vCand_.px() << endl;
            cout << "Jet py = "<< vCand_.py() << endl;
            cout << "Jet pz = "<< vCand_.pz() << endl;
            cout << "Jet e = "<< vCand_.energy() << endl;
	    
	    cout << "W Transverse Mass =  " << wCand_.mt() << endl;
	    cout << "Neutrino Min pz = " << wzCand_.neutrinoPz(kMinPz) << endl;
	    cout << "Neutrino Max pz = " << wzCand_.neutrinoPz(kMaxPz) << endl;

	    cout << "WZ Mass" << endl;
	    cout << "WZMass(MinPZ) = " << wzCand_(kMinPz).mass() << endl; 
            cout << "WZMass(MaxPZ) = " << wzCand_(kMaxPz).mass() << endl;
	    */

	    //h_MET_AllCuts->Fill(met_.et(), weight_);
	  
	  }

	}
    }

  vCand_ = ZCandidate();
  wCand_ = WCandidate();
  hadVZ_ = VZCandidate();
  int iCut=0;
  if( !passNoCut() ) return;
  tabulateEvent(iCut, weight_); ++iCut;

  if( !passMinNLeptonsCut(looseElectrons_, looseMuons_, minNLeptons_) ) return;
  tabulateEvent(iCut, weight_); ++iCut;

  if( !passMinNJetsCut(looseJets_, minNJets_) ) return;
  tabulateEvent(iCut, weight_); ++iCut;

  //////////////////////////////
  ///////  Z With Muons  ///////
  //////////////////////////////
  
  // Make a Z candidate out of the loose leptons. 
  ZCandV zCands = getZCands(looseElectrons_, looseMuons_, 100., false);
  removeLowLepPtCands(zCands, 45., 10);
  sort(zCands.begin(), zCands.end(), closestToZMass());
  removeOverlapping(zCands);

  zCand_ = zCands.size() ? zCands[0] : ZCandidate();
  if (debug_) cout << "Made zCand" << endl;

  if( !passValidZCut(zCand_) ) return;
  tabulateEvent(iCut, weight_); ++iCut;

  if( !passZMassCut(zCand_, minZmass_, maxZmass_) ) return;
  tabulateEvent(iCut, weight_); ++iCut;

  if( !passZptCut  (zCand_, minZpt_) ) return;
  tabulateEvent(iCut, weight_); ++iCut;
  
  fillGoodZHistos();
  
  ///////////////////////////////////////
  //////// Make V from Jets  ////////////
  ///////////////////////////////////////
  
  fillJetMergingHistos();

  vCand_ = getVCand(looseJets_);
  if (debug_) cout << "Made vCand" << endl; 

  if( !passValidZCut(vCand_) ) return;
  tabulateEvent(iCut, weight_); ++iCut;

  LorentzVector j1j2 = looseJets_.at(0).p4();
  if(looseJets_.size()>1) j1j2 += looseJets_.at(1).p4();
  h_m1_vs_m12->Fill(looseJets_.at(0).mass(), j1j2.mass(), weight_);
  float bestMass = fabs(j1j2.mass() - 85.) < fabs(looseJets_.at(0).mass() - 85.) ? j1j2.mass() : looseJets_.at(0).mass();
  h_bestmass->Fill(bestMass, weight_);
  //if(bestMass < 60) cout<<" nJets: "<<looseJets_.size()<<" best: "<<bestMass<<" m1: "<<looseJets_.at(0).mass()<<" m12: "<<j1j2.mass()<<endl;

  if( !passZMassCut(vCand_, minVmass_, maxVmass_) ) return;
  tabulateEvent(iCut, weight_); ++iCut;

  /*if( !passZptCut(vCand_, minVpt_) ) return;
  tabulateEvent(iCut, weight_); ++iCut;

  if (debug_) cout << "Good V from jet" << endl;
  fillGoodHadVHistos();
  */
  ///////////////////////////////////////
  //////// Make VZ Candidate ////////////
  ///////////////////////////////////////
  if (debug_) cout << "Im just before boson cands" << endl;
  
  hadVZ_ = VZCandidate(zCand_, vCand_);
  if(debug_) cout << "Made my hadVZ" << endl;  

  if( !passValidVZCandCut() ) return;
  tabulateEvent(iCut, weight_); ++iCut;
  //h_MET_AllCuts->Fill(met_.et(), weight_);

  fillValidVZHistos();

  //get Trigger 
  //triggerEvent_ = getProduct<pat::TriggerEvent>(event,hltEventLabel_); 
  event.getByLabel(hltEventLabel_, triggerEventH_);

  if (!passTriggersCut())
    return;
  tabulateEvent(iCut, weight_); ++iCut;

  /*if( !passZptCut  (zCand_, minZpt_) ) return;
  tabulateEvent(iCut, weight_); ++iCut;

  fillGoodZHistos();
  */
  if( !passZptCut(vCand_, minVpt_) ) return;
  tabulateEvent(iCut, weight_); ++iCut;

  if (debug_) cout << "Good V from jet" << endl;
  fillGoodHadVHistos();


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
HadronicVZAnalyzer::passValidVZCandCut(){
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
void HadronicVZAnalyzer::printEventDetails() const{
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
HadronicVZAnalyzer::clearEvtVariables(){
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

}

//fill Histograms
void HadronicVZAnalyzer::fillHistos(const int& index, const float& weight){
  if(debug_) printf("filling Histos\n");

  if(hadVZ_){
    hVZMass[index]->Fill(hadVZ_.mass(), weight);
    if      (zCand_.flavor() == PDG_ID_ELEC) hVZeeMass[index]->Fill(hadVZ_.mass(), weight);
    else if (zCand_.flavor() == PDG_ID_MUON) hVZmmMass[index]->Fill(hadVZ_.mass(), weight);
    hVZpt[index]->Fill(hadVZ_.pt(), weight);
  }
  if(zCand_){
    hZMass[index]->Fill(zCand_.mass(), weight);
    if      (zCand_.flavor() == PDG_ID_ELEC) 
      {
	hZeeMass[index]->Fill(zCand_.mass(), weight);
	hMETee[index]->Fill(met_.et(), weight);
  hDRee[index]->Fill(deltaR(*zCand_.daughter(0), *zCand_.daughter(1)), weight);
      }
    else if (zCand_.flavor() == PDG_ID_MUON) 
      {
	hZmmMass[index]->Fill(zCand_.mass(), weight);
	hMETmm[index]->Fill(met_.et(), weight);
  hDRmm[index]->Fill(deltaR(*zCand_.daughter(0), *zCand_.daughter(1)), weight);
      }
    hZpt[index]->Fill(zCand_.pt(), weight);
    hEvtType[index]->Fill(2*(zCand_.flavor() == PDG_ID_MUON), weight);
  }
  if(vCand_){
   
    if (zCand_.flavor() == PDG_ID_MUON)
      {
	hmmVMass[index]->Fill(vCand_.mass(), weight);
	hmmVpt[index]->Fill(vCand_.pt(), weight);
      }
    else if (zCand_.flavor() == PDG_ID_ELEC)
      { 
        heeVMass[index]->Fill(vCand_.mass(), weight);
        heeVpt[index]->Fill(vCand_.pt(), weight);
      }


  }

  hNLLeps[index]->Fill(looseElectrons_.size()+looseMuons_.size(), weight);
  hNLJets[index]->Fill(looseJets_.size(), weight);

  hNVtxs[index]->Fill((*verticesH_).size(), weight);

  if(zCand_){
    if      (zCand_.flavor() == PDG_ID_ELEC)
      {
	heeNVtxs[index]->Fill((*verticesH_).size(), weight);
      }
    else if (zCand_.flavor() == PDG_ID_MUON)
      {
	hmmNVtxs[index]->Fill((*verticesH_).size(), weight);
      }
  }

  hMET[index]->Fill(met_.et(), weight);
  
  if(CutNames_[index] == "ValidVZ"){
    if(hadVZ_){
      VZMass_ = hadVZ_.mass();
      DeltaRVZ_ = deltaR(zCand_, vCand_);
      Q_ = hadVZ_.mass() - zCand_.mass() - vCand_.mass();
    }

    if(zCand_){
      evtType_ = 2*(zCand_.flavor() == PDG_ID_MUON);
      ZMass_ = zCand_.mass();
      Zpt_ = zCand_.pt();
      DeltaRll_ = deltaR(*zCand_.daughter(0), *zCand_.daughter(1));
    }
    if(vCand_){
      VMass_ = vCand_.mass();
      Vpt_ = vCand_.pt();
      weight_ = weight;
    }
    tVZCand->Fill();
  }
  
}//fillHistos

void HadronicVZAnalyzer::fillJetMultiplicityHists(){
  h_jet_mult->Fill(looseJets_.size(), weight_);
  uint max = min(10, (int)looseJets_.size());
  for(uint i=0; i<=max; ++i)
    h_jet_mult_inc->Fill(i, weight_);
}


void HadronicVZAnalyzer::fillGoodZHistos(){
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

void HadronicVZAnalyzer::fillGoodHadVHistos(){
  h_jet_HadV_pt->Fill(vCand_.pt(), weight_);
  h_jet_HadV_eta->Fill(vCand_.eta(), weight_);
  h_jet_HadV_phi->Fill(vCand_.phi(), weight_);
  if (debug_)
    cout << "filled my HadV histos" << endl;
}

void HadronicVZAnalyzer::fillValidVZHistos(){
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

void HadronicVZAnalyzer::fillJetMergingHistos(){

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



void HadronicVZAnalyzer::fillPOGMuonHists(){
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
