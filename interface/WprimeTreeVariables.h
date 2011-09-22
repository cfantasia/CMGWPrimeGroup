#include "TFile.h"
#include "TTree.h"
#include "Math/Vector4D.h"
#include "TLorentzVector.h"


struct WprimeVariables
{
  // tree definition
  TFile* m_outputRootFile;
  TTree* m_reducedTree;
  
  
  // input parameters
  float mW;
  int totEvents;
  float crossSection;
  int dataFlag; 
  int runId; 
  int lumiId; 
  int eventId; 
  int hltPrescale;
  
  TLorentzVector ele;
  TLorentzVector* p_ele;
  
  float ele_eSeed;
  float ele_timeSeed;
  float ele_eSC;
  float ele_e1x5;
  float ele_e2x5;
  float ele_e5x5;
  float ele_charge;
  float ele_dxy;
  float ele_dz;
  float ele_tkIso;
  float ele_emIso;
  float ele_hadIso;
  float ele_hadIso_d1;
  float ele_hadIso_d2;
  int ele_isEB;
  int ele_isEE;
  int ele_isEcalDriven;
  float ele_sigmaIetaIeta;
  float ele_DphiIn;
  float ele_DetaIn;
  float ele_HOverE;
  float ele_fbrem;
  float ele_EOverP;

  // met variables
  TLorentzVector met;
  TLorentzVector* p_met;
  
  TLorentzVector eleMet;
  float eleMet_mt;
  float eleMet_Dphi;
  float phoMet_mt;
  float phoMet_Dphi;
  
  
};



void InitializeTree(WprimeVariables& vars, TTree* tree, std::string& analysis);
void SetBranchAddresses(WprimeVariables& vars, TTree* tree);

void ClearWprimeVariables(WprimeVariables&, std::string& analysis);
void DeleteWprimeVariables(WprimeVariables&);
TTree* CloneTree(WprimeVariables& vars);

