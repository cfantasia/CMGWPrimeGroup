#include "UserCode/CMGWPrimeGroup/interface/WprimeTreeVariables.h"



void InitializeTree(WprimeVariables& vars, TTree* tree, std::string& analysis)
{
  //-------------
  // Reduced tree
  //-------------
  
  vars.m_reducedTree = tree;
  
  vars.m_reducedTree -> Branch("mW",           &vars.mW,                     "mW/F");
  vars.m_reducedTree -> Branch("totEvents",    &vars.totEvents,       "totEvents/I");
  vars.m_reducedTree -> Branch("crossSection", &vars.crossSection, "crossSection/F");
  vars.m_reducedTree -> Branch("dataFlag",     &vars.dataFlag,         "dataFlag/I");
  vars.m_reducedTree -> Branch("runId",        &vars.runId,               "runId/I");
  vars.m_reducedTree -> Branch("lumiId",       &vars.lumiId,             "lumiId/I");
  vars.m_reducedTree -> Branch("eventId",      &vars.eventId,           "eventId/I");
  vars.m_reducedTree -> Branch("hltPrescale",  &vars.hltPrescale,       "hltPrescale/I");
  

  if(analysis == "eleMET")
    {
      // electron variables
      vars.m_reducedTree -> Branch("ele", "TLorentzVector", &vars.p_ele);
      vars.m_reducedTree -> Branch("ele_eSC",  &vars.ele_eSC,   "ele_eSC/F");
      vars.m_reducedTree -> Branch("ele_eSeed",  &vars.ele_eSeed,   "ele_eSeed/F");
      vars.m_reducedTree -> Branch("ele_timeSeed",  &vars.ele_timeSeed,   "ele_timeSeed/F");
      vars.m_reducedTree -> Branch("ele_e1x5",  &vars.ele_e1x5,   "ele_ex5/F");
      vars.m_reducedTree -> Branch("ele_e2x5",  &vars.ele_e2x5,   "ele_e2x5/F");
      vars.m_reducedTree -> Branch("ele_e5x5",  &vars.ele_e5x5,   "ele_e5x5/F");
      vars.m_reducedTree -> Branch("ele_charge",  &vars.ele_charge,   "ele_charge/F");
      vars.m_reducedTree -> Branch("ele_dxy",     &vars.ele_dxy,         "ele_dxy/F");
      vars.m_reducedTree -> Branch("ele_dz",      &vars.ele_dz,           "ele_dz/F");
      vars.m_reducedTree -> Branch("ele_tkIso",   &vars.ele_tkIso,     "ele_tkIso/F");
      vars.m_reducedTree -> Branch("ele_emIso",   &vars.ele_emIso,     "ele_emIso/F");
      vars.m_reducedTree -> Branch("ele_hadIso",  &vars.ele_hadIso,   "ele_hadIso/F");
      vars.m_reducedTree -> Branch("ele_hadIso_d1",  &vars.ele_hadIso_d1,   "ele_hadIso_d1/F");
      vars.m_reducedTree -> Branch("ele_hadIso_d2",  &vars.ele_hadIso_d2,   "ele_hadIso_d2/F");
      vars.m_reducedTree -> Branch("ele_isEB",          &vars.ele_isEB,                   "ele_isEB/I");
      vars.m_reducedTree -> Branch("ele_isEE",          &vars.ele_isEE,                   "ele_isEE/I");
      vars.m_reducedTree -> Branch("ele_isEcalDriven",          &vars.ele_isEcalDriven,       "ele_isEcalDriven/I");
      vars.m_reducedTree -> Branch("ele_sigmaIetaIeta", &vars.ele_sigmaIetaIeta, "ele_sigmaIetaIeta/F");
      vars.m_reducedTree -> Branch("ele_DphiIn",        &vars.ele_DphiIn,               "ele_DphiIn/F");
      vars.m_reducedTree -> Branch("ele_DetaIn",        &vars.ele_DetaIn,               "ele_DetaIn/F");
      vars.m_reducedTree -> Branch("ele_HOverE",        &vars.ele_HOverE,               "ele_HOverE/F");
      vars.m_reducedTree -> Branch("ele_fbrem",         &vars.ele_fbrem,                 "ele_fbrem/F");
      vars.m_reducedTree -> Branch("ele_EOverP",        &vars.ele_EOverP,               "ele_EOverP/F");
      
      // met variables
      vars.m_reducedTree -> Branch("met", "TLorentzVector", &vars.p_met);
      vars.m_reducedTree -> Branch("eleMet_mt",    &vars.eleMet_mt,       "eleMet_mt/F");
      vars.m_reducedTree -> Branch("eleMet_Dphi",  &vars.eleMet_Dphi,   "eleMet_Dphi/F");  
    }
}

//in case you want to read the ttree after
void SetBranchAddresses(WprimeVariables& vars, TTree* tree)
{
  tree -> SetBranchAddress("runId",        &vars.runId);
  tree -> SetBranchAddress("lumiId",       &vars.lumiId);
  tree -> SetBranchAddress("eventId",      &vars.eventId);
  tree -> SetBranchAddress("eleMet_mt",    &vars.eleMet_mt);

  vars.p_ele = &vars.ele;
  tree -> SetBranchAddress("ele",          &vars.p_ele);
  tree -> SetBranchAddress("ele_eSC",      &vars.ele_eSC);
  tree -> SetBranchAddress("ele_eSeed",    &vars.ele_eSeed);
  tree -> SetBranchAddress("ele_timeSeed", &vars.ele_timeSeed);
  tree -> SetBranchAddress("ele_e1x5",     &vars.ele_e1x5);
  tree -> SetBranchAddress("ele_e2x5",     &vars.ele_e2x5);
  tree -> SetBranchAddress("ele_e5x5",     &vars.ele_e5x5);
  tree -> SetBranchAddress("ele_charge",   &vars.ele_charge);
  tree -> SetBranchAddress("ele_dxy",      &vars.ele_dxy);
  tree -> SetBranchAddress("ele_dz",       &vars.ele_dz);
  tree -> SetBranchAddress("ele_tkIso",    &vars.ele_tkIso);
  tree -> SetBranchAddress("ele_emIso",    &vars.ele_emIso);
  tree -> SetBranchAddress("ele_hadIso",   &vars.ele_hadIso);
  tree -> SetBranchAddress("ele_hadIso_d1",  &vars.ele_hadIso_d1);
  tree -> SetBranchAddress("ele_hadIso_d2",  &vars.ele_hadIso_d2);
  tree -> SetBranchAddress("ele_isEB",          &vars.ele_isEB);
  tree -> SetBranchAddress("ele_isEE",          &vars.ele_isEE);
  tree -> SetBranchAddress("ele_isEcalDriven",  &vars.ele_isEcalDriven);
  tree -> SetBranchAddress("ele_sigmaIetaIeta", &vars.ele_sigmaIetaIeta);
  tree -> SetBranchAddress("ele_DphiIn",        &vars.ele_DphiIn);
  tree -> SetBranchAddress("ele_DetaIn",        &vars.ele_DetaIn);
  tree -> SetBranchAddress("ele_HOverE",        &vars.ele_HOverE);
  tree -> SetBranchAddress("ele_fbrem",         &vars.ele_fbrem);
  tree -> SetBranchAddress("ele_EOverP",        &vars.ele_EOverP);
  
  // met variables
  vars.p_met = &vars.met;
  tree -> SetBranchAddress("met",          &vars.p_met);
  tree -> SetBranchAddress("eleMet_mt",    &vars.eleMet_mt);
  tree -> SetBranchAddress("eleMet_Dphi",  &vars.eleMet_Dphi);  
}


void ClearWprimeVariables(WprimeVariables& vars, std::string& analysis)
{
  //run variables
  vars.runId = -1;
  vars.lumiId = -1;
  vars.eventId = -1;
  vars.hltPrescale = -1;

  if(analysis == "eleMET")
    {
      vars.ele = TLorentzVector(0., 0., 0., 0.);
      vars.p_ele = NULL;
      
      vars.ele_eSC = -1.;
      vars.ele_eSeed = -1.;
      vars.ele_timeSeed = -99;
      vars.ele_e1x5 = -1.;
      vars.ele_e2x5 = -1.;
      vars.ele_e5x5 = -1.;
      vars.ele_charge = -1.;
      vars.ele_dxy = -1.;
      vars.ele_dz = -99.;
      vars.ele_tkIso = -1.;
      vars.ele_emIso = -1.;
      vars.ele_hadIso = -1.;
      vars.ele_hadIso_d1 = -1.;
      vars.ele_hadIso_d2 = -1.;
      vars.ele_isEB = -1;
      vars.ele_isEcalDriven = -1;
      vars.ele_sigmaIetaIeta = -1.;
      vars.ele_DphiIn = -99.;
      vars.ele_DetaIn = -99.;
      vars.ele_HOverE = -1.;
      vars.ele_fbrem = -1.;
      vars.ele_EOverP = -1.;
      
      // met variables 
      vars.met = TLorentzVector(0., 0., 0., 0.);
      vars.p_met = NULL;
      
      vars.eleMet = TLorentzVector(0., 0., 0., 0.);
      
      vars.eleMet_mt = -1.;
      vars.eleMet_Dphi = -1.;
    }
 
}



void DeleteWprimeVariables(WprimeVariables& vars)
{
  // save tree
  //not needed if we use the TFileService
}

TTree* CloneTree(WprimeVariables& vars)
{
  //clone tree
  return vars.m_reducedTree -> CloneTree(0);
}

