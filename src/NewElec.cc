#include "UserCode/CMGWPrimeGroup/interface/NewElec.h"

float
CalcElecCombRelIso(const pat::Electron & elec, const float pu){
  float num = elec.dr03TkSumPt();
  num += (elec.dr03HcalTowerSumEt() + 
          elec.hadronicOverEm() * elec.energy() * fabs(sin(elec.superCluster()->position().theta())));
  num -= pu;
  num += elec.isEB() ? std::max(0., elec.dr03EcalRecHitSumEt() - 1.) : elec.dr03EcalRecHitSumEt();
  return num / elec.p4().Pt();
}
