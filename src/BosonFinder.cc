#include "UserCode/CMGWPrimeGroup/interface/BosonFinder.h"

using std::vector; using std::cout; using std::endl;

//////// Useful Algorithms ///////////////////////////////////////////////////
/// Return a vector of non-overlapping Z candidates

ZCandV getZCands(const ElectronV & electrons, const MuonV & muons, 
                 float maxMassDiff, bool rmOverlap,  
                 const vector<bool> & eMask, const vector<bool> & mMask)
{

  ZCandV zCands;

  ZCandV zeeCands = getZCands(electrons, maxMassDiff, rmOverlap, eMask);
  zCands.insert(zCands.end(), zeeCands.begin(), zeeCands.end());

  ZCandV zmmCands = getZCands(    muons, maxMassDiff, rmOverlap, mMask);
  zCands.insert(zCands.end(), zmmCands.begin(), zmmCands.end());

  // Order by difference from Z mass
  sort(zCands.begin(), zCands.end(), closestToZMass());

  return zCands;
}

void removeOverlapping(ZCandV & zCands){
  // get rid of overlapping candidates
  for (ZCandV::iterator i = zCands.begin(); i != zCands.end(); ++i)
    for (ZCandV::iterator j = i + 1; j != zCands.end(); ++j)
      if (areOverlapping(*i, *j)) {
        zCands.erase(j);
        j = i;
      }
}

void removeWorstCands(ZCandV & zCands, const float & maxMassDiff){
  for (ZCandV::iterator i = zCands.begin(); i != zCands.end(); ++i)
    if( fabs(i->mass() - ZMASS) > maxMassDiff){
      zCands.erase(i);
      i--;
    }
}

void removeWorstCands(ZCandV & zCands, const float& min, const float& max){
  for (ZCandV::iterator i = zCands.begin(); i != zCands.end(); ++i)
    if(i->mass() > max  || i->mass() < min){
      zCands.erase(i);
      i--;
    }
}

void removeLowLepPtCands(ZCandV & zCands, const float& minPt1, const float& minPt2){
  for (ZCandV::iterator i = zCands.begin(); i != zCands.end(); ++i){
    float pt1 = std::max(i->daughter(0)->pt(), i->daughter(1)->pt());
    float pt2 = std::min(i->daughter(0)->pt(), i->daughter(1)->pt());
    if(pt1 < minPt1 || pt2 < minPt2){
      zCands.erase(i);
      i--;
    }
  }
}

WCandV getWCandidates(const ElectronV & electrons,
                      const MuonV & muons, const pat::MET & met){

  vector<WCandidate> eWcands = getWCandidates(electrons, met);
  vector<WCandidate> muWcands = getWCandidates(muons, met);
  vector<WCandidate>::const_iterator it;
  for (it = muWcands.begin(); it != muWcands.end(); ++it)
    eWcands.push_back(*it);
  return eWcands;
}

WCandidate getWCand(const ElectronV & electrons,
                    const MuonV & muons, 
                    const pat::MET & met){
  vector<WCandidate> wCands = getWCandidates(electrons, muons, met);
  
  sort(wCands.begin(), wCands.end(), highestPtLepton());
  if (wCands.size()) return wCands[0];
  
 return WCandidate();
}

/// Return a hadronic WCandidate
ZCandidate getVCand(const JetV & jets, const int mode)
{
  size_t bestPos = 0;
  if(mode==0){
    // Order by pt - do it the dumb way because we only have a const&
    double maxPt = -1.0;
    for (size_t i=0; i!= jets.size(); ++i) {
      if(jets[i].pt() > maxPt) {
        maxPt = jets[i].pt();
        bestPos = i;
      }
    }
  }else if(mode==1){
    //order by closest to 85GeV (VMASS)
    double minMassDiff = 999.0;
    for (size_t i=0; i!=jets.size(); ++i) {
      if(fabs(jets[i].mass()-VMASS) < minMassDiff ) {
        minMassDiff = fabs(jets[i].mass()-VMASS);
        bestPos = i;
      }
    }
  }
  return jets.size() ? ZCandidate(jets[bestPos]) : ZCandidate();
}

ZCandidate getVCand2(const JetV & jets)
{
  // Order by pt - do it the dumb way because we only have a const&
  if(jets.size() < 2) return ZCandidate();
  double maxPt = -1.0;
  size_t maxPtPosition1 = 0;
  size_t maxPtPosition2 = 1;
  for (size_t i=0; i!= jets.size(); ++i) {
    if(jets[i].pt() > maxPt) {
      maxPt = jets[i].pt();
      maxPtPosition1 = i;
    }
  }
  maxPt = -1.0;
  for (size_t i=0; i!= jets.size(); ++i) {
    if(i == maxPtPosition1) continue;
    if(jets[i].pt() > maxPt) {
      maxPt = jets[i].pt();
      maxPtPosition2 = i;
    }
  }
  
  return ZCandidate(jets[maxPtPosition1], jets[maxPtPosition2]);
}

/// Return a WCandidate using the highest-pT non-Z lepton
WCandidate getWCand(const ElectronV & electrons,
                    const MuonV & muons, 
                    const pat::MET & met,
                    const ZCandidate & zCand,
                    double minDeltaR)
{
  vector<WCandidate> wCands;
  
  for (ElectronV::const_iterator i = electrons.begin(); 
       i != electrons.end(); ++i)
    if (!zCand || 
        (!areOverlapping(i->patEle(), * zCand.daughter(0)) &&
         !areOverlapping(i->patEle(), * zCand.daughter(1)) &&
         reco::deltaR(* i, * zCand.daughter(0)) > minDeltaR &&
         reco::deltaR(* i, * zCand.daughter(1)) > minDeltaR)) {
      wCands.push_back(WCandidate(* i, met));
    }
  
  for (MuonV::const_iterator i = muons.begin(); 
       i != muons.end(); ++i)
    if (!zCand || 
        (!areOverlapping(* i, * zCand.daughter(0)) &&
         !areOverlapping(* i, * zCand.daughter(1)) &&
         reco::deltaR(* i, * zCand.daughter(0)) > minDeltaR &&
         reco::deltaR(* i, * zCand.daughter(1)) > minDeltaR)) {
      wCands.push_back(WCandidate(* i, met));
    }

  sort(wCands.begin(), wCands.end(), highestPtLepton());

  if (wCands.size()) return wCands[0];

  return WCandidate();
}

/// Return candidates using the highest-pT non-Z lepton and each of the METs
WCandV getWCands(const ElectronV & electrons, 
                 const MuonV & muons, 
                 const METV & mets, 
                 const ZCandidate & zCand,
                 double minDeltaR = 0.)
{
  WCandV wCands;
  for (size_t i = 0; i < mets.size(); i++) {
    WCandidate wCand = getWCand(electrons, muons, mets[i], zCand, minDeltaR);
    if (wCand) wCands.push_back(wCand);
  }
  return wCands;
}

 

