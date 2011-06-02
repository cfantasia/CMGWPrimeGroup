#include "UserCode/CMGWPrimeGroup/interface/BosonFinder.h"

//////// Useful Algorithms ///////////////////////////////////////////////////
/// Return a vector of non-overlapping Z candidates
ZCandV getZCands(const ElectronV & electrons, float maxMassDiff)
{
  ZCandV zCands;
  // Get opposite-charge pairs of electrons
  for (size_t i = 0; i < electrons.size(); i++)
    for (size_t j = i + 1; j < electrons.size(); j++)
      if (electrons[i].charge() != electrons[j].charge())
        zCands.push_back(ZCandidate(electrons[i], electrons[j]));

  // Order by difference from Z mass
  sort(zCands.begin(), zCands.end(), closestToZMass());

  removeOverlapping(zCands);

  return zCands;
}

ZCandV getZCands(const MuonV & muons, float maxMassDiff)
{
  ZCandV zCands;
  // Get opposite-charge pairs of muons
  for (size_t i = 0; i < muons.size(); i++)
    for (size_t j = i + 1; j < muons.size(); j++)
      if (muons[i].charge() != muons[j].charge())
        zCands.push_back(ZCandidate(muons[i], muons[j]));
  
  // Order by difference from Z mass
  sort(zCands.begin(), zCands.end(), closestToZMass());
  
  removeOverlapping(zCands);

  return zCands;
}

ZCandV getZCands(const ElectronV & electrons, 
                 const MuonV & muons, float maxMassDiff)
{

  ZCandV zCands;

  ZCandV zeeCands = getZCands(electrons, maxMassDiff);
  zCands.insert(zCands.end(), zeeCands.begin(), zeeCands.end());

  ZCandV zmmCands = getZCands(    muons, maxMassDiff);
  zCands.insert(zCands.end(), zmmCands.begin(), zmmCands.end());

  // Order by difference from Z mass
  sort(zCands.begin(), zCands.end(), closestToZMass());

  removeOverlapping(zCands);

  return zCands;
}

void removeOverlapping(ZCandV & zCands){
  // Get rid of overlapping candidates
  for (ZCandV::iterator i = zCands.begin(); i != zCands.end(); ++i)
    for (ZCandV::iterator j = i + 1; j != zCands.end(); ++j)
      if (areOverlapping(*i, *j)) {
        zCands.erase(j);
        j = i;
      }
}

WCandidate getWCand(const ElectronV & electrons,
                    const pat::MET & met){
  vector<WCandidate> wCands;
  
  for (ElectronV::const_iterator i = electrons.begin(); 
       i != electrons.end(); ++i)
    wCands.push_back(WCandidate(* i, met));

  sort(wCands.begin(), wCands.end(), highestPtLepton());

  if (wCands.size()) return wCands[0];

  return WCandidate();

}

WCandidate getWCand(const MuonV & muons, 
                    const pat::MET & met){
  vector<WCandidate> wCands;
  for (MuonV::const_iterator i = muons.begin(); 
       i != muons.end(); ++i)
    wCands.push_back(WCandidate(* i, met));

  sort(wCands.begin(), wCands.end(), highestPtLepton());

  if (wCands.size()) return wCands[0];

  return WCandidate();
}

WCandidate getWCand(const ElectronV & electrons,
                    const MuonV & muons, 
                    const pat::MET & met){
 vector<WCandidate> wCands;

 wCands.push_back(getWCand(electrons, met));
 wCands.push_back(getWCand(muons    , met));

 sort(wCands.begin(), wCands.end(), highestPtLepton());
 
 if (wCands.size()) return wCands[0];
 
 return WCandidate();
}

/// Return a WCandidate using the highest-pT non-Z lepton
WCandidate getWCand(const ElectronV & electrons,
                    const MuonV & muons, 
                    const pat::MET & met,
                    const ZCandidate & zCand,
                    double minDeltaR)
{
  if (!zCand) return WCandidate();

  vector<WCandidate> wCands;

  for (ElectronV::const_iterator i = electrons.begin(); 
       i != electrons.end(); ++i)
    if (!areOverlapping(i->patEle(), * zCand.daughter(0)) &&
        !areOverlapping(i->patEle(), * zCand.daughter(1)) &&
        reco::deltaR(* i, * zCand.daughter(0)) > minDeltaR &&
        reco::deltaR(* i, * zCand.daughter(1)) > minDeltaR) {
      wCands.push_back(WCandidate(* i, met));
    }

  for (MuonV::const_iterator i = muons.begin(); 
       i != muons.end(); ++i)
    if (!areOverlapping(* i, * zCand.daughter(0)) &&
        !areOverlapping(* i, * zCand.daughter(1)) &&
        reco::deltaR(* i, * zCand.daughter(0)) > minDeltaR &&
        reco::deltaR(* i, * zCand.daughter(1)) > minDeltaR) {
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

TVector2 getPtDiff(const heep::Ele & e){ 
  TVector2  chosenAlgo( e.p4().px(), e.p4().py() );
  LorentzVector p4Def = e.patEle().p4();
  TVector2 defaultAlgo( p4Def.Px(), p4Def.Py() );
  return chosenAlgo - defaultAlgo;
}

TVector2
adjustPt(const ElectronV & electrons){
  TVector2 diff(0.,0.);
  for (ElectronV::const_iterator i = electrons.begin(); 
       i != electrons.end(); ++i){
    diff = diff + getPtDiff(*i);
  }
  return diff;
}

TVector2
adjustPt(const MuonV & muons){
  TVector2 diff(0.,0.);
  for (MuonV::const_iterator i = muons.begin(); 
       i != muons.end(); ++i){
    diff = diff + i->getPtDiff();
  }
  return diff;
}

pat::MET AdjustedMET(const ElectronV & electrons,
                     const pat::MET & met){
  TVector2 adj = adjustPt(electrons);
  TVector2 newmet(met.px()-adj.Px(), met.py()-adj.Py());
  return pat::MET(reco::MET(LorentzVector(newmet.Px(), newmet.Py(), 0., newmet.Mod()), reco::MET::Point(0,0,0)));
//Note: This should include a change of sumET for significance measurements
//  return pat::MET(reco::MET(met.sumEt()+dSumEt, LorentzVector(newmet.Px(), newmet.Py(), 0., newmet.Mod()), reco::MET::Point(0,0,0)));
//Note: Should the new met be wrt beamspot??, what is old met wrt?
//  pat::MET scaledMET(reco::MET(met.sumEt()+dSumEt, reco::MET::LorentzVector(scaledMETPx, scaledMETPy, 0, sqrt(scaledMETPx*scaledMETPx+scaledMETPy*scaledMETPy)), reco::MET::Point(0,0,0)));
}

pat::MET AdjustedMET(const MuonV & muons,
                     const pat::MET & met){
  TVector2 adj = adjustPt(muons);
  TVector2 newmet(met.px()-adj.Px(), met.py()-adj.Py());
  return pat::MET(reco::MET(LorentzVector(newmet.Px(), newmet.Py(), 0., newmet.Mod()), reco::MET::Point(0,0,0)));
//Note: This should include a change of sumET for significance measurements
//  return pat::MET(reco::MET(met.sumEt()+dSumEt, LorentzVector(newmet.Px(), newmet.Py(), 0., newmet.Mod()), reco::MET::Point(0,0,0)));
}
 
pat::MET
AdjustedLeptonMET(const ElectronV & electrons,
                  const MuonV & muons,
                  const pat::MET & met){
  pat::MET met1 = AdjustedMET(electrons, met);
  pat::MET met2 = AdjustedMET(muons    , met1);

  return met2;
}

