#ifndef BOSONFINDER_H
#define BOSONFINDER_H

#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <TROOT.h>
#include <TH1F.h>
#include <TVector3.h>
#include <TFile.h>
#include <TSystem.h>

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"

#include "UserCode/CMGWPrimeGroup/interface/TeVMuon.h"
#include "SHarper/HEEPAnalyzer/interface/HEEPEleSelector.h"
#include "SHarper/HEEPAnalyzer/interface/HEEPEle.h"


#include "UserCode/CMGWPrimeGroup/interface/util.h"

class BosonCandidate : public reco::CompositeCandidate {
 public:
  int flavor() const {
    if (numberOfDaughters() > 0)
      return abs(daughter(0)->pdgId());
    return 0;
  }
  operator bool() const {
    return (numberOfDaughters() > 0);
  }

  protected:
  void addDaughters(const reco::Candidate & p1, const reco::Candidate & p2) {
    addDaughter(p1);
    addDaughter(p2);
    AddFourMomenta addP4;
    addP4.set(* this);
  }
  int findGenMotherId(const reco::Candidate * p) const {
    const reco::Candidate * mom = findMother(p);
    if (mom) return mom->pdgId();
    return 0;
  }
};

class ZCandidate : public BosonCandidate {
 public:
  ZCandidate() : 
    elec1_(NULL), elec2_(NULL), muon1_(NULL), muon2_(NULL)
  {genLepton1_ = genLepton2_ = 0; leptonic_ = false;}
  ZCandidate(const heep::Ele & p1, const heep::Ele & p2) :
    elec1_(&p1), elec2_(&p2), muon1_(NULL), muon2_(NULL){
    genLepton1_ = p1.patEle().genLepton();
    genLepton2_ = p2.patEle().genLepton();
    addDaughters(p1.patEle(), p2.patEle());
    leptonic_ = true;
  }
  ZCandidate(const TeVMuon & p1, const TeVMuon & p2) :
    elec1_(NULL), elec2_(NULL), muon1_(&p1), muon2_(&p2){
    genLepton1_ = p1.genLepton();
    genLepton2_ = p2.genLepton();
    addDaughters(p1, p2);
    leptonic_ = true;
  }
  ZCandidate(const pat::Jet & jet) :
    elec1_(NULL), elec2_(NULL), muon1_(NULL), muon2_(NULL){
    addDaughter(jet);
    AddFourMomenta addP4;
    addP4.set(* this);
    leptonic_ = false;
 }
 ZCandidate(const pat::Jet & jet1, const pat::Jet & jet2) :
    elec1_(NULL), elec2_(NULL), muon1_(NULL), muon2_(NULL){
    addDaughters(jet1,jet2);
    AddFourMomenta addP4;
    addP4.set(* this);
    leptonic_ = false;
 }
  const reco::Candidate * genLepton1() const {return genLepton1_;}  
  const reco::Candidate * genLepton2() const {return genLepton2_;}
  const reco::Candidate * jet() const {return daughter(0);}
  int genMotherId1() const {return findGenMotherId(genLepton1_);}
  int genMotherId2() const {return findGenMotherId(genLepton2_);}
  const heep::Ele * elec1() const {return elec1_;}
  const heep::Ele * elec2() const {return elec2_;}
  const TeVMuon * muon1() const {return muon1_;}
  const TeVMuon * muon2() const {return muon2_;}

  bool isLeptonic() const {return leptonic_;}
 private:
  const reco::Candidate * genLepton1_;
  const reco::Candidate * genLepton2_;
  const heep::Ele * elec1_;
  const heep::Ele * elec2_;
  const TeVMuon * muon1_;
  const TeVMuon * muon2_;

  bool leptonic_;
};

//This class is for W's decaying leptonically
//anything with a inv. mass should be a ZCandidate
class WCandidate : public BosonCandidate {
 public:
  WCandidate() {genLepton_ = 0; mt_ = -999.9; elec_=NULL; muon_=NULL;}
  WCandidate(const heep::Ele & lepton, const reco::Candidate & met) :
    elec_(&lepton), muon_(NULL){
    genLepton_ = lepton.patEle().genLepton();
    addDaughters(lepton.patEle(), met);
    mt_ = calcMT();
  }
  WCandidate(const TeVMuon & lepton, const reco::Candidate & met) :
    elec_(NULL), muon_(&lepton)  {
    genLepton_ = lepton.genLepton();
    addDaughters(lepton, met);
    mt_ = calcMT();
  }
    WCandidate(const reco::Candidate & GenLepton, const reco::Candidate & GenNeutrino) :
      elec_(NULL), muon_(NULL) {
      addDaughters(GenLepton, GenNeutrino);
      mt_ = calcMT();
    }

  const reco::Candidate * lepton() const {return daughter(0);}
  const reco::Candidate * met() const {return daughter(1);}
  const reco::Candidate * genLepton() const {return genLepton_;}
  
  const heep::Ele * elec() const {return elec_;}
  const TeVMuon * muon() const {return muon_;}

  int genMotherId() const {return findGenMotherId(genLepton_);}
  double mt() const { return mt_;}
  double mass() const {std::cout<<" You asked for the mass of a leptonic W!!!\nWe don't know its value"<<std::endl; return -1;}
  double calcDPhi(){return reco::deltaPhi(daughter(0)->phi(), daughter(1)->phi());}
 private:
  const reco::Candidate * genLepton_;
  double calcMT(){return sqrt(2 * daughter(0)->et() * daughter(1)->et() * OneMinusCosine());}
  double OneMinusCosine(){return 1 - cos(calcDPhi());}
  double mt_;
  const heep::Ele * elec_;
  const TeVMuon * muon_;
};

enum NuAlgos { kMinPz, kMaxPz, kMinDR, kMaxDR, kMinDRW, kMaxDRW, kMinTheta, kMaxTheta };
#define kNuAlgos 8
class XWLeptonic {
public:
  XWLeptonic(){initialize_();}
  template<class T>
  XWLeptonic(const T& X, const WCandidate & W) {
    initialize_();
    if(!X.mass() || !W) return;
    if(!setP4Solns(X.p4(), W.daughter(0)->p4(), W.daughter(1)->p4())) 
      initialize_();
  }

  template<class T>
  XWLeptonic(const T& X, const XWLeptonic & XW) {
    initialize_();
    if(!X.mass() || !XW) return;
    p4_[0] = X.p4() + XW.p4_[0];
    p4_[1] = X.p4() + XW.p4_[1];
    neutrinoPz_ = XW.neutrinoPz_;
    soln_ = XW.soln_;
  }
  
  inline double neutrinoPz(const NuAlgos& type) const{
    return neutrinoPz_[soln_[type]];
  }

  inline const LorentzVector & p4(const NuAlgos& type) const{
    return p4_[soln_[type]];
  }

  inline operator bool() const {
    return p4_[0].mass()>0;
  }

  inline const LorentzVector & operator ()(const NuAlgos& type=kMinPz) const{
    return p4(type);
  }

  double discriminant() const{
    return discriminant_;
  }
  double discriminantFrac() const{
    return discriminantFrac_;
  }
  double discriminantAngle() const{
    return discriminantAngle_;
  }
  double discriminantReal() const{
    return discriminantReal_;
  }
  double discriminantImag() const{
    return discriminantImag_;
  }
  
protected:

  std::vector<double> neutrinoPz_;
  std::vector<LorentzVector> p4_;
  std::vector<bool> soln_;
  double discriminant_;
  double discriminantFrac_;
  double discriminantAngle_;
  double discriminantReal_;
  double discriminantImag_;

  void initialize_() {
    neutrinoPz_ = std::vector<double>(2, 0.);
    p4_ = std::vector<LorentzVector>(2, LorentzVector());
    soln_ = std::vector<bool>(kNuAlgos, false);
  }

  bool setNuSolns(const LorentzVector & wLep, const LorentzVector & met){
    float mt = sqrt(2 * wLep.Et() * met.Et() * ( 1-cos(reco::deltaPhi(wLep.phi(),met.phi())) ) );
    double motherMass = std::max(WMASS, mt);

    double dPhi         = deltaPhi(wLep.phi(), met.phi());
    double g            = (motherMass * motherMass / 2. + 
                           wLep.pt() * met.pt() * cos(dPhi));
    double a            = - wLep.Perp2();
    double b            = 2 * g * wLep.pz();
    double c            = pow(g, 2) - wLep.P2() * met.Perp2();
    double term1        = b/(2*a);
    double discriminant = (b * b) - (4 * a * c);
    double discFrac     = term1*term1 - c/a;

    discriminant_ = discriminant;
    discriminantFrac_ = discFrac;
    discriminantReal_ = -b + (discriminant > 0 ? sqrt(discriminant) : 0);
    discriminantImag_ = discriminant > 0 ? 0 :  sqrt(-discriminant);

    discriminantAngle_ = atan( discriminantImag_ / discriminantReal_ );

    //Need a tolerace in here to deal with floating point math
    if(discFrac < 0 && discFrac > -1.) discFrac = 0.;
    if (discFrac >= 0) {
      
      TVector3 leptonP3(wLep.px(), wLep.py(), wLep.pz());
    
      double   pz1 = -term1 + sqrt(discFrac);
      TVector3 p1  = TVector3(met.px(), met.py(), pz1);
      TVector3 w1  = leptonP3 + p1;
      double   dr1 = p1.DeltaR(leptonP3);
      double   drW1= p1.DeltaR(w1);
      double   dt1 = fabs(p1.Theta() - leptonP3.Theta());
      
      double   pz2 = -term1 - sqrt(discFrac);
      TVector3 p2  = TVector3(met.px(), met.py(), pz2);
      TVector3 w2  = leptonP3 + p2;
      double   dr2 = p2.DeltaR(leptonP3);
      double   drW2= p2.DeltaR(w2);
      double   dt2 = fabs(p2.Theta() - leptonP3.Theta());

      neutrinoPz_[0] = pz1;
      neutrinoPz_[1] = pz2;

      soln_[kMinPz] = fabs(pz1) > fabs(pz2); 
      soln_[kMaxPz] = !soln_[kMinPz];
      soln_[kMaxDR] = dr1 < dr2;
      soln_[kMinDR] = !soln_[kMaxDR];
      soln_[kMaxDRW]= drW1 < drW2;
      soln_[kMinDRW]= !soln_[kMaxDRW];
      soln_[kMaxTheta] = dt1 < dt2;
      soln_[kMinTheta] = !soln_[kMaxTheta];
      return true;
    }

/*
    std::cout<<" mother mass: "<<motherMass
             <<" mt:" <<mt
             <<" dphi:" <<dPhi
             <<" g:" <<g
             <<" a:" <<a
             <<" b:" <<b
             <<" c:" <<c
             <<" b*b:"<<b*b
             <<" 4ac: "<<4*a*c
             <<" first: "<<term1*term1
             <<" 2nd  : "<<c/a
             <<" disc:" <<discriminant
             <<" discFrac:" <<discFrac
             <<"\n";
*/
    return false;
  }

  bool setP4Solns(const LorentzVector & X, const LorentzVector & wLep, const LorentzVector & met){
    if( !setNuSolns(wLep, met) ) return false;
    for (size_t i = 0; i < neutrinoPz_.size(); i++) {
      double        px = met.px();
      double        py = met.py();
      double        pz = neutrinoPz_[i];
      double        E  = sqrt(pow(px,2) + pow(py,2) + pow(pz,2));
      LorentzVector neu(px, py, pz, E);
      p4_[i] = X + wLep + neu; 
    }
    return true;
  }
  
};

class VZCandidate : public reco::CompositeCandidate{
public:
  VZCandidate(){};
  VZCandidate(const ZCandidate & Z, const ZCandidate & V){
    if(!Z || !V) return;
    addDaughter(Z);
    addDaughter(V);
    AddFourMomenta addP4;
    addP4.set(* this);
  }
  VZCandidate(const pat::Jet & j1, const pat::Jet & j2){
    addDaughter(j1);
    addDaughter(j2);
    AddFourMomenta addP4;
    addP4.set(* this);
  }
  operator bool() const {
    return (numberOfDaughters() > 0);
  }
private:
};

typedef std::vector<ZCandidate > ZCandV;
typedef std::vector<WCandidate > WCandV;
typedef std::vector<XWLeptonic > XWLepV;

void removeWorstCands(ZCandV & zCands, const float & maxMassDiff);
void removeWorstCands(ZCandV & zCands, const float& min, const float& max);
void removeLowLepPtCands(ZCandV & zCands, const float& minPt1, const float& minPt2);
void removeOverlapping(ZCandV & zCands);

template<class L>
ZCandV ZCands(const std::vector<L> & leptons, const std::vector<bool> & mask){
  bool useMask = mask.size(); //Mask let's you ignore certain leptons
  ZCandV zCands;
  for (size_t i = 0; i < leptons.size(); i++){
    if(useMask && !mask[i]) continue;
    for (size_t j = i + 1; j < leptons.size(); j++){
      if(useMask && !mask[j]) continue;
      if (leptons[i].charge() != leptons[j].charge())  // get opposite-charge pairs of leptons
        zCands.push_back(ZCandidate(leptons[i], leptons[j]));
    }
  }
  return zCands;
}

template<class L>
ZCandV getZCands(const std::vector<L> & leptons, float maxMassDiff = ZMASS, bool rmOverlap=true, const std::vector<bool> & mask=std::vector<bool>())
{
  ZCandV zCands = ZCands(leptons, mask);
  removeWorstCands(zCands, maxMassDiff);
  // Order by difference from Z mass
  sort(zCands.begin(), zCands.end(), closestToZMass());

  if(rmOverlap) removeOverlapping(zCands);

  return zCands;
}
ZCandV getZCands(const ElectronV & electrons, const MuonV & muons, 
                 float maxMassDiff = ZMASS, bool rmOverlap=true,
                 const std::vector<bool> & eMask=std::vector<bool>(), 
                 const std::vector<bool> & mMask=std::vector<bool>());

template<class L>
WCandV getWCandidates(const std::vector<L> & leptons, const pat::MET & met){
  WCandV wCands;
  for (size_t i = 0; i < leptons.size(); i++)
    wCands.push_back(WCandidate(leptons[i], met));
  
  return wCands;
}

template<class L>
WCandidate getWCand(const std::vector<L> & leptons, const pat::MET & met){
  WCandV wCands = getWCandidates(leptons, met);

  sort(wCands.begin(), wCands.end(), highestPtLepton());
  if (wCands.size()) return wCands[0];

  return WCandidate();
}

WCandV getWCandidates(const ElectronV & electrons,
                      const MuonV & muons, 
                      const pat::MET & met);

WCandidate getWCand(const ElectronV & electrons,
                    const MuonV & muons, 
                    const pat::MET & met);
WCandidate getWCand(const ElectronV & electrons,
                    const MuonV & muons, 
                    const pat::MET & met,
                    const ZCandidate & zCand,
                    double minDeltaR = 0.);
ZCandidate getVCand(const JetV & jets, const int mode=0);
ZCandidate getVCand2(const JetV & jets);

#endif
