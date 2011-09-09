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
  bool isLeptonic() const {return leptonic_;}
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
  bool leptonic_;
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
  const reco::Candidate * genLepton1() const {return genLepton1_;}  
  const reco::Candidate * genLepton2() const {return genLepton2_;}
  const reco::Candidate * jet() const {return daughter(0);}
  int genMotherId1() const {return findGenMotherId(genLepton1_);}
  int genMotherId2() const {return findGenMotherId(genLepton2_);}
  const heep::Ele * elec1() const {return elec1_;}
  const heep::Ele * elec2() const {return elec2_;}
  const TeVMuon * muon1() const {return muon1_;}
  const TeVMuon * muon2() const {return muon2_;}
 private:
  const reco::Candidate * genLepton1_;
  const reco::Candidate * genLepton2_;
  const heep::Ele * elec1_;
  const heep::Ele * elec2_;
  const TeVMuon * muon1_;
  const TeVMuon * muon2_;
};

class WCandidate : public BosonCandidate {
 public:
  WCandidate() {genLepton_ = 0; leptonic_ = false; mt_ = -999.9; elec_=NULL; muon_=NULL;}
  WCandidate(const heep::Ele & lepton, const reco::Candidate & met) :
    elec_(&lepton), muon_(NULL){
    genLepton_ = lepton.patEle().genLepton();
    addDaughters(lepton.patEle(), met);
    leptonic_ = true;
    mt_ = calcMT();
  }
  WCandidate(const TeVMuon & lepton, const reco::Candidate & met) :
    elec_(NULL), muon_(&lepton)  {
    genLepton_ = lepton.genLepton();
    addDaughters(lepton, met);
    leptonic_ = true;
    mt_ = calcMT();
  }
  WCandidate(const pat::Jet & jet) :
    elec_(NULL), muon_(NULL){
    addDaughter(jet);
    AddFourMomenta addP4;
    addP4.set(* this);
    leptonic_ = false;
    mt_ = -999.9;
 }
  const reco::Candidate * lepton() const {return daughter(0);}
  const reco::Candidate * met() const {return daughter(1);}
  const reco::Candidate * genLepton() const {return genLepton_;}
  const reco::Candidate * jet() const {return daughter(0);}
  
  const heep::Ele * elec() const {return elec_;}
  const TeVMuon * muon() const {return muon_;}

  int genMotherId() const {return findGenMotherId(genLepton_);}
  double mt() const { return mt_;}
  double calcDPhi(){return reco::deltaPhi(daughter(0)->phi(), daughter(1)->phi());}
 private:
  const reco::Candidate * genLepton_;
  double calcMT(){return sqrt(2 * daughter(0)->et() * daughter(1)->et() * OneMinusCosine());}
  double OneMinusCosine(){return 1 - cos(calcDPhi());}
  double mt_;
  const heep::Ele * elec_;
  const TeVMuon * muon_;
};

class XWLeptonic {
public:
  XWLeptonic(){initialize_();}
  template<class T>
  XWLeptonic(const T& X, const WCandidate & W) {
    initialize_();
    if(!X || !W) return;
    metp4_ = W.daughter(1)->p4();
    setTransMass(X.p4(), W.daughter(0)->p4(), W.daughter(1)->p4());
    setPt(X.p4(), W.p4());//ToDo: Should take 2D Vectors for clarity
    setMassSolns(X.p4(), W.daughter(0)->p4(), W.daughter(1)->p4());
  }

  template<class T>
  XWLeptonic(const T& X, const XWLeptonic & XW) {
    initialize_();
    if(!X || !XW) return;
    LorentzVector tmp = XW.p4() - XW.metp4();
    setTransMass(X.p4(), tmp, XW.metp4());
    setPt(X.p4(), XW.p4());
    setMassSolns(X.p4(), tmp, XW.metp4());
  }
  
  double transMass() const {return transMass_;}
  double mt() const {return transMass_;}
  double pt() const {return pt_;}
  LorentzVector metp4() const {return metp4_;}

  double neutrinoPz(const std::string& type="minPz") const{
    return neutrinoPz_[index_(type)];
  }

  LorentzVector p4(const std::string& type="minPz") const{
    return p4_[index_(type)];
  }

  double mass(const std::string& type="minPz") const{
    return invariantMass_[index_(type)];
  }

  operator bool() const {
    return (invariantMass_[0] > 0);
  }

protected:

  double transMass_;
  double pt_;
  LorentzVector metp4_;
  std::vector<double> neutrinoPz_;
  std::vector<LorentzVector> p4_;
  std::vector<double> invariantMass_;

  void initialize_() {
    transMass_ = 0.;
    pt_ = 0;
    metp4_ = LorentzVector();
    neutrinoPz_ = std::vector<double>(4, 0.);
    p4_ = std::vector<LorentzVector>(4, LorentzVector());
    invariantMass_ = std::vector<double>(4, 0.);
  }

  size_t index_(const std::string& type) const{
    if (type == "minAngle") return 0;
    if (type == "maxAngle") return 1;
    if (type == "maxPz"   ) return 2;
    if (type == "minPz"   ) return 3;
    printf("Used a non-existent index in DiBoson Candidate\n");
    return 9;
  }

  void setTransMass(const LorentzVector & X, const LorentzVector & wLep, const LorentzVector & met){
    LorentzVector XplusLep = X + wLep;
    double term1 = sqrt(XplusLep.M2() + XplusLep.Perp2()) + met.pt();
    double term2 = (XplusLep + met).pt();
    transMass_ = sqrt(pow(term1, 2) - pow(term2, 2));
  }

  void setPt(const LorentzVector & X, const LorentzVector & W){
    pt_ = (X+W).pt();
  }

  bool setNuSolns(const LorentzVector & wLep, const LorentzVector & met){
    double dPhi         = deltaPhi(wLep.phi(), met.phi());
    double g            = (WMASS * WMASS / 2. + 
                           wLep.pt() * met.pt() * cos(dPhi));
    double a            = - wLep.Perp2();
    double b            = 2 * g * wLep.pz();
    double c            = pow(g, 2) - wLep.P2() * met.Perp2();
    double discriminant = (b * b) - (4 * a * c);
    
    if (discriminant > 0) {
      
      TVector3 leptonP3(wLep.px(), wLep.py(), wLep.pz());
    
      double   pz1 = -b/(2*a) + sqrt(discriminant)/(2*a);
      TVector3 p1  = TVector3(met.px(), met.py(), pz1);
      double   dr1 = p1.DeltaR(leptonP3);
      
      double   pz2 = -b/(2*a) - sqrt(discriminant)/(2*a);
      TVector3 p2  = TVector3(met.px(), met.py(), pz2);
      double   dr2 = p2.DeltaR(leptonP3);

      neutrinoPz_[0] = (dr1 < dr2) ? pz1 : pz2;
      neutrinoPz_[1] = (dr1 < dr2) ? pz2 : pz1;
      neutrinoPz_[2] = (fabs(pz1) > fabs(pz2)) ? pz1 : pz2;
      neutrinoPz_[3] = (fabs(pz1) > fabs(pz2)) ? pz2 : pz1;

      return true;
    }
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

  bool setMassSolns(const LorentzVector & X, const LorentzVector & wLep, const LorentzVector & met){
    if( !setP4Solns(X, wLep, met) ) return false;
    for (size_t i = 0; i < neutrinoPz_.size(); i++) {
      invariantMass_[i] = p4_[i].M();
    }
    return true;
  }
  
};

class VZCandidate : public reco::CompositeCandidate{
public:
  VZCandidate(){};
  VZCandidate(const ZCandidate & Z, const WCandidate & W){
    if(!Z || !W) return;
    addDaughter(Z);
    addDaughter(W);
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

void removeWorstCands(ZCandV & zCands, const float & maxMassDiff);
void removeWorstCands(ZCandV & zCands, const float& min, const float& max);
void removeLowLepPtCands(ZCandV & zCands, const float& minPt1, const float& minPt2);
void removeOverlapping(ZCandV & zCands);

template<class L>
ZCandV ZCands(const std::vector<L> & leptons){
  ZCandV zCands;
  for (size_t i = 0; i < leptons.size(); i++)
    for (size_t j = i + 1; j < leptons.size(); j++)
      if (leptons[i].charge() != leptons[j].charge())  // get opposite-charge pairs of leptons
        zCands.push_back(ZCandidate(leptons[i], leptons[j]));
  return zCands;
}

template<class L>
ZCandV getZCands(const std::vector<L> & leptons, float maxMassDiff, bool rmOverlap=true)
{
  ZCandV zCands = ZCands(leptons);
  removeWorstCands(zCands, maxMassDiff);
  // Order by difference from Z mass
  sort(zCands.begin(), zCands.end(), closestToZMass());

  if(rmOverlap) removeOverlapping(zCands);

  return zCands;
}
ZCandV getZCands(const ElectronV & electrons, const MuonV & muons, 
                 float maxMassDiff = ZMASS, bool rmOverlap=true);

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
WCandidate getWCand(const JetV & jets);

#endif
