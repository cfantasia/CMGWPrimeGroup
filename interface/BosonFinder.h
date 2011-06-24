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

typedef unsigned int uint;
typedef std::vector<int> vint;
typedef std::vector<std::string> vstring;
typedef math::XYZPoint Point;
typedef math::XYZTLorentzVector LorentzVector;
typedef edm::ParameterSet PSet;
typedef edm::MergeableCounter Counter;

//typedef vector<pat::Electron> ElectronV;
//typedef vector<pat::Muon    > MuonV;
typedef std::vector<heep::Ele > ElectronV;
typedef std::vector<TeVMuon  > MuonV;

typedef std::vector<pat::Jet     > JetV;
typedef std::vector<pat::MET     > METV;
typedef std::vector<reco::Track  > TrackV;
typedef std::vector<edm::InputTag> VInputTag;
typedef std::vector<reco::Candidate> CandV;
typedef std::vector<reco::GenParticle> GenParticleV;

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
  ZCandidate() {genLepton1_ = genLepton2_ = 0; leptonic_ = false;}
  ZCandidate(const heep::Ele & p1, const heep::Ele & p2) {
    genLepton1_ = p1.patEle().genLepton();
    genLepton2_ = p2.patEle().genLepton();
    addDaughters(p1.patEle(), p2.patEle());
    leptonic_ = true;
  }
  ZCandidate(const TeVMuon & p1, const TeVMuon & p2) {
    genLepton1_ = p1.genLepton();
    genLepton2_ = p2.genLepton();
    addDaughters(p1, p2);
    leptonic_ = true;
  }
  ZCandidate(const pat::Jet & jet) {
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
 private:
  const reco::Candidate * genLepton1_;
  const reco::Candidate * genLepton2_;
};

class WCandidate : public BosonCandidate {
 public:
  WCandidate() {genLepton_ = 0; leptonic_ = false; mt_ = -999.9;}
  WCandidate(const heep::Ele & lepton, const reco::Candidate & met) {
    genLepton_ = lepton.patEle().genLepton();
    addDaughters(lepton.patEle(), met);
    leptonic_ = true;
    mt_ = CalcMT();
  }
  WCandidate(const TeVMuon & lepton, const reco::Candidate & met) {
    genLepton_ = lepton.genLepton();
    addDaughters(lepton, met);
    leptonic_ = true;
    mt_ = CalcMT();
  }
  WCandidate(const pat::Jet & jet) {
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
  int genMotherId() const {return findGenMotherId(genLepton_);}
  double mt() const { return mt_;}
 private:
  const reco::Candidate * genLepton_;
  double CalcMT(){return sqrt(2 * daughter(0)->et() * daughter(1)->et() * CalcDPhi());}
  double CalcDPhi(){return 1 - cos(reco::deltaPhi(daughter(0)->phi(), daughter(1)->phi()));}
  double mt_;
};

class WZCandidate {

 public:

  WZCandidate() {initialize_();}

  WZCandidate(const ZCandidate & Z, const WCandidate & W) {

    initialize_();
    const reco::Candidate * wLep  = W.daughter(0);
    const reco::Candidate * met   = W.daughter(1);
    const reco::Candidate * zLep1 = Z.daughter(0);
    const reco::Candidate * zLep2 = Z.daughter(1);

    LorentzVector trilepP4 = zLep1->p4() + zLep2->p4() + wLep->p4();
    double term1 = sqrt(trilepP4.M2() + pow(trilepP4.pt(),2)) + met->pt();
    double term2 = (trilepP4 + met->p4()).pt();
    transMass_ = sqrt(pow(term1, 2) - pow(term2, 2));

    double dPhi         = deltaPhi(wLep->phi(), met->phi());
    double g            = (WMASS * WMASS / 2. + 
                           wLep->pt() * met->pt() * cos(dPhi));
    double a            = - pow(wLep->pt(), 2);
    double b            = 2 * g * wLep->pz();
    double c            = pow(g, 2) - pow(wLep->p(), 2) * pow(met->pt(), 2);
    double discriminant = (b * b) - (4 * a * c);
    
    if (discriminant > 0) {
    
      TVector3 leptonP3(wLep->px(), wLep->py(), wLep->pz());
    
      double   pz1 = -b/(2*a) + sqrt(discriminant)/(2*a);
      TVector3 p1  = TVector3(met->px(), met->py(), pz1);
      double   dr1 = p1.DeltaR(leptonP3);
      
      double   pz2 = -b/(2*a) - sqrt(discriminant)/(2*a);
      TVector3 p2  = TVector3(met->px(), met->py(), pz2);
      double   dr2 = p2.DeltaR(leptonP3);

      neutrinoPz_[0] = (dr1 < dr2) ? pz1 : pz2;
      neutrinoPz_[1] = (dr1 < dr2) ? pz2 : pz1;
      neutrinoPz_[2] = (fabs(pz1) > fabs(pz2)) ? pz1 : pz2;
      neutrinoPz_[3] = (fabs(pz1) > fabs(pz2)) ? pz2 : pz1;

      for (size_t i = 0; i < neutrinoPz_.size(); i++) {
        double        px = met->px();
        double        py = met->py();
        double        pz = neutrinoPz_[i];
        double        E  = sqrt(pow(px,2) + pow(py,2) + pow(pz,2));
        LorentzVector p4 = LorentzVector(px, py, pz, E);
        invariantMass_[i] = (Z.p4() + wLep->p4() + p4).M();
      }

    }

  }

  double transMass() {return transMass_;}
  double mt() {return transMass_;}

  double neutrinoPz(std::string type) {
    return neutrinoPz_[index_(type)];
  }

  double mass(std::string type) {
    return invariantMass_[index_(type)];
  }

 private:

  double transMass_;
  std::vector<double> neutrinoPz_;
  std::vector<double> invariantMass_;

  void initialize_() {
    transMass_ = 0.;
    neutrinoPz_ = std::vector<double>(4, 0.);
    invariantMass_ = std::vector<double>(4, 0.);
  }

  size_t index_(std::string type) {
    if (type == "minAngle") return 0;
    if (type == "maxAngle") return 1;
    if (type == "maxPz"   ) return 2;
    if (type == "minPz"   ) return 3;
    printf("Used a non-existent index in WZCandidate\n");
    return 9;
  }

};

typedef std::vector<ZCandidate > ZCandV;
typedef std::vector<WCandidate > WCandV;
typedef std::vector<WZCandidate> WZCandV;

ZCandV getZCands(const ElectronV & electrons, float maxMassDiff = ZMASS);
ZCandV getZCands(const MuonV & muons, float maxMassDiff = ZMASS);
ZCandV getZCands(const ElectronV & electrons,
                 const MuonV & muons, float maxMassDiff = ZMASS);
void removeOverlapping(ZCandV & zCands);

std::vector<WCandidate> getWCandidates(const ElectronV & electrons,
				       const pat::MET & met);
std::vector<WCandidate> getWCandidates(const ElectronV & muons,
				       const pat::MET & met);
std::vector<WCandidate> getWCandidates(const ElectronV & electrons,
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
WCandidate getWCand(const ElectronV & electrons, 
                    const pat::MET & met);
WCandidate getWCand(const MuonV & muons, 
                    const pat::MET & met);
WCandidate getWCand(const JetV & jets);

TVector2 getPtDiff(heep::Ele & e);

TVector2 adjustPt(const ElectronV & electrons);
TVector2 adjustPt(const MuonV & muons);

pat::MET AdjustedMET(const ElectronV & electrons,
                            const pat::MET & met);

pat::MET AdjustedMET(const MuonV & muons,
                             const pat::MET & met);

pat::MET AdjustedMET(const ElectronV & electrons,
                     const MuonV & muons,
                     const pat::MET & met);

#endif
