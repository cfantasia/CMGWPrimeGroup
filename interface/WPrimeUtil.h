#ifndef _wprime_util_h_
#define _wprime_util_h_

#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"

#include "UserCode/CMGWPrimeGroup/interface/util.h"
#include "UserCode/CMGWPrimeGroup/interface/BosonFinder.h"

#include <TVector2.h>
#include <TLorentzVector.h>

#include <fstream>
#include <iostream>
#include <string>

//class TFileService;
class TH1D;
class TH2D;


class WPrimeUtil
{
 public:
  WPrimeUtil(edm::InputTag genLabel, edm::InputTag pfLabel, std::string cross_sections);

  ~WPrimeUtil();

  // get input files (to be retrieved from samples_cross_sections.txt)
  void getInputFiles(std::vector<wprime::InputFile> & inputFiles);

  inline void setApplyHadronicRecoilCorrection(bool flag){applyHadronicRecoilCorrection_ = flag;}
  inline void setHadronicMETcalculated(bool flag){hadronicMETcalculated_ = flag;}

  inline bool shouldApplyHadronicRecoilCorrection(){return applyHadronicRecoilCorrection_;}

  inline void resetWarnings(){warningShown_ = false;}
  inline void setSampleName(std::string samplename){samplename_ = samplename;}
  inline void setSampleWeight(float weight){sampleweight_ = weight;}
  inline void setCurrentSample(std::vector<wprime::InputFile>::iterator sample){ currentSample_ = sample;}
  inline void setWeight(float weight){weight_ = weight;}
  inline void setEventsToDebug(const std::vector<edm::EventID>& vEvents){vEventsToDebug_ = vEvents;}

  inline std::string getSampleName() const{return samplename_;}
  inline float getSampleWeight() const{return sampleweight_;}
  inline float getWeight() const{return weight_;}
  inline float getLumiWeight   (const int   & nInt){ return LumiWeights_.weight   (nInt);}
  inline float getLumiWeight3BX(const float & nInt){ return LumiWeights_.weight3BX(nInt);}
  int   getPU1BX(const std::vector< PileupSummaryInfo > & PupInfo);
  float getPUWeight1BX(const std::vector< PileupSummaryInfo > & PupInfo);
  float getPUWeight3BX(edm::EventBase const & event, const std::string& label);
  float getPU3BX(const std::vector< PileupSummaryInfo > & PupInfo);
  float getPUWeight3BX(const std::vector< PileupSummaryInfo > & PupInfo);
  float getPUWeight3D(const std::vector< PileupSummaryInfo > & PupInfo);

  static void getEff(float & eff, float & deff,float Num,float Denom);

  inline std::vector<wprime::InputFile>::iterator getCurrentSample(){ return currentSample_;}

  // true if current file under processing contains "data" in name
  inline bool runningOnData() const{return runningOnData_;};
  inline void setRunningOnData(){ runningOnData_ = (getSampleName().find("data") != std::string::npos);};
  // true if current file under processing is signal sample
  inline bool isSignalSample() const {return isSignalSample_;}
  inline void setIsSignalSample(bool flag){isSignalSample_ = flag;}

  //Check if Run/Evt is in Debug list
  bool DebugEvent(edm::EventBase const& event) const;

  // integrated luminosity in pb^-1
  inline float getLumi_ipb(){return lumi_ipb;}

  static inline bool SameTrigger(const std::string & versionedName, const std::string & wildcardedName){
    size_t k = wildcardedName.find("*");
    return (k == std::string::npos) ? 
      !versionedName.compare(wildcardedName) : //No '*', early triggers
      !versionedName.compare(0, k, wildcardedName, 0, k);
  }

  // return pointer to gen-particle with pdgId and mother pdgId_mother
  const reco::Candidate * getGenParticle(edm::EventBase const & event, int pdgId, int pdgId_mother); 

  static void printEvent(edm::EventBase const & event);

  static bool passTriggersCut(edm::EventBase const & event, std::string label,const std::vector<std::string>& triggerNames);
  static bool passTriggersCut(const pat::TriggerEvent & triggerEvent,const std::vector<std::string>& triggerNames);
  static void printPassingTriggers(const pat::TriggerEvent & triggerEvent,const std::vector<std::string>& triggerNames);
  static bool FoundAndpassed(const pat::TriggerEvent & triggerEvent, const pat::TriggerPath& path, const std::vector<std::string>& triggerNames);
  static bool passed(const pat::TriggerEvent & triggerEvent, const pat::TriggerPath& path);
  static unsigned L1Prescale(const pat::TriggerEvent & triggerEvent, const pat::TriggerPath& path);
  static unsigned MaxL1Prescale(const pat::TriggerEvent & triggerEvent, const pat::TriggerPath& path);
  static bool FindTrigger(const pat::TriggerPath& path, const std::vector<std::string>& triggerNames);

  // get hadronic MET component (that needs to be corrected 
  // if applyMETCorrection=true)from Z data; this will be done according to hadronic 
  // activity from Z->mumu reconstructed events
  TVector2 getHadronicMET(edm::EventBase const & event);

  void setLumiWeights(const std::string & MCFile, const std::string & DataFile, 
                      const std::string & MCHist, const std::string & DataHist);
  static void CheckStream(const ofstream& stream, const std::string & s);
  
  /////////////////////
  ////Adjust Pt Fns////
  /////////////////////
  
  template<class T>
    static TVector2
    adjustPt(const std::vector<T>& leptons, const std::vector<pat::PFParticle>& pfCands ){
    TVector2 diff(0.,0.);
    for (uint i=0; i<leptons.size(); ++i){
      int pfCandIdx = FindPFCand(leptons[i], pfCands);
      if(pfCandIdx == -1) {continue;} //PF Obj not found
      diff = diff + getPtDiff(leptons[i],pfCands[pfCandIdx]);
    }
    return diff;
  }
  
  template<class T>
    static int 
    FindPFCand(const T & lep, const std::vector<pat::PFParticle>& pfCands){
    for(unsigned i=0; i<pfCands.size(); ++i){
      if(Match(lep,pfCands[i])) return i;
    }
    return -1;
  }
  
  static TVector2 getPtDiff(const heep::Ele & e, const pat::PFParticle & pfCand){
    TVector2  chosenAlgo( e.p4().px(), e.p4().py() );
    TVector2 defaultAlgo( pfCand.px(), pfCand.py() );
    return chosenAlgo - defaultAlgo;
  }

  template<class T>
    static TVector2 getPtDiff(const T & p, const pat::PFParticle & pfCand){
    TVector2  chosenAlgo(      p.px(),      p.py() );
    TVector2 defaultAlgo( pfCand.px(), pfCand.py() );
    return chosenAlgo - defaultAlgo;
  }

  // adjust MET <met> by subtracting <subtract>; 
  // takes MET significance properly into account by adjusting sumEt
  static void AdjustMET(pat::MET & met, const TVector2 & subtract)
    {
      TVector2 newmet(met.px()-subtract.Px(), met.py()-subtract.Py());
      //      std::cout<<"Before met et: "<<met.et()<<" met phi: "<<met.phi()<<std::endl;
      //Note: Should the new met be wrt beamspot??, what is old met wrt?
      met = pat::MET(reco::MET(met.sumEt()+subtract.Mod(), 
			       LorentzVector(newmet.Px(), newmet.Py(), 0., newmet.Mod()), 
			       reco::MET::Point(0,0,0)));
      //      std::cout<<"after met et: "<<met.et()<<" met phi: "<<met.phi()<<std::endl;
    }
  
/////////////////////
////Adjust MET Fns///
/////////////////////
  template<class T>
    static void AdjustMET(const std::vector<T> & leptons,
                          const std::vector<pat::PFParticle> & pfCands,
                          pat::MET & met){
    AdjustMET(met, adjustPt(leptons, pfCands));
  }
  
  template<class T>
    static void AdjustMET(edm::EventBase const & event, const std::vector<T>& leptons, const edm::InputTag& pfLabel,  pat::MET & met){
    std::vector<pat::PFParticle> pfCands;
    getPFCands(event, pfLabel, pfCands);
    AdjustMET(leptons, pfCands, met);
  }

  static void AdjustMET(edm::EventBase const & event, 
                        const ElectronV & electrons, const MuonV & muons,
                        const edm::InputTag& pfLabel,  pat::MET & met);

/*
//These functions don't work as written, bc the patEle go out of scope
  static void getElectronsMET(edm::EventBase const & event,
                              const edm::InputTag& eLabel, ElectronV & electrons,
                              const edm::InputTag& metLabel, const bool & adjMET, pat::MET & met,
                              const edm::InputTag& pfLabel);
  static void getMuonsMET(edm::EventBase const & event,
                          const edm::InputTag& muLabel, const uint& muAlgo, MuonV & muons,
                          const edm::InputTag& metLabel, const bool & adjMET, pat::MET & met,
                          const edm::InputTag& pfLabel);
  static void getLeptonsMET(edm::EventBase const & event, 
                            const edm::InputTag& eLabel, ElectronV & electrons,
                            const edm::InputTag& muLabel, const uint& muAlgo, MuonV & muons,
                            const edm::InputTag& metLabel, const bool & adjMET, pat::MET & met, 
                            const edm::InputTag& pfLabel);
*/
  static void getElectronsMET(edm::EventBase const & event,
                              const std::vector<pat::Electron>& patElectrons, ElectronV & electrons,
                              const edm::InputTag& metLabel, const bool & adjMET, pat::MET & met,
                              const edm::InputTag& pfLabel);
  static void getMuonsMET(edm::EventBase const & event,
                          const std::vector<pat::Muon    >& patMuons, const uint& muAlgo, MuonV & muons,
                          const edm::InputTag& metLabel, const bool & adjMET, pat::MET & met,
                          const edm::InputTag& pfLabel);
  static void getLeptonsMET(edm::EventBase const & event, 
                            const std::vector<pat::Electron>& patElectrons, ElectronV & electrons,
                            const std::vector<pat::Muon    >& patMuons, const uint& muAlgo, MuonV & muons,
                            const edm::InputTag& metLabel, const bool & adjMET, pat::MET & met,
                            const edm::InputTag& pfLabel);
    
  static void getElectronsMET(const PatElectronVH & patElectronsH, ElectronV & electrons,
                              const METVH & metH, const bool & adjMET, pat::MET & met,
                              const PFCandidateVH & pfCandidatesH);
  static void getMuonsMET(const PatMuonVH & patMuonsH, const uint& muAlgo, MuonV & muons,
                          const METVH & metH, const bool & adjMET, pat::MET & met,
                          const PFCandidateVH & pfCandidatesH);
  static void getLeptonsMET(const PatElectronVH & patElectronsH, ElectronV & electrons,
                            const PatMuonVH & patMuonsH, const uint& muAlgo, MuonV & muons,
                            const METVH & metH, const bool & adjMET, pat::MET & met,
                            const PFCandidateVH & pfCandidatesH);


  static void convertElectrons(const std::vector<pat::Electron>& patElectrons, ElectronV & electrons);
  static void convertMuons(const std::vector<pat::Muon>& patMuons, const uint& muAlgo, MuonV & muons);
  static void getElectrons(const edm::EventBase & event, const edm::InputTag& label, ElectronV & electrons);
  static void getMuons    (const edm::EventBase & event, const edm::InputTag& label, const uint&  muonAlgo, MuonV & muons);
  static void getPFCands  (const edm::EventBase & event, const edm::InputTag& label, std::vector<pat::PFParticle> & pfCands);
  static void getMET      (const edm::EventBase & event, const edm::InputTag& label, pat::MET & met)
    {met = getProduct<METV>(event, label)[0];}

  // MET adjustments, designed for lepton+MET signatures. There are two corrections to be made:
  // (a) the hadronic recoil component (that needs to be adjusted in simulated W->lepton samples 
  // if applyHadronicRecoilCorrection=true) from Z data; this will be done according to hadronic 
  // activity from Z->ll reconstructed events, with the addition of the lepton pt; if enabled, 
  // this ignores completely the measured (pf)MET in the event
  // (b) the lepton-pt component that needs to be updated if we switch to one
  // of the dedicated high-pt TeV or HEEP reconstructors
  template<class T>
    void getNewMET(const edm::EventBase & event, const T & lepton, pat::MET & met, const edm::InputTag & metLabel)
    {
      if(shouldApplyHadronicRecoilCorrection())
	{ // this is correction (a)
	  TVector2 hadronMET = getHadronicMET(event);
	  // initialize event's MET to hadronic-MET
	  met = pat::MET(reco::MET(LorentzVector(hadronMET.Px(), hadronMET.Py(), 0, hadronMET.Mod()), reco::MET::Point(0,0,0)));
	  // adjust by subtracting the lepton pt
	  AdjustMET(met, TVector2(lepton.p4().px(), lepton.p4().py()));
	}
      else
	{ // this is correction (b)
	  std::vector<T> vl; vl.push_back(lepton);
	  getMET(event, metLabel, met);
	  // correct (pf)MET by taking into account TeV/heep reconstruction for muons/electrons
	  AdjustMET(event, vl, pfLabel_, met);
	}
    }

  // calls getNewMET; returns W canidate from lepton and (adjusted) MET
  template<class T>
    WCandidate getNewMETandW(const edm::EventBase & event, const T & lepton, pat::MET & met, const edm::InputTag & metLabel)
    {
      getNewMET(event, lepton, met, metLabel);
      return WCandidate(lepton, met);
    }

/////////////
//Matching///
/////////////
template<class T1,class T2>
static bool Match(const T1 & p1, const T2 & p2){
  float tolerance = 0.01;
  if (p1.pdgId() == p2.pdgId() &&
      fabs(p1.eta() - p2.eta()) < tolerance &&
      fabs(reco::deltaPhi(p1.phi(),p2.phi())) < tolerance
    )
    return true;
  return false;
}

template<class T>
static bool Match(const heep::Ele & p1, const T & p2){
  return Match(p1.patEle(), p2);
}
template<class T>
static bool Match(const T & p1, const heep::Ele & p2){
  return Match(p1, p2.patEle());
}
static bool Match(const heep::Ele & p1, const heep::Ele & p2);

template<class T1, class T2>
static uint FindIndex(const T1 & p, const std::vector<T2>& vec){
  for(uint i=0; i<vec.size(); ++i){
    if(Match(vec[i], p)) return i;
  }
  std::cerr<<"Didn't find match for particle, returning random one!!!\n";
  return 0;
}

template<class T1, class T2>
static const T2 & Find(const T1 & p, const std::vector<T2>& vec){
  for(uint i=0; i<vec.size(); ++i){
    if(Match(vec[i], p)) return vec[i];
  }
  std::cerr<<"Didn't find match for particle, returning random one!!!\n";
  return vec[0];
}

template<class T1, class T2>
static bool Contains(const T1 & p, const std::vector<T2>& vec){
  for(uint i=0; i<vec.size(); ++i){
    if(Match(vec[i], p)) return true;
  }
  return false;
}

////////////////////
//Trigger Matching//
////////////////////

template<class P>
static bool 
passTriggerMatch(const P & p, const float cut, const vstring& triggers){
  for(uint i=0; i<p.triggerObjectMatches().size(); ++i){
    vstring names = p.triggerObjectMatches()[i].pathNames(true, false);
    for(uint j=0; j<names.size(); ++j)
      for (size_t k=0; k < triggers.size(); ++k)
        if(WPrimeUtil::SameTrigger(names[j], triggers[k]))
          if (p.triggerObjectMatchesByPath(names[j], true, false).size() > 0)
            if(p.triggerObjectMatchByPath(names[j], true, false)->pt() > cut) return true;
  }
  return false;
}

static bool 
passTriggerMatch(const heep::Ele & p, const float cut, const vstring& triggers){
  return passTriggerMatch(p.patEle(), cut, triggers);
}

// get GEN-level transverse mass for lepton + neutrino;
// using delta-R matching requirement between RECO-lepton and GEN-lepton
 float getGenWprimeMt(edm::EventBase const& event, int pdgId_lepton,
		      int pdgId_neutrino, const reco::Candidate * lepton);
 
private:

  // directory containing all input samples
  std::string top_level_dir; 

  // file with samples & cross-sections
  std::string sample_cross_sections;

  void setupZMETcorrection();
  void setRecoilProjections();
  
  TH1D * hRecoilPerp;
  TH2D * hRecoilParalvsVBPt;
  TH1D ** histRecoilParal;

  // keep track of input file name and weight (e.g. for scaling MC histograms);
  // values set at beginFile
  std::vector<wprime::InputFile>::iterator currentSample_;
  std::string samplename_;
  float sampleweight_;
  float weight_;
  bool runningOnData_;
  bool isSignalSample_;
  static bool warningShown_;
  std::vector<edm::EventID> vEventsToDebug_;

  float lumi_ipb; // in pb^-1, to be retrieved from samples_cross_sections.txt

  // used for the parsing of samples_cross_sections.txt
  void parseLine(const std::string & new_line, wprime::InputFile * in_file);

  bool applyHadronicRecoilCorrection_; // wether to apply hadronic recoil correction (aimed for W)
  bool hadronicMETcalculated_; // want to calculate this max. once per event

  TVector2 hadronicMETcached;

  edm::InputTag genLabel_;
  edm::InputTag pfLabel_;
 // Handle to GenParticleCollectiom>
  edm::Handle<reco::GenParticleCollection> genParticles;
  
  edm::LumiReWeighting LumiWeights_;
  edm::Lumi3DReWeighting LumiWeights3D_;
  
};

#endif //#define _wprime_util_h_
