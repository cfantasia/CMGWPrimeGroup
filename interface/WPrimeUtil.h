#ifndef _wprime_util_h_
#define _wprime_util_h_

#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "UserCode/CMGWPrimeGroup/interface/util.h"
#include "UserCode/CMGWPrimeGroup/interface/BosonFinder.h"

#include <TVector2.h>
#include <TLorentzVector.h>

#include <fstream>
#include <iostream>
#include <string>

class TFileService;
class TH1D;
class TH2D;


class WPrimeUtil
{
 public:
  WPrimeUtil(const char * out_filename, edm::InputTag genLabel, edm::InputTag pfLabel, std::string cross_sections);

  ~WPrimeUtil();

  fwlite::TFileService * getFileService(){return fs;}

  // get input files (to be retrieved from samples_cross_sections.txt)
  void getInputFiles(std::vector<wprime::InputFile> & inputFiles);

  void setApplyHadronicRecoilCorrection(bool flag){applyHadronicRecoilCorrection_ = flag;}
  void setHadronicMETCalculated(bool flag){hadronicMETcalculated_ = flag;}

  bool shouldApplyHadronicRecoilCorrection(){return applyHadronicRecoilCorrection_;}

  void setSampleName(std::string samplename){samplename_ = samplename;}
  void setWeight(float weight){weight_ = weight;}
  void setInputFile(float weight){weight_ = weight;}
  void SetEventsToDebug(const std::vector<edm::EventID>& vEvents){vEventsToDebug_ = vEvents;}

  std::string getSampleName() const{return samplename_;}
  float getWeight() const{return weight_;}
  float getTotalWeight3BX(edm::EventBase const & event, const std::string& label);
  float getLumiWeight   (const int   & nInt){ return LumiWeights_.weight   (nInt);}
  float getLumiWeight3BX(const float & nInt){ return LumiWeights_.weight3BX(nInt);}
  int   getPU1BX(const std::vector< PileupSummaryInfo > & PupInfo);
  float getPUWeight1BX(const std::vector< PileupSummaryInfo > & PupInfo);
  float getPUWeight3BX(edm::EventBase const & event, const std::string& label);
  float getPU3BX(const std::vector< PileupSummaryInfo > & PupInfo);
  float getPUWeight3BX(const std::vector< PileupSummaryInfo > & PupInfo);

  static void getEff(float & eff, float & deff,float Num,float Denom);

  // true if current file under processing contains "data" in description
  bool runningOnData() const{return runningOnData_;};
  void setRunningOnData();

  //Check if Run/Evt is in Debug list
  bool DebugEvent(edm::EventBase const& event) const;

  // integrated luminosity in pb^-1
  float getLumi_ipb(){return lumi_ipb;}

  static inline bool SameTrigger(const std::string & versionedName, const std::string & wildcardedName){
    return (wildcardedName.find("*") == std::string::npos) ? 
      !versionedName.compare(wildcardedName) : //No '*', early triggers
      !versionedName.compare(0, versionedName.size()-1, wildcardedName, 0, wildcardedName.size()-1);//This assumes v is under 10, need to fix
  }
  static void PrintEvent(edm::EventBase const & event);

  static bool PassTriggersCut(edm::EventBase const & event, std::string label,const std::vector<std::string>& triggerNames);
  static bool PassTriggersCut(const pat::TriggerEvent & triggerEvent,const std::vector<std::string>& triggerNames);
  static void PrintPassingTriggers(const pat::TriggerEvent & triggerEvent,const std::vector<std::string>& triggerNames);
  static bool FoundAndPassed(const pat::TriggerEvent & triggerEvent, const pat::TriggerPathRef path, const std::vector<std::string>& triggerNames);
  static bool Passed(const pat::TriggerEvent & triggerEvent, const pat::TriggerPathRef path);
  static unsigned L1Prescale(const pat::TriggerEvent & triggerEvent, const pat::TriggerPathRef path);
  static unsigned MaxL1Prescale(const pat::TriggerEvent & triggerEvent, const pat::TriggerPathRef path);
  static bool FindTrigger(const pat::TriggerPathRef path, const std::vector<std::string>& triggerNames);

  //transverse mass with a given MET object (TVector2)
  static float TMass(const TLorentzVector& lv, const TVector2& themet)
    {
      //------------------------------------------------------------------------
      float tmass = 0;
      float cdphi = TMath::Cos(lv.Phi()-themet.Phi());
      float tmass_sqr = 2*lv.Pt()*themet.Mod()*(1-cdphi);
      tmass = (tmass_sqr>0) ? sqrt(tmass_sqr) : 0;
    return tmass;
    }

  // get hadronic MET component (that needs to be corrected 
  // if applyMETCorrection=true)from Z data; this will be done according to hadronic 
  // activity from Z->mumu reconstructed events
  TVector2 getHadronicMET(edm::EventBase const & event);

  void SetLumiWeights(const std::string & MCFile, const std::string & DataFile, 
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
      if(pfCandIdx == -1) continue; //PF Obj not found
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
      //std::cout<<"Before met et: "<<met.et()<<" met phi: "<<met.phi()<<std::endl;
      //Note: Should the new met be wrt beamspot??, what is old met wrt?
      met = pat::MET(reco::MET(met.sumEt()+subtract.Mod(), 
			       LorentzVector(newmet.Px(), newmet.Py(), 0., newmet.Mod()), 
			       reco::MET::Point(0,0,0)));
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
  static void getMET      (const edm::EventBase & event, const edm::InputTag& label, pat::MET & met);

  static void tabulateSummary(wprime::EffV& results);
  static void printSummary(const std::string& dir, const std::string& description, const vstring & Cuts, const wprime::EffV& results, ofstream& out);
  
  // MET adjustments, designed for lepton+MET signatures. There are two corrections to be made:
  // (a) the hadronic recoil component (that needs to be adjusted in simulated W->lepton samples 
  // if applyHadronicRecoilCorrection=true) from Z data; this will be done according to hadronic 
  // activity from Z->ll reconstructed events, with the addition of the lepton pt; if enabled, 
  // this ignores completely the measured (pf)MET in the event
  // (b) the lepton-pt component that needs to be updated if we switch to one
  // of the dedicated high-pt TeV or HEEP reconstructors
  template<class T>
    void getNewMET(const edm::EventBase & event, const T & lepton, pat::MET & met)
    {
      if(shouldApplyHadronicRecoilCorrection())
	{ // this is correction (a)
	  TVector2 hadronMET = getHadronicMET(event);
	  // initialize event's MET to hadronic-MET
	  met = pat::MET(reco::MET(LorentzVector(hadronMET.Px(), hadronMET.Py(), 0, hadronMET.Mod()), reco::MET::Point(0,0,0)));
	  // adjust by subtracting the lepton pt
	  AdjustMET(met, TVector2(lepton.px(), lepton.py()));
	}
      else
	{ // this is correction (b)
	  std::vector<T> vl; vl.push_back(lepton);
	  // correct (pf)MET by taking into account TeV/heep reconstruction for muons/electrons
	  AdjustMET(event, vl, pfLabel_, met);
	}
    }

  // calls getNewMET; returns W canidate from lepton and (adjusted) MET
  template<class T>
    WCandidate getNewMETandW(const edm::EventBase & event, const T & lepton, pat::MET & met)
    {
      getNewMET(event, lepton, met);
      return WCandidate(lepton, met);
    }

/////////////
//Matching///
/////////////
template<class T1,class T2>
static bool Match(const T1 & p1, const T2 & p2){
  float tolerance = 0.0001;
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


private:
  fwlite::TFileService * fs;
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
  std::string samplename_;
  float weight_;
  bool runningOnData_;

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
  
};

#endif //#define _wprime_util_h_
