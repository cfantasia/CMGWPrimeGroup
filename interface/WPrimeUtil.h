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
  WPrimeUtil(const char * out_filename, edm::InputTag genParticles, std::string cross_sections);

  ~WPrimeUtil();

  fwlite::TFileService * getFileService(){return fs;}

  // get input files (to be retrieved from samples_cross_sections.txt)
  void getInputFiles(std::vector<wprime::InputFile> & inputFiles);

  void setApplyMETCorrection(bool flag){applyMETCorrection_ = flag;}
  void setHadronicMETCalculated(bool flag){hadronicMETcalculated_ = flag;}

  bool shouldApplyMETCorrection(){return applyMETCorrection_;}

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
  bool runningOnData() const;
  
  //Check if Run/Evt is in Debug list
  bool DebugEvent(edm::EventBase const& event) const;

  // integrated luminosity in pb^-1
  float getLumi_ipb(){return lumi_ipb;}

  static inline bool SameTrigger(const std::string & versionedName, const std::string & wildcardedName){
    return (wildcardedName.find("*") == std::string::npos) ? 
      !versionedName.compare(wildcardedName) : //No '*', early triggers
      !versionedName.compare(0, versionedName.size()-1, wildcardedName, 0, wildcardedName.size()-1);//This assumes v is under 10, need to fix
  }

  static bool PassTriggersCut(edm::EventBase const & event, std::string label,const std::vector<std::string>& triggerNames);
  static bool PassTriggersCut(const pat::TriggerEvent & triggerEvent,const std::vector<std::string>& triggerNames);
  static void PrintPassingTriggers(const pat::TriggerEvent & triggerEvent,const std::vector<std::string>& triggerNames);
  static bool FoundAndPassed(const pat::TriggerPathRef path, const std::vector<std::string>& triggerNames);
  static bool Passed(const pat::TriggerPathRef path);
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

  std::vector<edm::EventID> vEventsToDebug_;

  float lumi_ipb; // in pb^-1, to be retrieved from samples_cross_sections.txt

  // used for the parsing of samples_cross_sections.txt
  void parseLine(const std::string & new_line, wprime::InputFile * in_file);

  bool applyMETCorrection_; // wether to apply hadronic recoil correction (aimed for W)
  bool hadronicMETcalculated_; // want to calculate this max. once per event

  TVector2 hadronicMETcached;

  edm::InputTag genParticles_;
 // Handle to GenParticleCollectiom>
  edm::Handle<reco::GenParticleCollection> genParticles;
  
  edm::LumiReWeighting LumiWeights_;
  
};

#endif //#define _wprime_util_h_
