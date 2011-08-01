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
  void SetEventsToDebug(const std::vector<edm::EventID>& vEvents){vEventsToDebug_ = vEvents;}

  std::string getSampleName() const{return samplename_;}
  float getWeight() const{return weight_;}
  float getLumiWeight   (const int   & nInt){ return LumiWeights_.weight   (nInt);}
  float getLumiWeight3BX(const float & nInt){ return LumiWeights_.weight3BX(nInt);}
  int   GetPU1BX(const std::vector< PileupSummaryInfo > & PupInfo);
  float GetPUWeight1BX(const std::vector< PileupSummaryInfo > & PupInfo);
  float GetPU3BX(const std::vector< PileupSummaryInfo > & PupInfo);
  float GetPUWeight3BX(const std::vector< PileupSummaryInfo > & PupInfo);

  static void getEff(float & eff, float & deff,float Num,float Denom);

  // true if current file under processing contains "data" in description
  bool runningOnData() const;
  
  //Check if Run/Evt is in Debug list
  bool DebugEvent(edm::EventBase const& event) const;

  // integrated luminosity in pb^-1
  float getLumi_ipb(){return lumi_ipb;}

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

  void SetLumiWeights(std::string & MCFile, std::string & DataFile, 
                      std::string & MCHist, std::string & DataHist);
  static void CheckStream(ofstream& stream, std::string & s);

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
