#ifndef _TBAnalyzer_h_
#define _TBAnalyzer_h_

#include "UserCode/CMGWPrimeGroup/interface/AnalyzerBase.h"

class TBAnalyzer : public AnalyzerBase {
public:
  TBAnalyzer();                         // constructor; initialize the list to be empty
  TBAnalyzer(const edm::ParameterSet & cfg, WPrimeUtil * wprimeUtil);
  ~TBAnalyzer();

  //methods for stuff to be once per job
  void setupCutOrder();

  //methods for stuff to be done for each sample
  void defineHistos(const TFileDirectory& dir);

  //methods for stuff to be done for each event
  void eventLoop(edm::EventBase const & event);
  bool passCuts(const float& weight=1.);
  void clearEvtVariables();
  void fillHistos(const int& index, const float& weight=1.);

  //////Cuts//////////
  bool passValidBCut() const;
  bool passValidTCut() const;
  bool passValidTBCut() const;

//////Chosen Candidates
  pat::Jet bCand1_;
  pat::Jet bCand2_;
  WCandidate wCand_;
  XWLeptonic tCand_;
  XWLeptonic tbCand_;

// +++++++++++++++++++ Histogram Definitions

  //Cuts 
  typedef bool (TBAnalyzer::*CutFnPtr)() const; 
  std::map<std::string,CutFnPtr> mFnPtrs_;
  std::vector<CutFnPtr> CutFns_;

};

#endif//#define _TBAnalyzer_h_
