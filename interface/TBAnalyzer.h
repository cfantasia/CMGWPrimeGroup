#ifndef _TBAnalyzer_h_
#define _TBAnalyzer_h_

#include "UserCode/CMGWPrimeGroup/interface/AnalyzerBase.h"

class TBAnalyzer : public AnalyzerBase {
public:
  TBAnalyzer();                         // constructor; initialize the list to be empty
  TBAnalyzer(const edm::ParameterSet & cfg, int fileToRun);
  ~TBAnalyzer();

  //methods for stuff to be once per job
  void setupCutOrder();

  //methods for stuff to be done for each sample
  void defineHistos(const TFileDirectory& dir);

  //methods for stuff to be done for each event
  void eventLoop(edm::EventBase const & event);
  void clearEvtVariables();
  void fillHistos(const int& index, const float& weight=1.);

  //////Cuts//////////
  bool passValidBCut(const pat::Jet & b) const;
  bool passValidTCut(XWLeptonic & t) const;
  bool passValidTBCut(XWLeptonic & tb) const;

//////Chosen Candidates
  pat::Jet bCand1_;
  pat::Jet bCand2_;
  WCandidate wCand_;
  XWLeptonic tCand_;
  XWLeptonic tbCand_;

// +++++++++++++++++++ Histogram Definitions

};

#endif//#define _TBAnalyzer_h_
