/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooDLLSignificanceMCSModule2.h,v 1.1 2012/04/16 14:41:08 cleonido Exp $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *                                                                           *
 * Modified by Phil Hebda to add chi2,etc distributions for H0 hypothesis fit*
 *  (April 2012)                                                             *

 *****************************************************************************/


#ifndef ROO_DELTA_LL_SIGNIFICANCE_MCS_MODULE2
#define ROO_DELTA_LL_SIGNIFICANCE_MCS_MODULE2

#include "RooAbsMCStudyModule.h"
#include <string>

class RooDLLSignificanceMCSModule2 : public RooAbsMCStudyModule {
 public:

  RooDLLSignificanceMCSModule2(const RooRealVar& param, Double_t nullHypoValue=0) ;
  RooDLLSignificanceMCSModule2(const char* parName, Double_t nullHypoValue=0) ;
  RooDLLSignificanceMCSModule2(const RooDLLSignificanceMCSModule2& other) ;
  virtual ~RooDLLSignificanceMCSModule2() ;

  Bool_t initializeInstance() ; 

  Bool_t initializeRun(Int_t /*numSamples*/) ; 
  RooDataSet* finalizeRun() ;

  Bool_t processAfterFit(Int_t /*sampleNum*/)  ;
  
 private:

  std::string _parName ;  // Name of Nsignal parameter
  RooDataSet* _data ;     // Summary dataset to store results
  RooRealVar* _nll0h ;    // Container variable for NLL result on null hypothesis
  RooRealVar* _nll1h ;    // Container variable for NLL result on H1 hypothesis
  RooRealVar* _dllh ;    // Container variable for delta NLL 
  RooRealVar* _sig0h ;    // Container variable for NLL result with signal
  Double_t    _nullValue ;  // Numeric value of Nsignal parameter representing the null hypothesis

  /*
  RooRealVar* _chi2H0 ;    // Chi^2 of function w.r.t. data for H0
  RooRealVar* _ndofH0 ;    // Number of degrees of freedom for H0
  RooRealVar* _chi2redH0 ; // Reduced Chi^2 w.r.t data for H0
  RooRealVar* _probH0 ;    // Probability of chi^2,nDOF combination for H0
  RooRealVar* _chi2H1 ;    // Chi^2 of function w.r.t. data for H1
  RooRealVar* _ndofH1 ;    // Number of degrees of freedom for H1
  RooRealVar* _chi2redH1 ; // Reduced Chi^2 w.r.t data for H1
  RooRealVar* _probH1 ;    // Probability of chi^2,nDOF combination for H1
  */
  RooFitResult* frnull ;      // Fit results for H0

  ClassDef(RooDLLSignificanceMCSModule2,1) // MCStudy module to calculate Delta(-logL) significance w.r.t given null hypothesis
    } ;

#endif
