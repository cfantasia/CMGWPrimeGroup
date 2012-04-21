/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id: RooDLLSignificanceMCSModule2.cxx,v 1.2 2012/04/19 15:33:29 cleonido Exp $
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

//////////////////////////////////////////////////////////////////////////////
// 
// BEGIN_HTML
// RooDLLSignificanceMCSModule is an add-on modules to RooMCStudy that
// calculates the significance of a signal by comparing the likelihood of
// a fit fit with a given parameter floating with a fit with that given
// parameter fixed to a nominal value (usually zero). The difference in
// the -log(L) of those two fits can be interpreted as the probability
// that a statistical background fluctation may result in a signal as large
// or larger than the signal observed. This interpretation is contingent
// on underlying normal sampling distributions and a MC study is a good way
// to test that assumption.
// END_HTML
//
//
#include "Riostream.h"

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "TString.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooDLLSignificanceMCSModule2.h"
#include "RooMsgService.h"
#include "RooChi2Var.h"
#include "RooDataHist.h"
#include "TMath.h"
#include "RooGlobalFunc.h"


ClassImp(RooDLLSignificanceMCSModule2)
  ;



//_____________________________________________________________________________
RooDLLSignificanceMCSModule2::RooDLLSignificanceMCSModule2(const RooRealVar& param, Double_t nullHypoValue) : 
  RooAbsMCStudyModule(Form("RooDLLSignificanceMCSModule2_%s",param.GetName()),Form("RooDLLSignificanceMCSModule2_%s",param.GetName())),
  _parName(param.GetName()), 
  _data(0), _nll0h(0), _nll1h(0), _dllh(0), _sig0h(0), _nullValue(nullHypoValue)
  //  _chi2H0(0), _ndofH0(0), _chi2redH0(0), _probH0(0),
  // _chi2H1(0), _ndofH1(0), _chi2redH1(0), _probH1(0)
{
  // Constructor of module with parameter to be interpreted as nSignal and the value of the
  // null hypothesis for nSignal (usually zero)
}



//_____________________________________________________________________________
RooDLLSignificanceMCSModule2::RooDLLSignificanceMCSModule2(const char* parName, Double_t nullHypoValue) :
  RooAbsMCStudyModule(Form("RooDLLSignificanceMCSModule2_%s",parName),Form("RooDLLSignificanceMCSModule2_%s",parName)),
  _parName(parName), 
  _data(0), _nll0h(0), _nll1h(0), _dllh(0), _sig0h(0), _nullValue(nullHypoValue)
  //  _chi2H0(0), _ndofH0(0), _chi2redH0(0), _probH0(0),
  // _chi2H1(0), _ndofH1(0), _chi2redH1(0), _probH1(0)
{
  // Constructor of module with parameter name to be interpreted as nSignal and the value of the
  // null hypothesis for nSignal (usually zero)
}



//_____________________________________________________________________________
RooDLLSignificanceMCSModule2::RooDLLSignificanceMCSModule2(const RooDLLSignificanceMCSModule2& other) : 
  RooAbsMCStudyModule(other), 
  _parName(other._parName),
  _data(0), _nll0h(0), _nll1h(0), _dllh(0), _sig0h(0), _nullValue(other._nullValue)
  //  _chi2H0(0), _ndofH0(0), _chi2redH0(0), _probH0(0),
  // _chi2H1(0), _ndofH1(0), _chi2redH1(0), _probH1(0)
{
  // Copy constructor
}



//_____________________________________________________________________________
RooDLLSignificanceMCSModule2:: ~RooDLLSignificanceMCSModule2() 
{
  // Destructor

  if (_nll0h) {
    delete _nll0h ;
  }    
  if (_nll1h) {
    delete _nll1h ;
  }    
  if (_dllh) {
    delete _dllh ;
  }    
  if (_sig0h) {
    delete _sig0h ;
  }
  if (_data) {
    delete _data ;
  }
  /*
  if (_chi2H0) {
    delete _chi2H0 ;
  }
  if (_ndofH0) {
    delete _ndofH0 ;
  }
  if (_chi2redH0) {
    delete _chi2redH0 ;
  }
  if (_probH0) {
    delete _probH0 ;
  }
  if (_chi2H1) {
    delete _chi2H1 ;
  }
  if (_ndofH1) {
    delete _ndofH1 ;
  }
  if (_chi2redH1) {
    delete _chi2redH1 ;
  }
  if (_probH1) {
    delete _probH1 ;
  }
  */
}



//_____________________________________________________________________________
Bool_t RooDLLSignificanceMCSModule2::initializeInstance()
{
  // Initialize module after attachment to RooMCStudy object

  // Check that parameter is also present in fit parameter list of RooMCStudy object
  if (!fitParams()->find(_parName.c_str())) {
    coutE(InputArguments) << "RooDLLSignificanceMCSModule2::initializeInstance:: ERROR: No parameter named " << _parName << " in RooMCStudy!" << endl ;
    return kFALSE ;
  }

  // Construct variable that holds -log(L) fit with null hypothesis for given parameter
  TString nll0hName = Form("nll_nullhypo_%s",_parName.c_str()) ;
  TString nll0hTitle = Form("-log(L) with null hypothesis for param %s",_parName.c_str()) ;
  _nll0h = new RooRealVar(nll0hName.Data(),nll0hTitle.Data(),0) ;
  TString nll1hName = Form("nll_H1_%s",_parName.c_str()) ;
  TString nll1hTitle = Form("-log(L) with H1 hypothesis for param %s",_parName.c_str()) ;
  _nll1h = new RooRealVar(nll1hName.Data(),nll1hTitle.Data(),0) ;

  // Construct variable that holds -log(L) fit with null hypothesis for given parameter
  TString dllhName = Form("dll_nullhypo_%s",_parName.c_str()) ;
  TString dllhTitle = Form("-log(L) difference w.r.t null hypo for param %s",_parName.c_str()) ;
  _dllh = new RooRealVar(dllhName.Data(),dllhTitle.Data(),0) ;

  // Construct variable that holds significance corresponding to delta(-log(L)) w.r.t to null hypothesis for given parameter
  TString sig0hName = Form("significance_nullhypo_%s",_parName.c_str()) ;
  TString sig0hTitle = Form("Gaussian signficiance of Delta(-log(L)) w.r.t null hypo for param %s",_parName.c_str()) ;
  _sig0h = new RooRealVar(sig0hName.Data(),sig0hTitle.Data(),-10,100) ;

  _nbgdH0 = new RooRealVar("nbgd_H0","nbgd from fit under H0",0);

  /*
  _chi2H0     = new RooRealVar("chi2H0","chi^2 for H0",0) ;
  _ndofH0     = new RooRealVar("ndofH0","number of degrees of freedom for H0",0) ;   
  _chi2redH0  = new RooRealVar("chi2redH0","reduced chi^2 for H0",0) ; 
  _probH0     = new RooRealVar("probH0","prob(chi2H0,ndofH0)",0) ;
  _chi2H1     = new RooRealVar("chi2H1","chi^2 for H1",0) ;
  _ndofH1     = new RooRealVar("ndofH1","number of degrees of freedom for H1",0) ;   
  _chi2redH1  = new RooRealVar("chi2redH1","reduced chi^2 for H1",0) ; 
  _probH1     = new RooRealVar("probH1","prob(chi2H1,ndofH1)",0) ;
  */

  // Construct variable that holds -log(L) fit with null hypothesis for given parameter
  frnull      = new RooFitResult("frnull","Fit result for H0");

  // Create new dataset to be merged with RooMCStudy::fitParDataSet
  _data = new RooDataSet("DeltaLLSigData","Additional data for Delta(-log(L)) study",RooArgSet(*_nll0h,*_nll1h, *_dllh,*_sig0h,*_nbgdH0)) ;
  //did not include *_chi2H0,*_ndofH0,*_chi2H1,*_ndofH1,*frnull due to space constraints with RooArgSet
  return kTRUE ;
}



//_____________________________________________________________________________
Bool_t RooDLLSignificanceMCSModule2::initializeRun(Int_t /*numSamples*/) 
{
  // Initialize module at beginning of RooCMStudy run

  _data->reset() ;
  return kTRUE ;
}



//_____________________________________________________________________________
RooDataSet* RooDLLSignificanceMCSModule2::finalizeRun() 
{
  // Return auxiliary dataset with results of delta(-log(L))
  // calculations of this module so that it is merged with
  // RooMCStudy::fitParDataSet() by RooMCStudy

  return _data ;
}



//_____________________________________________________________________________
Bool_t RooDLLSignificanceMCSModule2::processAfterFit(Int_t /*sampleNum*/)  
{
  // Bin dataset and calculate chi2 of p.d.f. w.r.t. binned dataset for H1

  RooAbsData* data = genSample() ;
  RooDataHist* binnedData = dynamic_cast<RooDataHist*>(data) ;
  Bool_t deleteData(kFALSE) ;
  if (!binnedData) {
    deleteData = kTRUE ;
    binnedData = ((RooDataSet*)data)->binnedClone() ;
  }

  /*
  RooChi2Var chi2VarH1("chi2VarH1","chi2VarH1",*fitModel(),*binnedData,RooFit::Extended(extendedGen()),RooFit::DataError(RooAbsData::SumW2)) ;

  RooArgSet* floatParsH1 = (RooArgSet*) fitParams()->selectByAttrib("Constant",kFALSE) ;  

  _chi2H1->setVal(chi2VarH1.getVal()) ;
  _ndofH1->setVal(binnedData->numEntries()-floatParsH1->getSize()-1) ; 
  _chi2redH1->setVal(_chi2H1->getVal()/_ndofH1->getVal()) ;
  _probH1->setVal(TMath::Prob(_chi2H1->getVal(),static_cast<int>(_ndofH1->getVal()))) ;
  */

  // Save likelihood from nominal fit, fix chosen parameter to its
  // null hypothesis value and rerun fit Save difference in likelihood
  // and associated Gaussian significance in auxilary dataset

  RooRealVar* par = static_cast<RooRealVar*>(fitParams()->find(_parName.c_str())) ;
  par->setVal(_nullValue) ;
  par->setConstant(kTRUE) ;
  frnull = refit() ;

  //Now repeat chi2 calculatation for H0 before continuing the likelihood calculations
  /*
  RooChi2Var chi2VarH0("chi2VarH0","chi2VarH0",*fitModel(),*binnedData,RooFit::Extended(extendedGen()),RooFit::DataError(RooAbsData::SumW2)) ;

  RooArgSet* floatParsH0 = (RooArgSet*) fitParams()->selectByAttrib("Constant",kFALSE) ;  

  _chi2H0->setVal(chi2VarH0.getVal()) ;
  _ndofH0->setVal(binnedData->numEntries()-floatParsH0->getSize()-1) ; 
  _chi2redH0->setVal(_chi2H0->getVal()/_ndofH0->getVal()) ;
  _probH0->setVal(TMath::Prob(_chi2H0->getVal(),static_cast<int>(_ndofH0->getVal()))) ;
  */

  if (deleteData) {
    delete binnedData ;
  }
  //  delete floatParsH1 ;
  // delete floatParsH0 ;

  //continue likelihood step
  par->setConstant(kFALSE) ;
  
  _nll0h->setVal(frnull->minNll()) ;
  _nll1h->setVal(nllVar()->getVal());
  Double_t deltaLL = (frnull->minNll() - nllVar()->getVal()) ;
  Double_t signif = deltaLL>0 ? sqrt(2*deltaLL) : -sqrt(-2*deltaLL) ;
  _sig0h->setVal(signif) ;
  _dllh->setVal(deltaLL) ;
  _nbgdH0->setVal( ((RooRealVar*)(frnull->floatParsFinal().at(frnull->floatParsFinal().index("nbgd"))))->getVal() ) ;
  //_nbgdH0->setVal(0) ;


  _data->add(RooArgSet(*_nll0h,*_nll1h,*_dllh,*_sig0h,*_nbgdH0)) ;

  return kTRUE ;
}
