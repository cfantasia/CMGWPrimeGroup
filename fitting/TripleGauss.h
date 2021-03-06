/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef TRIPLEGAUSS
#define TRIPLEGAUSS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class TripleGauss : public RooAbsPdf {
public:
  TripleGauss() {} ; 
  TripleGauss(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _m1,
	      RooAbsReal& _s1,
	      RooAbsReal& _f1,
	      RooAbsReal& _m2,
	      RooAbsReal& _s2,
	      RooAbsReal& _f2,
	      RooAbsReal& _m3,
	      RooAbsReal& _s3);
  TripleGauss(const TripleGauss& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new TripleGauss(*this,newname); }
  inline virtual ~TripleGauss() { }

protected:

  RooRealProxy x ;
  RooRealProxy  m1 ;
  RooRealProxy  s1 ;
  RooRealProxy  f1 ;
  RooRealProxy  m2 ;
  RooRealProxy  s2 ;
  RooRealProxy  f2 ;
  RooRealProxy  m3 ;
  RooRealProxy s3 ;
  
  Double_t evaluate() const ;

private:

  ClassDef(TripleGauss,1) // Your description goes here...
};
 
#endif
