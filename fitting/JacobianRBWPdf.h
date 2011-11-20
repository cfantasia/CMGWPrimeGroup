/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef JACOBIANRBWPDF
#define JACOBIANRBWPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class JacobianRBWPdf : public RooAbsPdf {
public:
  JacobianRBWPdf() {} ; 
  JacobianRBWPdf(const char *name, const char *title,
	      RooAbsReal& _mt,
	      RooAbsReal& _mass,
	      RooAbsReal& _width);
  JacobianRBWPdf(const JacobianRBWPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new JacobianRBWPdf(*this,newname); }
  inline virtual ~JacobianRBWPdf() { }

  static float fastInvSqrt(float x);

protected:

  RooRealProxy mt ;
  RooRealProxy mass ;
  RooRealProxy width ;
  
  Double_t evaluate() const ;

private:

  ClassDef(JacobianRBWPdf,1) // Your description goes here...
};
 
#endif
