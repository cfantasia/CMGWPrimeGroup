// MassPoint: should be between 1 and Nsignal_points-1
// channel: should be 0 (MuMET) or 1 (ElMET)
int Make_singleMassPoint(int MassPoint, int channel)
{
  /* find location of include file by issuing
     scramv1 tool info roofitcore
   */
  TString incpath = gSystem->GetIncludePath();
  incpath.Append(" -I/afs/cern.ch/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms10/include");
  gSystem->SetIncludePath(incpath.Data());
  

  gSystem->Load("libRooFit") ;
  using namespace RooFit ;

  /*
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/fitting/JacobianRBWPdf.cxx+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/fitting/RooBgdPdf.cxx+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/fitting/RooBgdPdf2.cxx+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/fitting/TripleGauss.cxx+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/fitting/WprimeFitter.cpp+");
    
  gROOT->ProcessLine(".x UserCode/CMGWPrimeGroup/fitting/fit_wprime.C+");

  */
  
  gROOT->ProcessLine(".L JacobianRBWPdf.cxx+");
  gROOT->ProcessLine(".L RooBgdPdf.cxx+");
  gROOT->ProcessLine(".L RooBgdPdf2.cxx+");
  gROOT->ProcessLine(".L TripleGauss.cxx+");
  gROOT->ProcessLine(".L WprimeFitter.cpp+");
  
  char command[1024];
  sprintf(command, ".x fit_wprime_singleMassPoint.C+(%d, %d)", MassPoint, channel);
  // cout << command << endl;
  gROOT->ProcessLine(command);
  return 0;
}
