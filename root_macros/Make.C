{
  // Update the include path so that we can find wprimeEvent.cc
  TString incpath = gSystem->GetIncludePath();
  incpath.Append(" -I$CMSSW_BASE/src");
  gSystem->SetIncludePath(incpath.Data());

  // compile code
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/src/wprimeEvent.cc+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/read_wprime.C+");
  gROOT->ProcessLine("read_wprime()");
}
