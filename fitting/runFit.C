{
  TString incpath = gSystem->GetIncludePath();
  string pwd = " -I\""; pwd += gSystem->pwd(); pwd += "\"";
  incpath.Append(pwd.c_str());
  gSystem->SetIncludePath(incpath.Data());

  // see common_fit.h for the mass points to which mass_option corresponds
  const unsigned mass_option = 2;
  const unsigned N_EXP = 1; // # of pseudo-experiments
  const bool bgdOnlyFit = false;

  // compile code
  gROOT->ProcessLine(".L Results.C+");
  gROOT->ProcessLine(".L fit_wprime.cpp+");
  gROOT->ProcessLine(".L fitSigBgd_eventLoop.C+");
  gROOT->ProcessLine(".L fitSigBgd.C+");
  gROOT->ProcessLine("fitSigBgd(mass_option, N_EXP, bgdOnlyFit)"); //

}
