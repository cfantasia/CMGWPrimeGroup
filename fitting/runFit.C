{
  TString incpath = gSystem->GetIncludePath();
  string pwd = " -I\""; pwd += gSystem->pwd(); pwd += "\"";
  incpath.Append(pwd.c_str());
  gSystem->SetIncludePath(incpath.Data());

  const unsigned mass_option = 1; //  1: 1.0 TeV, 2: 1.5 TeV
  const unsigned N_EXP = 1; // # of pseudo-experiments

  // compile code
  gROOT->ProcessLine(".L Results.C+");
  gROOT->ProcessLine(".L fit_wprime.cpp+");
  gROOT->ProcessLine(".L fitSigBgd_eventLoop.C+");
  gROOT->ProcessLine(".L fitSigBgd.C+");
  gROOT->ProcessLine("fitSigBgd(mass_option, N_EXP)"); //

}
