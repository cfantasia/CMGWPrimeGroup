{
  TString incpath = gSystem->GetIncludePath();
  string pwd = " -I\""; pwd += gSystem->pwd(); pwd += "\"";
  incpath.Append(pwd.c_str());
  gSystem->SetIncludePath(incpath.Data());

  // valid mass_option parameters:
  // 0->signal-free, 1->0.8 TeV, 2->1.0 TeV, 3->1.1TeV, 4->1.2TeV, 5->1.3TeV, 
  // 6->1.4TeV, 7->1.5TeV, 8->2.0 TeV
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
