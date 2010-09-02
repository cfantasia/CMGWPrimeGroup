{
  TString incpath = gSystem->GetIncludePath();
  string pwd = " -I\""; pwd += gSystem->pwd(); pwd += "\"";
  incpath.Append(pwd.c_str());
  gSystem->SetIncludePath(incpath.Data());

  // compile code
  gROOT->ProcessLine(".L Results.C+");
  gROOT->ProcessLine(".L fit_wprime.cpp+");
  gROOT->ProcessLine(".L fitSigBgd_eventLoop.C+");
  gROOT->ProcessLine(".L fitSigBgd_data.C+");
  gROOT->ProcessLine("fitSigBgd_data()"); //

}
