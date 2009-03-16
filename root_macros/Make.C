{

  gROOT->Reset();

  // Update the include path so that we can find wprimeEvent.cc
  TString incpath = gSystem->GetIncludePath();
  incpath.Append(" -I$CMSSW_BASE/src");
  gSystem->SetIncludePath(incpath.Data());

  // compile code
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/src/wprimeEvent.cc+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/loadInputFiles.C+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/loadCrossSections.C+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/loadCuts.C+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/GetMuonPtDistribution.C+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/GetMuonPtDistribution_JetIso.C+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/GetDistributionGeneric.C+");

  float lumiPb = 100;

  TFile *fout = new TFile("Wprime_analysis_V54.root","recreate");

  vector<wprime::InputFile> qcd_files; vector<wprime::InputFile> z_files;
  vector<wprime::InputFile> w_files; vector<wprime::InputFile> top_files;
  vector<wprime::InputFile> wprime10_files;
  vector<wprime::InputFile> wprime15_files;
  vector<wprime::InputFile> wprime20_files;

  gROOT->ProcessLine("loadCrossSections(qcd_files, z_files, w_files,top_files, wprime10_files, wprime15_files, wprime20_files)");

  
  string dir = "QCD";
  gROOT->ProcessLine("loadInputFiles(dir, qcd_files, lumiPb)");
  gROOT->ProcessLine("GetDistributionGeneric(qcd_files, fout, dir)");
  
  dir = "Z";
  gROOT->ProcessLine("loadInputFiles(dir, z_files, lumiPb)");
  gROOT->ProcessLine("GetDistributionGeneric(z_files, fout, dir)");

  dir = "W";
  gROOT->ProcessLine("loadInputFiles(dir, w_files, lumiPb)");
  gROOT->ProcessLine("GetDistributionGeneric(w_files, fout, dir)");
  
  dir = "Top";
  gROOT->ProcessLine("loadInputFiles(dir, top_files, lumiPb)");
  gROOT->ProcessLine("GetDistributionGeneric(top_files, fout, dir)");  

  dir = "wprime10";
  gROOT->ProcessLine("loadInputFiles(dir, wprime10_files, lumiPb)");
  gROOT->ProcessLine("GetDistributionGeneric(wprime10_files, fout, dir)");  

  dir = "wprime15";
  gROOT->ProcessLine("loadInputFiles(dir, wprime15_files, lumiPb)");
  gROOT->ProcessLine("GetDistributionGeneric(wprime15_files, fout, dir)");  

  dir = "wprime20";
  gROOT->ProcessLine("loadInputFiles(dir, wprime20_files, lumiPb)");
  gROOT->ProcessLine("GetDistributionGeneric(wprime20_files, fout, dir)");  

  fout->Close();
}
