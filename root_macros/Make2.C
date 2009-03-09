{
  
gROOT->Reset();
  // Update the include path so that we can find wprimeEvent.cc

  TString incpath = gSystem->GetIncludePath();
  incpath.Append(" -I$CMSSW_BASE/src");
  gSystem->SetIncludePath(incpath.Data());

  // compile code
//  gROOT->ProcessLine(".L wprimeEvent.cc+");
gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/src/wprimeEvent.cc+");
gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/GetDistributionGeneric.C+");


TFile *fout = new TFile("test.root","recreate");
float lumiPb = 100;

cout<<"Processing QCD "<<endl<<endl; 
const int NfilesQCD = 9;
 string low[NfilesQCD]= {"100", "150", "200", "300", "400", "600", "800", "1200", "1600"};
 string high[NfilesQCD]={"150", "200", "300", "400", "600", "800","1200", "1600", "up"};
 
 string fileNameQCD[NfilesQCD];
 
 for(int i =0;i<NfilesQCD;i++){fileNameQCD[i]=string("UserCode/CMGWPrimeGroup/QCD_")+low[i]+string("_")+high[i]+string("_212_Ideal_Minv_ptGlobMu.root");}
 


float weightQCD[NfilesQCD]={
  131.11111,
  19.52941,
  5.45455,
  0.75949,
  0.54074,
  0.02129,
  0.00923,
  0.00210,
  0.00216
};
for(int i =0;i<NfilesQCD;i++){
  weightQCD[i]= weightQCD[i]*lumiPb/100;
}

string dir = "QCD";
gROOT->ProcessLine("GetDistributionGeneric( NfilesQCD, fileNameQCD, weightQCD, fout, dir)");

cout<<endl<<"Processing Z "<<endl<<endl; 

const int NfilesZ = 11;
string lowZ[NfilesZ]={"30", "110","200", "300", "400", "500", "600", "700", "800", "900", "1000"}


string fileNameZ[NfilesZ];

for(int i =0;i<NfilesZ;i++){fileNameZ[i]=string("UserCode/CMGWPrimeGroup/Z_212_Ideal_Minv_ptGlobalMu_pt")+lowZ[i]+string(".root");}

float weightZ[NfilesZ]={
7.329843,
0.463415,
0.120000,
0.026667,
0.008333,
0.003222,
0.001444,
0.000900,
0.000471,
0.000211,
0.000483
};

for(int i =0;i<NfilesZ;i++){
  weightZ[i]=weightZ[i]*lumiPb/100.;
}
dir = "Z";
gROOT->ProcessLine("GetDistributionGeneric( NfilesZ, fileNameZ, weightZ, fout, dir)");

cout<<endl<<"Processing W "<<endl<<endl; 
const int NfilesW = 22;
string lowW[NfilesW]={"0", "200", "250", "300", "350", "400", "450", "500", "600", "700", "800", "900", "1000", "1100", "1200", "1300", "1400", "1500", "1600", "1700", "1800", "1900"}
string highW[NfilesW]={"200", "250", "300", "350", "400", "450", "500", "600", "700", "800", "900", "1000", "1100", "1200", "1300", "1400", "1500", "1600", "1700", "1800", "1900","2000"}

string fileNameW[NfilesW];

for(int i =0;i<NfilesW;i++){fileNameW[i]=string("UserCode/CMGWPrimeGroup/Wmunu_219_Ideal_Minv_")+lowW[i]+string("_")+highW[i]+string(".root");}

float weightW[NfilesW]={
1.16800000000,
0.00086790000,
0.00032530000,
0.00014310000,
0.00006946000,
0.00003653000,
0.00002026000,
0.00003730000,
0.00001426800,
0.00000589800,
0.00000262400,
0.00000121600,
0.00000096833,
0.00000047267,
0.00000024017,
0.00000012053,
0.00000006143,
0.00000003196,
0.00000001644,
0.00000000851,
0.00000000435,
0.00000000221
};

dir = "W";
gROOT->ProcessLine("GetDistributionGeneric( NfilesW, fileNameW, weightW, fout, dir)");

cout<<endl<<"Processing Top "<<endl<<endl; 
const int NfilesTop = 1;
string fileNameTop[NfilesTop];
fileNameTop[0]="UserCode/CMGWPrimeGroup/TTbar_Minv_PtMuGlobal.root";

float weightTop[NfilesTop];
weightTop[0]= 0.173745*lumiPb/100;

dir = "Top";
gROOT->ProcessLine("GetDistributionGeneric( NfilesTop, fileNameTop, weightTop, fout, dir)");


}
