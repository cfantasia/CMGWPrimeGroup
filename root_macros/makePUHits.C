#include <iostream>
#include "TFile.h"
#include "TH1.h"

void
makePUHits(){
  // Run this script with root -b -l -q makePUHits.C+

  TFile MCFile("MCPUDist.root", "recreate");

  // Distribution used for Summer2012 MC.
  Double_t Summer2012[60] = {
    2.344E-05,
    2.344E-05,
    2.344E-05,
    2.344E-05,
    4.687E-04,
    4.687E-04,
    7.032E-04,
    9.414E-04,
    1.234E-03,
    1.603E-03,
    2.464E-03,
    3.250E-03,
    5.021E-03,
    6.644E-03,
    8.502E-03,
    1.121E-02,
    1.518E-02,
    2.033E-02,
    2.608E-02,
    3.171E-02,
    3.667E-02,
    4.060E-02,
    4.338E-02,
    4.520E-02,
    4.641E-02,
    4.735E-02,
    4.816E-02,
    4.881E-02,
    4.917E-02,
    4.909E-02,
    4.842E-02,
    4.707E-02,
    4.501E-02,
    4.228E-02,
    3.896E-02,
    3.521E-02,
    3.118E-02,
    2.702E-02,
    2.287E-02,
    1.885E-02,
    1.508E-02,
    1.166E-02,
    8.673E-03,
    6.190E-03,
    4.222E-03,
    2.746E-03,
    1.698E-03,
    9.971E-04,
    5.549E-04,
    2.924E-04,
    1.457E-04,
    6.864E-05,
    3.054E-05,
    1.282E-05,
    5.081E-06,
    1.898E-06,
    6.688E-07,
    2.221E-07,
    6.947E-08,
    2.047E-08
  };  

  TH1D hSummer12Dist("Summer12Dist", "Summer 12 True PU Dist", 60, 0.5, 60.5);
  for(int i=0; i<60; ++i) hSummer12Dist.Fill(i+1, Summer2012[i]);

  //Intended distribution
//   Double_t probdistFlat10[25] = {
//     0.0698146584,
//     0.0698146584,
//     0.0698146584,
//     0.0698146584,
//     0.0698146584,
//     0.0698146584,
//     0.0698146584,
//     0.0698146584,
//     0.0698146584,
//     0.0698146584,
//     0.0698146584,
//     0.0630151648,
//     0.0526654164,
//     0.0402754482,
//     0.0292988928,
//     0.0194384503,
//     0.0122016783,
//     0.007207042,
//     0.004003637,
//     0.0020278322,
//     0.0010739954,
//     0.0004595759,
//     0.0002229748,
//     0.0001028162,
//     4.58337152809607E-05,
//   };
//   TH1D hprobdistFlat10("probdistFlat10", "probdistFlat10", 25, -0.5, 24.5);
//   for(int i=0; i<25; ++i) hprobdistFlat10.Fill(i, probdistFlat10[i]);

//   //Average of +1, 0, -1 BX
//   Double_t PoissonIntDist[25] = {
//     0.104109,
//     0.0703573,
//     0.0698445,
//     0.0698254,
//     0.0697054,
//     0.0697907,
//     0.0696751,
//     0.0694486,
//     0.0680332,
//     0.0651044,
//     0.0598036,
//     0.0527395,
//     0.0439513,
//     0.0352202,
//     0.0266714,
//     0.019411,
//     0.0133974,
//     0.00898536,
//     0.0057516,
//     0.00351493,
//     0.00212087,
//     0.00122891,
//     0.00070592,
//     0.000384744,
//     0.000219377
//   };
//   TH1D hPoissonIntDist("PoissonIntDist", "PoissonIntDist", 25, -0.5, 24.5);
//   for(int i=0; i<25; ++i) hPoissonIntDist.Fill(i, PoissonIntDist[i]);
  
//   // Summer11 PU_S4, distribution obtained by only looking at the in-time crossing.  This is the "spike+smear" distribution, RECOMMENDED FOR REWEIGHTING.
//   Double_t PoissonOneXDist[35] = {
//     1.45346E-01,
//     6.42802E-02,
//     6.95255E-02,
//     6.96747E-02,
//     6.92955E-02,
//     6.84997E-02,
//     6.69528E-02,
//     6.45515E-02,
//     6.09865E-02,
//     5.63323E-02,
//     5.07322E-02,
//     4.44681E-02,
//     3.79205E-02,
//     3.15131E-02,
//     2.54220E-02,
//     2.00184E-02,
//     1.53776E-02,
//     1.15387E-02,
//     8.47608E-03,
//     6.08715E-03,
//     4.28255E-03,
//     2.97185E-03,
//     2.01918E-03,
//     1.34490E-03,
//     8.81587E-04,
//     5.69954E-04,
//     3.61493E-04,
//     2.28692E-04,
//     1.40791E-04,
//     8.44606E-05,
//     5.10204E-05,
//     3.07802E-05,
//     1.81401E-05,
//     1.00201E-05,
//     5.80004E-06
//   };
// /*
//   //Only in-time pileup
//   Double_t PoissonOneXDist[25] = {
//     0.14551,
//     0.0644453,
//     0.0696412,
//     0.0700311,
//     0.0694257,
//     0.0685655,
//     0.0670929,
//     0.0646049,
//     0.0609383,
//     0.0564597,
//     0.0508014,
//     0.0445226,
//     0.0378796,
//     0.0314746,
//     0.0254139,
//     0.0200091,
//     0.0154191,
//     0.0116242,
//     0.00846857,
//     0.00614328,
//     0.00426355,
//     0.00300632,
//     0.00203485,
//     0.00133045,
//     0.000893794
//   };
// */
//   TH1D hPoissonOneXDist("PoissonOneXDist", "PoissonOneXDist", 35, -0.5, 34.5);
//   for(int i=0; i<35; ++i) hPoissonOneXDist.Fill(i, PoissonOneXDist[i]);


//   //////////////////
//   ////Fall 11///////
//   //////////////////

//   Double_t Fall2011[50] = {
//     0.003388501,
//     0.010357558,
//     0.024724258,
//     0.042348605,
//     0.058279812,
//     0.068851751,
//     0.072914824,
//     0.071579609,
//     0.066811668,
//     0.060672356,
//     0.054528356,
//     0.04919354,
//     0.044886042,
//     0.041341896,
//     0.0384679,
//     0.035871463,
//     0.03341952,
//     0.030915649,
//     0.028395374,
//     0.025798107,
//     0.023237445,
//     0.020602754,
//     0.0180688,
//     0.015559693,
//     0.013211063,
//     0.010964293,
//     0.008920993,
//     0.007080504,
//     0.005499239,
//     0.004187022,
//     0.003096474,
//     0.002237361,
//     0.001566428,
//     0.001074149,
//     0.000721755,
//     0.000470838,
//     0.00030268,
//     0.000184665,
//     0.000112883,
//     6.74043E-05,
//     3.82178E-05,
//     2.22847E-05,
//     1.20933E-05,
//     6.96173E-06,
//     3.4689E-06,
//     1.96172E-06,
//     8.49283E-07,
//     5.02393E-07,
//     2.15311E-07,
//     9.56938E-08
//   };

//   TH1D hFall11Dist("Fall11Dist", "Fall 11 Dist", 50, -0.5, 49.5);
//   for(int i=0; i<50; ++i) hFall11Dist.Fill(i, Fall2011[i]);

  MCFile.Write();
  MCFile.Close();

}

/*
hadd -f Cert_160404-180252_7TeV_Collisions11_JSON.pileupTruth_v2_finebin.root \
/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.pileupTruth_v2_finebin.root \
/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/Cert_165088-167913_7TeV_PromptReco_JSON.pileupTruth_v2_finebin.root \
/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v2.pileupTruth_v2_finebin.root \
/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/Cert_172620-173692_PromptReco_JSON.pileupTruth_v2_finebin.root \
/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/Cert_175832-177515_PromptReco_JSON.pileupTruth_v2_finebin.root \
/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/Cert_177718_178078_7TeV_PromptReco_Collisons11_JSON.pileupTruth_v2_finebin.root \
/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/Cert_178098-180252_7TeV_PromptReco_Collisions11_JSON.pileupTruth_v2_finebin.root
*/
