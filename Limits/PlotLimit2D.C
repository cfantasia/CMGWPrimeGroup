//Author: Cory Fantasia 2011
//Purpose: Plot limit of of TC on 2D plane of MPi vs MRho
//Usage: root -b -l -q 'PlotLimit2D.C+()'

#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "../root_macros/CMSStyle.C"
#include "../root_macros/setTDR_modified.C"
#include "consts.h"
#include "TGraph2D.h"

const float minPI = 0;
const float maxPI = 1300;
const float maxRHO = 1500;
const float minRHO = 100;
const int RHO_INC = 10;
const int PI_INC = 10;

bool isInside(float x, float y);
int findIntersection(Double_t A[], Double_t key, int imin, int imax);

void PlotLimit2D() {
  //gErrorIgnoreLevel = kWarning;
  //CMSstyle();
  setTDRStyle();
  gROOT->ForceStyle();

  TCanvas* c1 = new TCanvas("TC-PiVsRho","");
  c1->DrawFrame(minRHO,minPI,maxRHO,maxPI);

  TTree* tLimit = new TTree("tLimit", "Limits");
  tLimit->ReadFile("../combined_limits/nLimit_WprimeWZ_MarkovChainMC.txt");
  tLimit->Draw("Mass:Lumi:ObsLimit:ExpLimit", "Mass%100==0", "para goff");
  int n = tLimit->GetSelectedRows(); assert(n); 
  const Double_t* mass = tLimit->GetVal(0);
  const Double_t* lumi = tLimit->GetVal(1);
  const Double_t* ObsLimit = tLimit->GetVal(2);
  const Double_t* ExpLimit = tLimit->GetVal(3);
  TGraph* gExpLim = new TGraph(n, mass, ExpLimit);
  TGraph* gObsLim = new TGraph(n, mass, ObsLimit);
  cout<<"Done with importing limits.  There were "<<n<<endl;

  //Fill Xsec 2D graph
  TTree* tXsec = new TTree("tXsec", "X Sec");
  tXsec->ReadFile("xSec_TCWZ.dat");
  tXsec->Draw("Rho:Pi:Xsec", "", "para goff");
  int nXsec = tXsec->GetSelectedRows(); assert(nXsec); 
  Double_t* gRho = tXsec->GetVal(0);
  Double_t* gPi =  tXsec->GetVal(1);
  Double_t* gXsec = tXsec->GetVal(2);
  TGraph2D *gXsec2D = new TGraph2D(nXsec, gRho, gPi, gXsec);
  cout<<"Done with 2D xsec graph. There were "<<nXsec<<endl;

  vector<float> masses, expPiLims, obsPiLims;
  masses.reserve((maxRHO-minRHO));
  expPiLims.reserve((maxRHO-minRHO));
  obsPiLims.reserve((maxRHO-minRHO));
  for(int rho=minRHO; rho<=maxRHO; rho+=RHO_INC){
    int startPI = 0;
    masses.push_back(rho);
    bool debug = false;
    if(masses.size() % 10 == 0 || rho%100==0) debug = true;
    if(rho > 300 && rho < 375) debug = true;

    if(debug) cout<<"Working on rho = "<<rho<<endl;
    
    int nRho = tXsec->Draw("Pi:Xsec", Form("Rho==%i",rho), "para goff");
    if(debug) cout<<"nRho is "<<nRho<<" for "<<Form("Rho==%i",rho)<<endl;
    TGraph*  gXsec1D;
    if(nRho>0) gXsec1D = new TGraph(nRho, tXsec->GetVal(0), tXsec->GetVal(1));
    else{      gXsec1D = new TGraph((maxPI-minPI)/PI_INC + 1);
      for(int pi=minPI, i=0; pi<=maxPI; pi+=PI_INC, ++i){
        gXsec1D->SetPoint(i, pi, gXsec2D->Interpolate(rho,pi));
      }
    }
    if(debug) cout<<"Done with 1D xsec"<<endl;

    float expXsecLim = max(0., gExpLim->Eval(rho));
    if(debug) cout<<"Exp lim is "<<expXsecLim<<endl;
    //float limExp = findIntersection(gXsec1D->GetX(), expXsecLim, 0, gXsec1D->GetN());
    float limExp = maxPI;
    startPI = minPI; //(gXsec1D->Eval(700) > expXsecLim) ? 700 : minPI;
    for(int pi=startPI; pi<=maxPI; pi+=PI_INC){
      //cout<<"old: "<<gXsec2D->Interpolate(rho,pi)<<" new: "<<gXsec1D->Eval(pi)<<endl;
      //if(pi%100 == 0) cout<<" Working Exp Lim on pi = "<<pi<<endl;
      //if(gXsec2D->Interpolate(rho,pi) > expXsecLim){
      if(gXsec1D->Eval(pi) > expXsecLim){
        limExp = pi;
        break;
      }
    }
    //printf("mine: %.4f old %.4f\n", 
    //       findIntersection(minPI, maxPI, rho, gXsec2D, expXsecLim), limExp);
    if(debug) cout<<" exp limit is "<<limExp<<endl;
    expPiLims.push_back(limExp);

    float obsXsecLim = max(0., gObsLim->Eval(rho));
    if(debug) cout<<"Obs lim is "<<obsXsecLim<<endl;
    //float limObs = findIntersection(minPI, maxPI, rho, gXsec2D, obsXsecLim);
    float limObs = maxPI;
    startPI = minPI; //(gXsec1D->Eval(700) > obsXsecLim) ? 700 : minPI;
    for(int pi=startPI; pi<=maxPI; pi+=PI_INC){
      //if(pi%100 == 0) cout<<" Working Obs Lim on pi = "<<pi<<endl;
      //if(gXsec2D->Interpolate(rho,pi) > obsXsecLim){
      if(gXsec1D->Eval(pi) > obsXsecLim){
        limObs = pi;
        break;
      }
    }
    if(debug) cout<<" obs limit is "<<limObs<<endl;
    obsPiLims.push_back(limObs);

    delete gXsec1D;
  }

  cout<<"Done with looping "<<endl;
  TGraph* exp = new TGraph(expPiLims.size()+1);//n+2 to close the curve
  TGraph* obs = new TGraph(expPiLims.size()+1);

  float minExpRho(-1), minObsRho(-1);
  for(unsigned i=0; i<expPiLims.size(); ++i){//Find min excluded
    if(minExpRho < 0 && expPiLims[i] < maxPI) minExpRho = masses[i];
    if(minObsRho < 0 && obsPiLims[i] < maxPI) minObsRho = masses[i];
  }
  float maxExpRho(-1), maxObsRho(-1);
  for(int i=expPiLims.size()-1; i>=0; --i){//Find max excluded
    if(maxExpRho < 0 && expPiLims[i] < maxPI) maxExpRho = masses[i];
    if(maxObsRho < 0 && obsPiLims[i] < maxPI) maxObsRho = masses[i];
  }

  for(unsigned i=0; i<expPiLims.size(); ++i){//Fill Graph
    float xExp = masses[i];
    //if(xExp < minExpRho) xExp = minExpRho;
    //if(xExp > maxExpRho) xExp = maxExpRho;
    //if(expPiLims[i] == maxPI) continue;
    exp->SetPoint(i+1, xExp, expPiLims[i]);
  }

  for(unsigned i=0; i<expPiLims.size(); ++i){//Fill Graph
    float xObs = masses[i];
    //if(xObs < minObsRho) xObs = minObsRho;
    //if(xObs > maxObsRho) xObs = maxObsRho;
    //if(xObs > 299. && xObs <= 301 ) obsPiLims[i] += 10.;//Cory: hack to smooth plot
    //if(obsPiLims[i] == maxPI) continue;
    obs->SetPoint(i+1, xObs, obsPiLims[i]);

    //cout<<"mass:expLim:obsLim:ExpPi:ObsPi = "<<masses[i]<<":"<<gExpLim->Eval(masses[i])<<":"<<gObsLim->Eval(masses[i])<<"\t"<<expPiLims[i]<<"\t\t"<<obsPiLims[i]<<endl;  
  }

  //Close the curve
  //exp->SetPoint(expPiLims.size()+1, maxExpRho, maxPI);
  //obs->SetPoint(expPiLims.size()+1, maxObsRho, maxPI);
  exp->SetPoint(0                 , minExpRho, maxPI);
  obs->SetPoint(0                 , minObsRho, maxPI);

  ///////Find Optimistic Limits///////////
  //Expected
  for(int i=expPiLims.size()-1; i>=0; --i){///////
    if(masses[i] > 350) continue;
    if(expPiLims[i] > masses[i] - 80.4){
      cout<<"Lower Optimistic Exp Limit is "<<masses[i]<<endl;
      break;
    }
  }
  for(unsigned i=0; i<expPiLims.size(); ++i){///////
    if(masses[i] < 350) continue;
    if(expPiLims[i] > masses[i] - 80.4){
      cout<<"Upper Optimistic Exp Limit is "<<masses[i]<<endl;
      break;
    }
  }
  //Observed
  for(int i=obsPiLims.size()-1; i>=0; --i){///////
    if(masses[i] > 350) continue;
    if(obsPiLims[i] > masses[i] - 80.4){
      cout<<"Lower Optimistic Obs Limit is "<<masses[i]<<endl;
      break;
    }
  }
  for(unsigned i=0; i<obsPiLims.size(); ++i){///////
    if(masses[i] < 350) continue;
    if(obsPiLims[i] > masses[i] - 80.4){
      cout<<"Upper Optimistic Obs Limit is "<<masses[i]<<endl;
      break;
    }
  }

  /////Les Houches Param/////////////
  //Expected
  for(int i=expPiLims.size()-1; i>=0; --i){///////
    if(masses[i] > 350) continue;
    if(expPiLims[i] > 0.75 * masses[i] - 25){
      cout<<"Lower Les Houches Exp Limit is "<<masses[i]<<endl;
      break;
    }
  }
  for(unsigned i=0; i<expPiLims.size(); ++i){///////
    if(masses[i] < 350) continue;
    if(expPiLims[i] > 0.75 * masses[i] - 25){
      cout<<"Upper Les Houches Exp Limit is "<<masses[i]<<endl;
      break;
    }
  }
  //Observed
  for(int i=obsPiLims.size()-1; i>=0; --i){///////
    if(masses[i] > 350) continue;
    if(obsPiLims[i] > 0.75 * masses[i] - 25){
      cout<<"Lower Les Houches Obs Limit is "<<masses[i]<<endl;
      break;
    }
  }
  for(unsigned i=0; i<obsPiLims.size(); ++i){///////
    if(masses[i] < 350) continue;
    if(obsPiLims[i] > 0.75 * masses[i] - 25){
      cout<<"Upper Les Houches Obs Limit is "<<masses[i]<<endl;
      break;
    }
  }



  /////CDF Bump Param/////////

  for(unsigned i=0; i<obsPiLims.size(); ++i){///////
    if(masses[i] == 290) cout<<"For rho="<<masses[i]<<", pi limit is "<<obsPiLims[i]
                             <<" (obs) and "<<expPiLims[i]<<" (exp)"<<endl;
  }

  //For PRL
  obs->SetFillStyle(1001);
  exp->SetFillStyle(0);
  obs->SetFillColor(kCyan-8);
  exp->SetFillColor(0);
  obs->SetLineColor(kCyan-8);
  exp->SetLineColor(kBlack);
  exp->SetLineStyle(3);
  exp->SetLineWidth(3);
/*//From FR
  obs->SetFillStyle(3003);
  exp->SetFillStyle(3004);
  obs->SetFillColor(kRed);
  exp->SetFillColor(4);
  obs->SetLineColor(kRed);
  exp->SetLineColor(4);
*/
  TMultiGraph *mg = new TMultiGraph("mg", ";M(#rho_{TC}) (GeV);M(#pi_{TC}) (GeV)");
  mg->Add(obs,"F");
  mg->Add(exp,"F");
  mg->Add(obs,"L");//Don't need line anymore
  mg->Add(exp,"L");
  mg->SetMinimum(minPI);
  mg->SetMaximum(maxPI);
  mg->Draw("");
  mg->GetXaxis()->SetNdivisions(505);
  mg->GetYaxis()->SetNdivisions(505);

  TGraph* cdfPoint = new TGraph(1);
  cdfPoint->SetMarkerColor(kBlack);
  cdfPoint->SetMarkerStyle(kFullStar);
  cdfPoint->SetMarkerSize(2);
  cdfPoint->SetPoint(0,290, 160);
  //cdfPoint->Draw("p");

  TLegend* leg = new TLegend(0.20, 0.75, 0.60, 0.90);
  leg->SetTextSize(0.05);
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(exp, "Exp. 95% C.L.", "f");
  leg->AddEntry(obs, "Obs. 95% C.L.", "f");
  //leg->AddEntry(cdfPoint, "CDF Anomoly", "p");
  leg->Draw();

  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextSize(0.05);
  latexLabel.SetTextFont(42);
  latexLabel.DrawLatex(0.33, 0.96, "CMS Preliminary 2012");
  latexLabel.DrawLatex(0.59, 0.30, "#sqrt{s} = 8 TeV");
  latexLabel.DrawLatex(0.57, 0.20, Form("#intL dt = %.1f fb^{-1}",lumi[0]/1000.));

  ////
  float x1, y1, x2, y2;
  if(isInside(minRHO,minRHO - 80.4)){
    x1 = minRHO; y1 = minRHO-80.4;
  }else{
    x1 = minPI+80.4; y1 = minPI;
  }
  if(isInside(maxRHO,maxRHO - 80.4)){
    x2 = maxRHO;     y2 = maxRHO-80.4;
  }else{
    x2 = maxPI+80.4; y2 = maxPI;
  }
  TLine* line1 = new TLine(x1, y1, x2, y2);
  line1->SetLineStyle(2);
  line1->SetLineWidth(2);
  line1->Draw();

  TLatex* text3 = new TLatex(317, 265, "M(#pi_{TC}) = M(#rho_{TC}) - M(W)");
  text3->SetTextSize(0.05);
  //text3->SetTextAngle(TMath::ATan((y2-y1)/(x2-x1))*TMath::RadToDeg());
  text3->SetTextAngle(50);
  text3->Draw(); 

  c1->SetTitle(";M(#rho_{TC}) (GeV);M(#pi_{TC}) (GeV)");

  
  //pi = .75 * rho - 25
  if(isInside(minRHO, .75*minRHO-25.)){
    x1 = minRHO; y1=.75*minRHO-25.;
  }else{
    x1 = 4./3.*(maxPI+25); y1 = maxPI;
  }
  if(isInside(maxRHO, .75*maxRHO-25.)){
    x2 = maxRHO;           y2=.75*maxRHO-25.;
  }else{
    x2 = 4./3.*(maxPI+25); y2 = maxPI;
  }
  TLine* line2 = new TLine(x1, y1, x2, y2);
  line2->SetLineStyle(1);
  line2->SetLineWidth(2);
  line2->Draw();

  TLatex* text4 = new TLatex(357, 265, "M(#pi_{TC}) = #frac{3}{4}M(#rho_{TC}) - 25 GeV");
  text4->SetTextSize(0.05);
  text4->SetTextAngle(40);
  text4->Draw(); 
  
  c1->RedrawAxis();
  c1->Print("tcLimit-2D.pdf");
  c1->Print("tcLimit-2D.C");
  //c1->Print("tcLimit-2D.png");

  TFile *fin = TFile::Open("WprimeWZ-Limits.root", "update"); assert(fin);
  c1->Write();
  fin->Close();

  /*
  cout<<"Extra points:\n";
  for(int i=150; i<=2000; i+=25){
    cout<<i<<" "<<i*3./4. - 25<<" 0.333333 "<<gXsec2D->Interpolate(i, i*3./4. - 25)<<endl;
  }
  */
}


bool
isInside(float x, float y){
  return (minRHO <= x  && x <= maxRHO &&
          minPI  <= y  && y <= maxPI);
}

int
findIntersection(Double_t A[], Double_t key, int imin, int imax){
  // test if array is empty
  if (imax < imin){
    // set is empty, so return value showing not found
    return -1;//KEY_NOT_FOUND;
  }else{
    // calculate midpoint to cut set in half
    int imid = (imin+imax)/2;
    
    // three-way comparison
    if (A[imid] > key)
      // key is in lower subset
      return findIntersection(A, key, imin, imid-1);
    else if (A[imid] < key)
      // key is in upper subset
      return findIntersection(A, key, imid+1, imax);
    else
      // key has been found
      return imid;
  }
}
