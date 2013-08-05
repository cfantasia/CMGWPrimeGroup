//Usage: root -b -q 'OptimizeCuts.C+()'
#include "consts.h"

void
OptimizeCuts(){
//Read in limits file and find best one for each mass
  
  TTree* tLims = new TTree("tLims", "Lims");
  tLims->ReadFile(Form("nLimit_%i.txt", mass));
                  
  TTree* tCuts = new TTree("tCuts", "Cuts");
  tCuts->ReadFile("cutValues.wzFull.dat");
  
  tLims->Draw("ExpLimit", Form("SignalCode==%f",code), "para goff");//Cpry put in singal code check with maps
  float nLims = tLims->GetSelectedRows(); 

  tCuts->Draw("HtCut:minWindow:maxWindow", Form("SignalCode==%f",code), "para goff");
  int nCuts = tCuts->GetSelectedRows();

  assert(nLims == nCuts);

  Double_t  *limits    = tLims->GetVal(0);
  Double_t  *minHt     = tCuts->GetVal(0);
  Double_t  *minWind   = tCuts->GetVal(1);
  Double_t  *maxWind   = tCuts->GetVal(2);

  Double_t *windSize   = new Double_t[nCuts];
  for(int i=0; i< nCuts; ++i){
    windSize[i] = maxWind[i] - minWind[i];
  }

  TGraph2D* lim2D = new TGraph2D(nLims, minHt, windSize, limits);
  
  gStyle->SetPalette(1);
  lim2D->Draw("pcol");

  float xmin = lim2D->GetXmin();
  float xmax = lim2D->GetXmax();
  float ymin = lim2D->GetYmin();
  float ymax = lim2D->GetYmax();

  float zmin = lim2D->GetZmin();

  for(int x=xmin; x<=xmax; x+=10){
    for(int y=ymin; y<=ymax; y+=10){
      if(lim2D->Interpolate(x,y) == zmin){
        cout<<" Ht= "<<x
            <<" windSize = "<<y
            <<endl;
      }
    }
  }
  
/*
  int bestIdx = 0;
  for(int iLim=0; iLim<nLims; ++iLim){

    

//    if(bestIdx == 0) cout<<" For code "<<codes[iLim]<<" the best window was the smaller, consider trying even smaller\n";
//    if(bestIdx == nLims-1) cout<<" For code "<<Lims[iLim]<<" the best window was the biggest, consider trying even bigger\n";

    int iVar=0;
    const Double_t  HtCut     = tCuts->GetVal(iVar++)[bestIdx];
    const Double_t  minWind   = tCuts->GetVal(iVar++)[bestIdx];
    const Double_t  maxWind   = tCuts->GetVal(iVar++)[bestIdx];

  }
*/
}
