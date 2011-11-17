//Macro to estimate the non-zpeaking background in a WZ resonance
//Author: Edgar Carrera, Cory Fantasia with bits from root examples
#include "consts.h"
const double PI=3.14159265358979323846;

//----------------------------------------------------------
Double_t zfuncbreit(Double_t *x, Double_t *par) {
//-----------------------------------------------------------
   //[0]=half width
   //[1]=mean
   //[2]=norm
    return par[0]*par[2]/( PI * ((x[0]-par[1])*(x[0]-par[1]) + par[0]*par[0]) );
}

// Sum of background and peak function
//----------------------------------------------------------
Double_t fitFunctionAll(Double_t *x, Double_t *par) {
//----------------------------------------------------------
  return fline(x,par) + zfuncbreit(x,&par[2]);
}

//-----------------------------------------------------------
TF1* shape_nozpeakfunc(TH1F* myh){
//-----------------------------------------------------------
    
    TF1 *Pfunc = new TF1("Pfunc",fline,FitWind_low,FitWind_high,2);
        
    //Pfunc->SetParLimits(0,0.0001,100.);
    //Pfunc->SetParameter(0,4.24);
    myh->Fit(Pfunc,"MRLL");
    //myh->Fit(Pfunc,"MRELL");
    return Pfunc;
    
}

//-----------------------------------------------------------
TF1* shape_zpeakfunc(TH1F* myh ){
//-----------------------------------------------------------
    
    TF1 *BGfunc = new TF1("BGfunc",zfuncbreit,ZWind_low,ZWind_high,3);

    BGfunc->SetParameter(0,5);
    BGfunc->SetParameter(1,91);
    BGfunc->SetParameter(2,10);
    
    BGfunc->SetParLimits(1,85.,95.);

    myh->Fit(BGfunc,"MRELL");

    return BGfunc;
    
}

//-----------------------------------------------------------
TF1* shape_allfunc(TH1F* myh, Double_t polpars[]){
//-----------------------------------------------------------
    TF1* FitAllfunc;
    FitAllfunc = new TF1("FitAllfunc",fitFunctionAll,FitWind_low,FitWind_high,5);
    
    FitAllfunc->SetParameter(0,polpars[0]); FitAllfunc->SetParName(0,"Const Coef");
    FitAllfunc->SetParameter(1,polpars[1]); FitAllfunc->SetParName(1,"Linear Coef");
    FitAllfunc->SetParameter(2,5); FitAllfunc->SetParName(2,"Half Width");
    FitAllfunc->SetParameter(3,91); FitAllfunc->SetParName(3,"Z peak");
    FitAllfunc->SetParameter(4,10); FitAllfunc->SetParName(4,"Breit Norm");
    
    //cout<<"finished constructing FitAllfunc"<<endl;        
    myh->Fit(FitAllfunc,"MRELL");

    return FitAllfunc;

}

//----------------------------------------------------------
void SidebandSubtraction(bool useData=false) {
//----------------------------------------------------------
    //gROOT->ProcessLine(".L tdrstyle.C");
    //gROOT->ProcessLine("setTDRStyle()");
    gROOT->ProcessLine("gStyle->SetOptStat(0)");
    
    TFile *f = TFile::Open("../Wprime_analysis.root", "read");

    map<string, float> allshapesamples;
    if(!useData){
      //1.get the samples to shape the Z
      //cout<<"$$$$$$$$$$$$$$$$$  get the genuine-Z background shape"<<endl;
      map<string, float> zshapesamples;
//      zshapesamples["Wprime400"]=KFactor("Wprime400");
      zshapesamples["WZ3l"]=KFactor("WZ");
      zshapesamples["ZGamma2l"]=KFactor("ZGamma");
      zshapesamples["DYJets2l"]=KFactor("ZJets");
      zshapesamples["ZZ"]=KFactor("ZZ");
//      zshapesamples["ZZ4l"]=KFactor("ZZ4l");
      //TH1F* zpeakhist = get_sum_of_hists(f, zshapesamples, "hZMass_AllCuts", 1);  

      //cout<<"$$$$$$$$$$$$$$$$$  get the non-genuine-Z background shape"<<endl;
      //1.get the samples to shape the non-genuine Z background
      map<string, float> nozshapesamples;
      nozshapesamples["TTbar2l"]=KFactor("TTbar2l");
      nozshapesamples["WenuJets"]=KFactor("WenuJets");
      nozshapesamples["WmunuJets"]=KFactor("WmunuJets");
//      nozshapesamples["WlnuJetsMadgraph"]  =KFactor("WlnuJets");
      //TH1F* nozpeakhist = get_sum_of_hists(f, nozshapesamples, "hZMass_AllCuts", 1);

      allshapesamples.insert(zshapesamples.begin(),zshapesamples.end());
      allshapesamples.insert(nozshapesamples.begin(),nozshapesamples.end());
    }else{
      allshapesamples["Run2010"]=1.;
      allshapesamples["Run2011"]=1.;
    }
    //TH1F* hAll = get_sum_of_hists(f, allshapesamples, "hZMass_AllCuts", 1);//Cory: Original
    TH1F* hAll = get_sum_of_hists(f, allshapesamples, "hZMass_WLepIso", 1);

    ////Fit Background Outside Z Peak
    TCanvas* cBkg = new TCanvas("cBkg", "Background Fit");
    Double_t bkgFuncPars[2];
    TF1* fBkg = new TF1("fBkg", fline, FitWind_low, FitWind_high, 2);
    reject = true;
    hAll->Fit(fBkg, "MRELL", "e1");
    reject = false;
    fBkg->GetParameters(bkgFuncPars);
    cBkg->Print(Form("fitBkg_useData%i.eps", useData));

    /////Fit the Whole Thing
    cout<<"$$$$$$$$$$$$$$ fit signal + background"<<endl;
    //2.shape the function
    TCanvas *cGlobal = new TCanvas("cGlobal", "Global Fit");
    TF1* fAll = shape_allfunc(hAll,bkgFuncPars);
    fAll = new TF1("fAll",fitFunctionAll,FitWind_low,FitWind_high,5);
    fAll->SetParameter(0,bkgFuncPars[0]); fAll->SetParName(0,"Const Coef");
    fAll->SetParameter(1,bkgFuncPars[1]); fAll->SetParName(1,"Linear Coef");
    fAll->SetParameter(2,1.); fAll->SetParName(2,"Half Width");
    fAll->SetParameter(3,91.); fAll->SetParName(3,"Z peak");
    fAll->SetParameter(4,1.); fAll->SetParName(4,"Breit Norm");

    fAll->SetLineColor(kMagenta);
    fAll->SetLineWidth(2);
    fAll->SetNpx(500);
    hAll->Fit(fAll,"MRELL", "e1");
    cGlobal->Print(Form("fitGlobal_useData%i.eps", useData));

    //$$$$$$$$$$$$$$$$$$$Draw and extract the components
    cout<<"  Final Steps "<<endl;

    TCanvas *cComponent = new TCanvas("cComponent", "Component Fit");
    //hAll->Fit(fAll,"MRELL", "e1");
    hAll->Draw("e1");
    hAll->SetMinimum(0);
    hAll->SetMarkerStyle(1);
    
    //Get the final results of the fit
    fAll->Draw("same");

    fBkg->SetLineColor(kRed);
    fBkg->SetNpx(500);
    
   // write the fit results into the par array
    Double_t* allFuncPars = fAll->GetParameters();
    Double_t* allFuncErrs = fAll->GetParErrors ();

    fBkg->SetParameters(allFuncPars);
    fBkg->SetParErrors (allFuncErrs);
    fBkg->Draw("same");

    TF1 *fSignal = new TF1("fSignal",zfuncbreit,ZWind_low,ZWind_high,3);
    fSignal->SetLineColor(kBlue);
    fSignal->SetNpx(500);
    fSignal->SetParameters(&allFuncPars[2]);
    fSignal->SetParErrors (&allFuncErrs[2]);
    fSignal->Draw("same"); 

    // draw the legend
    TLegend *legend=new TLegend(0.6,0.45,0.88,0.65);
    legend->AddEntry(hAll,"All MC","lpe");
    legend->AddEntry(fBkg,"Background fit","l");
    legend->AddEntry(fSignal,"Signal fit (BW)","l");
    legend->AddEntry(fAll,"Global Fit (Bkg+Sig)","l");
    legend->Draw("same");
    legend->SetFillColor(0);
    legend->SetLineColor(0);
    legend->SetLineWidth(0);
    legend->SetBorderSize(0);

    cComponent->Print(Form("fitComponents_useData%i.eps", useData));

/*
    allfunc->Print();
    fBkg->Print();
    double N_bk       = fBkg->Integral     (ZWind_low,ZWind_high);
    double sigma_N_bk = fBkg->IntegralError(ZWind_low,ZWind_high);

    cout<<"Number of background events in range = "<<N_bk
        <<" +/- "<<sigma_N_bk
        <<endl;
*/
    ///////////////////////
    double sumbkg = fBkg->Integral     (ZWind_low, ZWind_high, allFuncPars);
    double sumerr = fBkg->IntegralError(ZWind_low, ZWind_high, allFuncPars, allFuncErrs);
    cout<<"sumbkg = "<<sumbkg<<" +/- "<<sumerr<<endl;

    /*
    double sumbkg = 0;
    double sumerr_sq = 0;
    for(int bin=0;bin<nozpeakhist->GetXaxis()->GetNbins();bin++) {
        float x = nozpeakhist->GetBinCenter(bin);
        if(x>=ZWind_low && x<=ZWind_high) {
            sumbkg += max(0,nozpeakhist->GetBinContent(bin));
            sumerr_sq += TMath::Power( nozpeakhist->GetBinError(bin),2 );
            cout<<x<<"\t"<<nozpeakhist->GetBinContent(bin)<<"\t"<<nozpeakhist->GetBinError(bin)<<endl;
        }
        
    }
    */

    /*
    //binned
    double bkgerrs[2];
    for (int i = 0;i<2;i++) bkgerrs[i] = allfunc->GetParError(i);
    
    double sumbkg = 0;
    double sumerr_sq = 0;
    TF1* errlinearfunc = new TF1("errlinearfunc",error_linear,FitWind_low,FitWind_high,2);
    errquadfunc->SetParameters(bkgerrs);
    for(int bin=0;bin<allhist->GetXaxis()->GetNbins();bin++) {
        float x = allhist->GetBinCenter(bin);
        if(x>=78.67 && x<=103.68) {
            sumbkg += fBkg->Eval(x);
            sumerr_sq += TMath::Power(errquadfunc->Eval(x),2);
            cout<<x<<"\t"<<fBkg->Eval(x)<<"\t"<<errquadfunc->Eval(x)<<endl;
        }
        
    }
    */

} 
