//Usage: root -b -l -q 'CalcFakeRate.C+(file, useData)'
//e.g. root -b -l -q 'CalcFakeRate.C+("../../../WZFakeRateMuon.root", 1)'

#include <vector>
#include "THStack.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "../Limits/consts.h"

void GetNDC(TCanvas & c, double xin, double yin, double & x, double & y) {
  c.Update();//this is necessary!
  //printf("x %g, x1 %g, x2 %g, xr %g\n", xin, c.GetX1(), c.GetX2(), c.GetX2()-c.GetX1());
  //printf("y %g, y1 %g, y2 %g, yr %g\n", yin, c.GetY1(), c.GetY2(), c.GetY2()-c.GetY1());
  x = (c.XtoPad(xin) - c.GetX1())/(c.GetX2()-c.GetX1());
  y = (c.YtoPad(yin) - c.GetY1())/(c.GetY2()-c.GetY1());

}

void
CalcFakeRate(string infile, bool useData=true){  
  TFile *f = TFile::Open(infile.c_str(), "read"); assert(f);
  TH2F* hFakeRateNum = NULL;
  TH2F* hFakeRateAll = NULL;
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetTextFont(132);
  gStyle->SetTextSize(1.2);
  //gROOT->ForceStyle();

  bool useElectrons = infile.find("Muon") == string::npos;

  vector<string> allsamples;
  if(!useData){
    allsamples.push_back("WJetsToLNu");//corrupted file
    allsamples.push_back("TTJets");
    allsamples.push_back("ZZ");
    allsamples.push_back("GVJets");
    allsamples.push_back("WWTo2L2Nu");
    allsamples.push_back("WZJetsTo3LNu");
    allsamples.push_back("DYJetsToLL");
  }else{
    allsamples.push_back("data");
  }
  
  float offset = 0.01;
  float eta[] = {0., 1.5, 2.5}; int neta = 3;
  if(!useElectrons){
    eta[1] = 2.4;
    neta = 2;//eta[2] = 2.4;
  }

  float pt [] = {0., 20., 25., 30., 50., 1000}; const int npt  = 6;
  hFakeRateNum = new TH2F("hFakeRateNum", "hFakeRateNum;p_{T} (GeV);#eta", npt-1, pt, neta-1, eta);
  hFakeRateAll = new TH2F("hFakeRateAll", "hFakeRateAll;p_{T} (GeV);#eta", npt-1, pt, neta-1, eta);
  string title = useElectrons ? "Electron Fake Rate" : "Muon Fake Rate";
  title += useData ? " from Data" : " from MC";
  hFakeRateNum->SetTitle(title.c_str());

/////////This is for eta, pt agnotistic Fake Rate Calc
  float in_tot(0), out_tot(0);
  for(unsigned i=0; i<allsamples.size(); ++i){
    string hist_name = allsamples[i] + "/hNumEvts";
    TH1F* hist = (TH1F*) f->Get(hist_name.c_str()); assert(hist);
    int lastbin = hist->GetNbinsX();
    float in(0), out(0);

    in  = hist->GetBinContent(lastbin);
    out = hist->GetBinContent(lastbin-1);

    in_tot += in;
    out_tot += out;

    //cout<<"Sample: "<<allsamples[i]<<" = "<<in/out<<" = pass/total = "<<in<<"/"<<out<<endl;
    printf("  Sample: %s = %.2f : pass/total = %.2f / %.2f \n",allsamples[i].c_str(), in/out, in, out);
  }
  float  eff = in_tot/out_tot;
  float deff = TMath::Sqrt(eff * (1-eff)/out_tot);
  //cout<<"Total: "<<eff*100<<"% +/- "<<deff*100<<"% = in_tot/out_tot*100% = "<<in_tot<<"/"<<out_tot<<"*100%\n";
  printf("Total: %.2f%% +/- %.2f%% = in_tot/out_tot*100%% = %.2f/%.2f\n", eff*100, deff*100, in_tot, out_tot);


/////////This is for 2D Fake Rate Calc
  for(unsigned i=0; i<allsamples.size(); ++i){
    string hist_name = allsamples[i] + "/hNumEvts";
    TH1F* hist = (TH1F*) f->Get(hist_name.c_str()); assert(hist);
    int lastbin = hist->GetNbinsX();
    float in(0), out(0);
    string binNameNum   = hist->GetXaxis()->GetBinLabel(lastbin);
    string binNameDenom = hist->GetXaxis()->GetBinLabel(lastbin-1);
      
    string hist_nameNum = allsamples[i] + "/hEtaVsPt_" + binNameNum;
    TH2F* histNum = (TH2F*) f->Get(hist_nameNum.c_str()); assert(histNum);
      
    string hist_nameDenom = allsamples[i] + "/hEtaVsPt_" + binNameDenom;
    TH2F* histDenom = (TH2F*) f->Get(hist_nameDenom.c_str()); assert(histDenom);

    //cout<<" Total from 2D Plot is "<<histNum->Integral()<<" / "<<histDenom->Integral()<<endl; //By default, integral doesn't count over/underflow
    //cout<<" Total from 2D Plot is "<<histNum->GetEntries()<<" / "<<histDenom->GetEntries()<<endl;

    for(int ieta=0; ieta<neta-1; ++ieta){
      float ymin = eta[ieta]; 
      float ymax = eta[ieta+1]-offset;
      for(int ipt=0; ipt<npt-1; ++ipt){
        float xmin = pt[ipt]; 
        float xmax = pt[ipt+1];
          
        int xminbin,xmaxbin,yminbin,ymaxbin;
        xminbin = histNum->GetXaxis()->FindBin(xmin) + (ipt!=0);//Avoid overlap except for first bin
        xmaxbin = histNum->GetXaxis()->FindBin(xmax);
        yminbin = histNum->GetYaxis()->FindBin(ymin);
        ymaxbin = histNum->GetYaxis()->FindBin(ymax);
        
/*
        cout<<"("<<pt[ipt]<<", "<<eta[ieta]<<")\t"
            <<xmin<<"-"<<xmax<<":"
            <<ymin<<"-"<<ymax<<":"
            <<"\t";
        cout<<"("<<ipt<<", "<<ieta<<")\t"
            <<xminbin<<"-"<<xmaxbin<<":"
            <<yminbin<<"-"<<ymaxbin
            <<endl;
*/        
        in  = histNum  ->Integral(xminbin, xmaxbin, yminbin, ymaxbin);
        out = histDenom->Integral(xminbin, xmaxbin, yminbin, ymaxbin);
          
        //Cory: Deal with negative eta (Y axis)
        yminbin = histNum->GetYaxis()->FindBin(-1*ymax);
        ymaxbin = histNum->GetYaxis()->FindBin(-1*ymin)-1;//avoid overlap
        /*
        cout<<"("<<pt[ipt]<<", "<<eta[ieta]<<")\t"
            <<xmin<<"-"<<xmax<<":"
            <<ymin<<"-"<<ymax<<":"
            <<"\t";
        cout<<"("<<ipt<<", "<<ieta<<")\t"
            <<xminbin<<"-"<<xmaxbin<<":"
            <<yminbin<<"-"<<ymaxbin
            <<endl;
        */
        in  += histNum  ->Integral(xminbin, xmaxbin, yminbin, ymaxbin);
        out += histDenom->Integral(xminbin, xmaxbin, yminbin, ymaxbin);

        //cout<<"("<<pt[ipt]<<","<<eta[ieta]<<") "<<in<<"/"<<out<<endl;
        //cout<<"("<<(pt[ipt]+pt[ipt+1])/2<<","<<(eta[ieta]+eta[ieta+1])/2<<") "<<in<<"/"<<out<<endl;
        hFakeRateNum   ->Fill( (pt[ipt]+pt[ipt+1])/2, (eta[ieta]+eta[ieta+1])/2, in);
        hFakeRateAll->Fill( (pt[ipt]+pt[ipt+1])/2, (eta[ieta]+eta[ieta+1])/2, out);

      }
    }

  }

  if(hFakeRateNum){
    //TGraphAsymmErrors* hFakeRateEff = new TGraphAsymmErrors(hFakeRateNum->GetArray(), hFakeRateAll->GetArray(), "");
    //TGraphAsymmErrors* hFakeRateEff = new TGraphAsymmErrors(hFakeRateNum, hFakeRateAll, "");
    TH2F* hFakeRate = (TH2F*) hFakeRateNum->Clone("hFakeRate");
    //hFakeRate->Divide(hFakeRateAll);
    //hFakeRate->Scale(100.); //make is a percent
    TCanvas c1;
    c1.SetLogx();


    //gStyle->SetTextFont(132);
    //gStyle->SetTextSize(1.2);
    hFakeRate->SetMarkerSize(3);
    //gStyle->SetPaintTextFormat("3.0f m");
    gStyle->SetPaintTextFormat("4.0f");
    hFakeRate->Draw("colz");

    TLatex latexLabel;
    latexLabel.SetNDC();
    latexLabel.SetTextSize(0.04);
    latexLabel.SetTextFont(42);
    latexLabel.SetTextAngle(90);
    latexLabel.SetTextAlign(22);
    for(int ybin = 1; ybin<= hFakeRate->GetYaxis()->GetNbins(); ++ybin){
      for(int xbin = 1; xbin<= hFakeRate->GetXaxis()->GetNbins(); ++xbin){
        float xcen = hFakeRate->GetXaxis()->GetBinCenter(xbin);
        float ycen = hFakeRate->GetYaxis()->GetBinCenter(ybin);
        
        double xpos, ypos;
        GetNDC(c1, xcen, ycen, xpos, ypos);
        int bin = hFakeRate->GetBin(xbin, ybin);
        float pass = hFakeRateNum->GetBinContent(bin);
        float tot  = hFakeRateAll->GetBinContent(bin);
        float mean    = pass / tot;//hFakeRateEff->GetY()       [bin-1];
        float errUp   = TEfficiency::ClopperPearson(tot, pass, 0.68, true);//hFakeRateEff->GetErrorYhigh(bin-1);
        float errDown = TEfficiency::ClopperPearson(tot, pass, 0.68, false);//hFakeRateEff->GetErrorYlow(bin-1);
        hFakeRate->SetBinContent(bin, mean);
        //hFakeRate->SetBinError(bin, err);
        float err = (errUp + errDown)/2;
        printf("pt %3.0f eta %.1f, bin %2i, mean %.4f err %.4f (%.2f / %.2f)\n", 
               hFakeRate->GetXaxis()->GetBinLowEdge(xbin), hFakeRate->GetYaxis()->GetBinLowEdge(ybin), bin, mean, err, pass, tot);
        if(xbin != 1) latexLabel.DrawLatex(xpos, ypos, Form("%.0f^{+%.0f}_{-%.0f}%% (%.0f/%.0f)(%.f-%.f)GeV",mean*100, errUp*100, errDown*100, pass, tot, pt[xbin-1],pt[xbin]));
      }
    }
    
    string outName = "FakeRate";
    outName += useElectrons ? "Elec" : "Muon";
    outName += useData ? "Data" : "MC";
    outName += ".png";
    c1.Print(outName.c_str());

    //Now make profiles
    c1.Clear();
    c1.SetLogx(0);
    c1.SetGrid();

    TH1D* hFakeRateNumEta = hFakeRateNum->ProjectionY("hFakeRateNumEta", 0, -1, "e");
    TH1D* hFakeRateAllEta = hFakeRateAll->ProjectionY("hFakeRateAllEta", 0, -1, "e");
    //TH1D* hFakeRateEta = (TH1D*) hFakeRateNumEta->Clone("hFakeRateEta");

    int n = hFakeRateAllEta->GetXaxis()->GetNbins();
    TGraphAsymmErrors* hFakeRateEta = new TGraphAsymmErrors(n);//hFakeRateNumEta, hFakeRateAllEta);
    hFakeRateEta->SetMaximum(1.);
    hFakeRateEta->SetMinimum(0.);
    hFakeRateEta->SetTitle((title + ";#eta;Rate").c_str());

    //hFakeRateEta->SetMarkerStyle(21);    
    for(int xbin = 1; xbin <= n; ++xbin){
      float x = hFakeRateAllEta->GetXaxis()->GetBinCenter(xbin);
      float pass = hFakeRateNumEta->GetBinContent(xbin);
      float tot  = hFakeRateAllEta->GetBinContent(xbin);
      float mean    = pass / tot;
      float errUp   = TEfficiency::ClopperPearson(tot, pass, 0.68, true);//hFakeRateEff->GetErrorYhigh(bin-1);
      float errDown = TEfficiency::ClopperPearson(tot, pass, 0.68, false);//hFakeRateEff->GetErrorYlow(bin-1);
      hFakeRateEta->SetPoint      (xbin, x, mean);
      hFakeRateEta->SetPointEYhigh(xbin, errUp);
      hFakeRateEta->SetPointEYlow (xbin, errDown);
      //printf("bin: %i, x=%.0f, %.0f/%.0f = %.0f +%.0f -%.0f\n", xbin, x, pass, tot, mean*100, errUp*100, errDown*100);
    }
    
    //hFakeRateEta->Divide(hFakeRateAllEta);
    //hFakeRateEta->Scale(100.);
    hFakeRateEta->Draw("ap*");
    outName = "FakeRateEta";
    outName += useElectrons ? "Elec" : "Muon";
    outName += useData ? "Data" : "MC";
    outName += ".png";
    c1.Print(outName.c_str());

    //Pt projections (split these by eta?)
    c1.Clear();
    c1.SetLogx(1);
    c1.SetGrid();

    TH1D* hFakeRateNumPt = hFakeRateNum->ProjectionX("hFakeRateNumPt", 0, -1, "e");
    TH1D* hFakeRateAllPt = hFakeRateAll->ProjectionX("hFakeRateAllPt", 0, -1, "e");
    //TH1D* hFakeRatePt = (TH1D*) hFakeRateNumPt->Clone("hFakeRatePt");

    n = hFakeRateAllPt->GetXaxis()->GetNbins();
    TGraphAsymmErrors* hFakeRatePt = new TGraphAsymmErrors(n);//hFakeRateNumPt, hFakeRateAllPt);
    hFakeRatePt->SetMaximum(1.);
    hFakeRatePt->SetMinimum(0.);
    hFakeRatePt->SetTitle((title + ";p_{T} (GeV);Rate").c_str());

    for(int xbin = 1; xbin <= n; ++xbin){
      float x = hFakeRateAllPt->GetXaxis()->GetBinCenter(xbin);
      float pass = hFakeRateNumPt->GetBinContent(xbin);
      float tot  = hFakeRateAllPt->GetBinContent(xbin);
      float mean    = pass / tot;
      float errUp   = TEfficiency::ClopperPearson(tot, pass, 0.68, true);//hFakeRateEff->GetErrorYhigh(bin-1);
      float errDown = TEfficiency::ClopperPearson(tot, pass, 0.68, false);//hFakeRateEff->GetErrorYlow(bin-1);
      hFakeRatePt->SetPoint      (xbin, x, mean);
      hFakeRatePt->SetPointEYhigh(xbin, errUp);
      hFakeRatePt->SetPointEYlow (xbin, errDown);
      //printf("bin: %i, x=%.0f, %.0f/%.0f = %.0f +%.0f -%.0f\n", xbin, x, pass, tot, mean*100, errUp*100, errDown*100);
    }
    


    //hFakeRatePt->Divide(hFakeRateAllPt);
    //hFakeRatePt->Scale(100.);
    hFakeRatePt->Draw("ap*");
    outName = "FakeRatePt";
    outName += useElectrons ? "Elec" : "Muon";
    outName += useData ? "Data" : "MC";
    outName += ".png";
    c1.Print(outName.c_str());

    
  }


}
