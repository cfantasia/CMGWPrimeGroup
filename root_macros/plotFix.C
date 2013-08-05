{

  TGraphErrors* gBil = new TGraphErrors(10); gBil->SetLineColor(kGreen);
  TGraphErrors* gHal = new TGraphErrors(10); gHal->SetLineColor(kBlue);
  TGraphErrors* gFil = new TGraphErrors(10); gFil->SetLineColor(kRed);
  TGraphErrors* gOld = new TGraphErrors(10);
  TMultiGraph* mg = new TMultiGraph("mg", ";p_{T}^{Gen} (GeV);Resolution=RMS(p_{T}^{Reco} - p_{T}^{Gen}) (GeV)");
  //TFile *_file0 = TFile::Open("WprimeWZ-Resolution.root");
  //TFile *_file0 = TFile::Open("WprimeWZ-Resolution-Skimmed.root");
  TFile *_file0 = TFile::Open("WprimeWZ-Resolution-test.root");
  _file0->ls();

  for(int i=1; i<=10; i++){
    float pt = i*100. - 50.;
    float res, reserr;
    res    = hbilResW2D->ProjectionX("", i,i)->GetRMS();
    reserr = hbilResW2D->ProjectionX("", i,i)->GetRMSError();
    gBil->SetPoint(i, pt, res);
    gBil->SetPointError(i, 50, reserr);

    res    = hhalResW2D->ProjectionX("", i,i)->GetRMS();
    reserr = hhalResW2D->ProjectionX("", i,i)->GetRMSError();
    gHal->SetPoint(i, pt, res);
    gHal->SetPointError(i, 50, reserr);

    res    = hfilResW2D->ProjectionX("", i,i)->GetRMS();
    reserr = hfilResW2D->ProjectionX("", i,i)->GetRMSError();
    gFil->SetPoint(i, pt, res);
    gFil->SetPointError(i, 50, reserr);

    res    = holdResW2D->ProjectionX("", i,i)->GetRMS();
    reserr = holdResW2D->ProjectionX("", i,i)->GetRMSError();
    gOld->SetPoint(i, pt, res);
    gOld->SetPointError(i, 50, reserr);
  }
  
  mg->Add(gBil);
  //mg->Add(gFil);
  mg->Add(gHal);
  mg->Add(gOld);
  mg->Draw("alp");

  TLegend *leg = new TLegend(0.2, 0.43,0.3, 0.89,"");
  leg->SetTextSize(0.05);
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(gBil, "Bil", "le");
  leg->AddEntry(gHal, "Hal", "le");
  leg->AddEntry(gOld, "Old", "le");
  leg->Draw();

}
