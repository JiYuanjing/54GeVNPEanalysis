#include "myStyle.h"
void plotDecayLength()
{
  SetMyStyle();
  TFile* f100 = TFile::Open("output/fdecatLforPlots_1.root");
  TFile* f85 =  TFile::Open("output/fdecatLforPlots_0.7.root");
  TFile* f115 =  TFile::Open("output/fdecatLforPlots_1.3.root");

  TH1F* hdata = (TH1F*)f100->Get("hdata")->Clone("hdata");
  hdata->SetDirectory(0);
  TH1F* hMC100 = (TH1F*)f100->Get("hrc")->Clone("hrc");
  hMC100->SetDirectory(0);
  TH1F* hMC85 = (TH1F*)f85->Get("hrc")->Clone("hrclow");
  hMC85->SetDirectory(0);
  hMC85->SetMarkerColor(kRed);
  TH1F* hMC115 = (TH1F*)f115->Get("hrc")->Clone("hrchigh");
  hMC115->SetDirectory(0);
  hMC115->SetMarkerColor(kBlue);
  f100->Close();
  f85->Close();
  f115->Close();
  // int scalebin = 3;
  // hMC100->Rebin(scalebin);
  // hMC115->Rebin(scalebin);
  // hMC85->Rebin(scalebin);
  // hMC100->Scale(1.0/scalebin);
  // hMC115->Scale(1.0/scalebin);
  // hMC85->Scale(1.0/scalebin);

  hMC100->Scale(0.9/hdata->GetMaximum());
  hMC115->Scale(0.9/hdata->GetMaximum());
  hMC85->Scale(0.9/hdata->GetMaximum());
  hdata->Scale(0.9/hdata->GetMaximum());

  hMC100 = reBinHist(2,0,70,hMC100);
  hMC115 = reBinHist(2,0,70,hMC115);
  hMC85 = reBinHist(2,0,70,hMC85);
  hdata = reBinHist(8,0,70,hdata);

  int const tmp  = hMC100->GetNbinsX();
  int const nbins = tmp;
  double x[nbins],y[nbins],errl[nbins],errh[nbins];
  for (int i=1;i<=hMC100->GetNbinsX();i++) 
  {
      x[i-1] = hMC100->GetBinCenter(i);
      y[i-1] = hMC100->GetBinContent(i);
      double e1 = hMC115->GetBinContent(i);
      double e2 = hMC85->GetBinContent(i);
      double estat = hMC100->GetBinError(i);;
      if (e1>e2) {  
        errh[i-1] = e1-y[i-1];
        errl[i-1] = y[i-1]-e2;
        // errh[i-1] = e1;
        // errl[i-1] = e2;
      }
      else {
        errl[i-1] = y[i-1]-e1;
        errh[i-1] = e2-y[i-1];
        // errl[i-1] = e1;
        // errh[i-1] = e2;
      }
      errh[i-1] = sqrt(errh[i-1]*errh[i-1]+estat*estat);
      errl[i-1] = sqrt(errl[i-1]*errl[i-1]+estat*estat);
  }
  TGraphAsymmErrors* gasymm = new TGraphAsymmErrors(nbins,x,y,0,0,errl,errh);
  gasymm->SetName("gDecayL");
  hdata->SetMarkerSize(1);
  hdata->GetYaxis()->SetTitle("Arb. unit");
  hdata->GetXaxis()->SetTitle("Pair decay length [cm]");
  hdata->Draw("p");
  gasymm->SetFillStyle(3004);
  gasymm->SetFillColor(kBlue);
  gasymm->Draw("same3");
  hMC85->SetMarkerColor(kBlue);
  // hMC85->Draw("same");

  hMC100->SetMarkerSize(1.2);
  hMC100->SetMarkerColor(kRed);
  hMC100->SetLineColor(kRed+2);

  hMC100->SetFillColor(kRed);
  hMC100->SetFillStyle(3009);
  // hMC100->DrawCopy("sameE3");
  // hMC100->DrawCopy("sameC");
  // hMC100->DrawCopy("samep");
  hMC115->SetMarkerColor(kGreen);
  // hMC115->DrawCopy("samep");
  
  TLegend* leg = new TLegend(0.65,0.67,0.87,0.85);
  // leg->AddEntry(hMC100,"MC","lef");
  leg->AddEntry(gasymm,"MC","lef");
  leg->AddEntry(hdata,"Data","ep");
  leg->Draw();
  drawLatex(0.6,0.6,"0.2 GeV<p_{T}<2.5 GeV");
  drawSTAR(0.6,0.5);
  
  gPad->SaveAs("pdf/decayL.pdf");
  gPad->SaveAs("pdf/decayL.png");
  TFile* fout = new TFile("output/gDecayL.root","recreate");
  hMC100->Write();
  hMC115->Write();
  hMC85->Write();
  hdata->Write();
  gasymm->Write();
  fout->Close();

}
void drawSTAR(double x,double y)
{
  TLatex lat;
  lat.SetTextSize(0.05);
  lat.SetTextFont(72);
  lat.SetTextColor(kRed);
  lat.DrawLatexNDC ( x, y, "STAR Preliminary");

}
TH1F* reBinHist(double xmid,double xlow, double xhigh, TH1F* h)
{
  int const nMax = 1000;
  double xedge[nMax];
  int nbin=0;
  xedge[0]=xlow;
  while (xedge[nbin]<=xmid)
  {
    nbin++;
    xedge[nbin]=xedge[nbin-1]+h->GetBinWidth(1);
  }
  int turningpoint = nbin+1;
  double rebin = 3;
  while (xedge[nbin]<xhigh)
  {
    nbin++;
    xedge[nbin]=xedge[nbin-1]+h->GetBinWidth(1)*rebin;
    // h->SetBinContent(h->FindBin(xedge[nbin]-h->GetBinWidth(1)*1.5),h->GetBinContent(h->FindBin(xedge[nbin]-h->GetBinWidth(1)*1.5))*0.5);
    // h->SetBinContent(h->FindBin(xedge[nbin]-h->GetBinWidth(1)*0.5),h->GetBinContent(h->FindBin(xedge[nbin]-h->GetBinWidth(1)*0.5))*0.5);
  }
  // return (TH1F*)h->Rebin(nbin,h->GetName(),xedge);
  h = (TH1F*)h->Rebin(nbin,h->GetName(),xedge);
  for (int ibin=turningpoint;ibin<h->GetNbinsX();ibin++)
  {
    h->SetBinContent(ibin,h->GetBinContent(ibin)/(1.0*rebin));
  }
  return h;
}
