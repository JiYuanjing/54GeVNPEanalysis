#include "myStyle.h"
void plotDecayLength()
{
  SetMyStyle();
  // TFile* f110 = TFile::Open("output/fdecayL_smear1.0.root");
  TFile* f110 = TFile::Open("output/fdecayLforPlots_1.root");
  TFile* f100 = TFile::Open("output/fdecayL_smear0.5.root");
  TFile* f85 =  TFile::Open("output/fdecayL_def.root");

  TH1F* hdata = (TH1F*)f100->Get("hdata");
  hdata->SetDirectory(0);
  TH1F* hMC100 = (TH1F*)f100->Get("hrc")->Clone("hrc100");
  hMC100->SetDirectory(0);
  TH1F* hMC110 = (TH1F*)f100->Get("hrc")->Clone("hrc100");
  hMC110->SetDirectory(0);
  TH1F* hMC85 = (TH1F*)f85->Get("hrc")->Clone("hrc85");
  hMC85->SetDirectory(0);
  hMC85->SetMarkerColor(kRed);
  f110->Close();
  f100->Close();
  f85->Close();

  // TGraphAsymmErrors* gasymm = new TGraphAsymmErrors(nbins,x,y,0,0,errl,errh);
  // for (int i=1;i<=hMC100->GetNbinsX();i++) 
  // {
  //     // x[i-1] = hMC100->GetBinCenter(i);
  //     // y[i-1] = hMC100->GetBinContent(i);
  //     double y = hMC100->GetBinContent(i);
  //     // double e1 = hMC115->GetBinContent(i);
  //     // double e2 = hMC85->GetBinContent(i);
  //     // double err1 = fabs(e1-y)>fabs(e2-y)?fabs(e1-y):fabs(e2-y);
  //     // double err2 = hMC100->GetBinError(i);;
  //     hMC100->SetBinError(i,0);
  // }
   
  // hdata->SetMarkerStyle(24);
  hdata->SetMarkerSize(1);
  hdata->GetYaxis()->SetTitle("Arb. unit");
  hdata->GetXaxis()->SetTitle("Pair decay length [cm]");
  hMC85->Draw();
  hdata->Draw("psame");
  // gasymm->SetFillStyle(3004);
  // gasymm->SetFillColor(kRed);
  // gasymm->Draw("AsameF");
  hMC85->Draw("same");

  hMC100->SetMarkerColor(kBlue);
  hMC100->SetLineColor(kBlue);
  hMC100->SetMarkerStyle(20);
  hMC100->SetMarkerSize(1.2);
  // hMC100->SetFillStyle(3009);
  // hMC100->DrawCopy("sameE3");
  hMC100->DrawCopy("samep");
  
  hMC110->SetLineColor(kGreen+2);
  hMC110->SetMarkerColor(kGreen+2);
  hMC110->SetLineColor(kGreen+2);
  hMC110->SetMarkerStyle(20);
  hMC110->SetMarkerSize(1.2);
  hMC110->DrawCopy("samep"); 

  TLegend* leg = new TLegend(0.65,0.67,0.87,0.85);
  leg->AddEntry(hMC100,"MC 0.5% smear","lefp");
  leg->AddEntry(hMC85,"MC default","lefp");
  leg->AddEntry(hdata,"Data","ep");
  leg->Draw();
  drawLatex(0.6,0.6,"0.2 GeV<p_{T}<2.5 p_{T}");
  drawSTAR(0.6,0.5);
  // TLegend* leg = new TLegend(0.7,0.7,0.88,0.88);
  // leg->AddEntry();
  // leg->AddEntry();
  // leg->AddEntry();
  // leg->Draw();
  
  gPad->SaveAs("pdf/decayL.pdf");
  gPad->SaveAs("pdf/decayL.png");
}
void drawSTAR(double x,double y)
{
  TLatex lat;
  lat.SetTextSize(0.05);
  lat.SetTextFont(72);
  lat.SetTextColor(kRed);
  lat.DrawLatexNDC ( x, y, "STAR Preliminary");

}
