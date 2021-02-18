#include "myStyle.h"
void plotPartPt()
{
  SetMyStyle();
  TFile* f100 = TFile::Open("output/fPartnerPt_1.root");
  TFile* f85 =  TFile::Open("output/fPartnerPt_0.7.root");
  TFile* f115 =  TFile::Open("output/fPartnerPt_1.3.root");

  TH1F* hdata = (TH1F*)f100->Get("hdata_x")->Clone("hdata");
  hdata->SetDirectory(0);
  TH1F* hMC100 = (TH1F*)f100->Get("hrc_x")->Clone("hrc");
  hMC100->SetDirectory(0);
  TH1F* hMC85 = (TH1F*)f85->Get("hrc_x")->Clone("hrclow");
  hMC85->SetDirectory(0);
  hMC85->SetMarkerColor(kRed);
  TH1F* hMC115 = (TH1F*)f115->Get("hrc_x")->Clone("hrchigh");
  hMC115->SetDirectory(0);
  hMC115->SetMarkerColor(kBlue);
  f100->Close();
  f85->Close();
  f115->Close();
 
  hMC100->Scale(0.9/hdata->GetMaximum());
  hMC115->Scale(0.9/hdata->GetMaximum());
  hMC85->Scale(0.9/hdata->GetMaximum());
  hdata->Scale(0.9/hdata->GetMaximum());

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
      // cout<< estat<<endl;
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
  int startpoint = hMC100->GetXaxis()->FindBin(0.25);
  TGraphAsymmErrors* gasymm = new TGraphAsymmErrors(nbins-startpoint+1,x+startpoint-1,y+startpoint-1,0,0,errl+startpoint-1,errh+startpoint-1);
  gasymm->SetName("gPartPt");
  hdata->SetMarkerSize(1);
  hdata->GetYaxis()->SetTitle("Arb. unit");
  hdata->GetXaxis()->SetTitle("Partner e p_{T} [GeV/c]");
  hdata->Draw("p");
  gasymm->SetFillStyle(3004);
  gasymm->SetFillColor(kBlue);
  gasymm->Draw("same3L");
  hMC85->SetMarkerColor(kBlue);
  // hMC85->Draw("same");

  hMC100->SetMarkerSize(1.2);
  hMC100->SetMarkerColor(kRed);
  hMC100->SetLineColor(kRed+2);

  // hMC100->SetFillColor(kRed);
  // hMC100->SetFillStyle(3009);
  // hMC100->DrawCopy("sameE3");
  // hMC100->DrawCopy("sameC");
  // hMC100->DrawCopy("samepH");
  hMC115->SetMarkerColor(kGreen);
  // hMC115->DrawCopy("samep");
  
  TLegend* leg = new TLegend(0.65,0.67,0.87,0.85);
  // leg->AddEntry(hMC100,"MC","lef");
  leg->AddEntry(gasymm,"MC","lef");
  leg->AddEntry(hdata,"Data","ep");
  leg->Draw();
  drawLatex(0.6,0.6,"0.2 GeV<p_{T}<2.5 GeV");
  drawSTAR(0.6,0.5);
  
  gPad->SaveAs("pdf/PartePt.pdf");
  gPad->SaveAs("pdf/PartePt.png");
  TFile* fout = new TFile("output/gPartPt.root","recreate");
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
