#include "sPhenixStyle.h"
void drawDe()
{
  SetsPhenixStyle();
  gStyle->SetPalette(1);
  TFile * f = new TFile("noresDecay25_fonll62.root");
  TCanvas* c = new TCanvas("c","c");
  TH2F* h = (TH2F*)f->Get("D0vsE");
  h->GetXaxis()->SetRangeUser(0,7.5);
  h->GetYaxis()->SetRangeUser(0,4);
  h->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
  h->GetYaxis()->SetTitle("e p_{T} (GeV/c)");
  h->GetXaxis()->SetTitleOffset(1.15);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetZaxis()->SetRangeUser(1e10,5e18);

  h->Draw("colz");

  gPad->Update();
  TPaletteAxis* palette = (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(1.2);
  palette->SetX2NDC(1.5);

  gPad->SetLogz(1);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  /* drawLatex(0.7,0.7,"D->e"); */
}
