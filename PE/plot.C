#include "myStyle.h"
void plot()
{
  // SetMyStyle();
  // color palette - manually define 'kBird' palette only available in ROOT 6
	const Int_t NRGBs = 5;
	const Int_t NCont = 255;

	Int_t fcol;

	Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
	Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
	fcol = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	gStyle->SetNumberContours(NCont);
  TFile* file = new TFile("PE_0.root");
  TH2F* h2 = (TH2F*)file->Get("hdEdx_tofpt");
  h2->Draw("colz");
  gPad->SetLogz();
}
