#include "sPhenixStyle.h"
void drawSTAR(double x,double y,double size,int font, int color)
{
  TLatex lat;
  // lat.SetTextSize(0.05);
  lat.SetTextSize(size);
  // lat.SetTextFont(72);
  lat.SetTextFont(font);
  // lat.SetTextColor(kRed);
  lat.SetTextColor(color);
  lat.DrawLatexNDC ( x, y, "STAR Preliminary");
}
void drawIncE()
{
   SetsPhenixStyle();
   double ptbin[9]={0.2,0.4,0.65,0.85,1,1.2,1.6,2,2.8};
   TFile* file = new TFile("out.root");
   TProfile* pf = (TProfile*)file->Get("pEv2");
   pf = (TProfile*)pf->Rebin(8,"pEv2",ptbin);
   pf->SetMarkerColor(kRed);
   pf->SetMarkerStyle(kFullSquare);
   pf->Draw();
   pf->GetYaxis()->SetRangeUser(0.02,0.16);
   drawLatex(0.5,0.35 , "Au+Au #sqrt{s_{NN}} = 54.4 GeV",0.06);
   drawLatex(0.65,0.25 , "0-60\%",0.06);
   TLegend* leg = new TLegend(0.2,0.88,0.5,0.75);
   leg->SetTextSize(0.06);
   leg->AddEntry( pf, "e^{Inc} v_{2}" , "pe");
   leg->Draw();
   drawSTAR(0.6,0.85,0.06,42,kGray+1);
   gPad->SaveAs("fig/IncEv2.pdf");
   gPad->SaveAs("fig/IncEv2.png");
}
