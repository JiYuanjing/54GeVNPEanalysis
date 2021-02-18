#include "myStyle.h"
void plotCorr()
{
  SetMyStyle();
  TFile* file = TFile::Open("test_200GeVNCR.root"); 
  TProfile* pf = (TProfile*)file->Get("D_e"); 
  pf->GetXaxis()->SetTitle("e^{D} p_{T} [GeV/c]");
  double xmax = 3.5,xlow=0;
  pf->GetXaxis()->SetRangeUser(xlow,xmax);
  /* pf->GetXaxis()->SetTitle(""); */
  pf->GetYaxis()->SetTitle("<cos(#phi_{D#rightarrow e}#minus #phi_{D})>");
  pf->Draw("c");
  TH1F* h = (TH1F*)pf->Clone();
  h->Draw("c");
  drawLine(xlow,0.95, xmax,0.95, 2,9,kRed);
  drawLatex(0.6,0.35,"PYTHIA",0.08);
  drawLatex(0.58,0.25,"p+p 54 GeV",0.08);
  drawLatex(0.22,0.8,"0.95",0.06,kRed);
  gPad->SaveAs("HFeCorr.pdf");
  gPad->SaveAs("HFeCorr.png");
}
