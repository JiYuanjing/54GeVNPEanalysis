#include "../sPhenixStyle.h"
void drawTheorySideLine(TGraphErrors *gr, int color)
{
    int np = gr->GetN();
    double *x = new double[np];
    double *y = new double[np];
    double *yL = new double[np];
    double *yH = new double[np];
    for(int i=0; i<np; i++) {
        double eyH, eyL;
        gr->GetPoint(i, x[i], y[i]);
        eyH = gr->GetErrorY(i);
        eyL = gr->GetErrorY(i);
        yL[i] = y[i] - eyL;
        yH[i] = y[i] + eyH;
    }
    TGraph *gL = new TGraph(np, x, yL);
    gL->SetLineColor(color);
    gL->SetLineWidth(2);
    gL->Draw("sameL");
    TGraph *gH = new TGraph(np, x, yH);
    gH->SetLineWidth(2);
    gH->SetLineColor(color);
    gH->Draw("sameL");
}
void compBD2e()
{
  SetsPhenixStyle();
  TFile* f = new TFile("fonllD2e.root");
  TGraphErrors* hD = (TGraphErrors*)f->Get(Form("gErrFonllD2e"));
  f->Close();
  f = new TFile("fonllB2e.root");
  TGraphErrors* hB = (TGraphErrors*)f->Get(Form("gErrFonllB2e"));
  f->Close();

  /* hD->SetMarkerSize(1.2); */
  /* hD->SetMarkerSize(0); */
  /* hD->SetMarkerStyle(kFullCircle); */
  /* hD->SetMarkerColor(kRed); */
  hD->SetFillColor(kRed-7);
  hD->SetLineColor(kRed-7);
  hD->SetFillStyle(3354);
  /* hD->SetLineWidth(2); */
  hD->GetYaxis()->SetTitleOffset(1.15);
  hD->GetYaxis()->SetRangeUser(1,1e8);
  hD->Draw("E3A");
  /* hD->GetXaxis()->SetRangeUser( 0, 3); */
  /* hB->SetMarkerSize(1.2); */
  /* hB->SetMarkerStyle(kFullCircle); */
  hB->SetFillColor(kTeal-1);
  hB->SetLineColor(kTeal-1);
  hB->SetFillStyle(3345);
  hB->Draw("E3same");
  drawTheorySideLine(hB,kTeal-1);
  drawTheorySideLine(hD,kRed-7);
  hD->Draw("E3same");
  hB->GetXaxis()->SetRangeUser( 0, 3);
  gPad->SetLogy(1);

  hD->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hD->GetYaxis()->SetTitle("d#sigma/dp_{T} (pb/GeV)");

  TLegend* leg = new TLegend(0.68,0.58,0.83,0.75);
  leg->SetTextSize(0.06);
  leg->AddEntry( hD, "D->e","f" );
  leg->AddEntry( hB, "B->e","f" );
  leg->Draw();
  drawLatex(0.65,0.77,"p+p 62.4 GeV",0.06);
  drawLatex(0.6,0.85,"FONLL calculation",0.06);
}
void compBD2eHist()
{
  SetsPhenixStyle();
  TFile* f = new TFile("fonllD2e.root");
  TH1F* hD = (TH1F*)f->Get(Form("hD2e"));
  hD->SetDirectory(0);
  f->Close();
  f = new TFile("fonllB2e.root");
  TH1F* hB = (TH1F*)f->Get(Form("hB2e"));
  hB->SetDirectory(0);
  f->Close();

  hD->SetMarkerSize(1.2);
  hD->SetMarkerStyle(kFullCircle);
  hD->SetMarkerColor(kRed);
  hD->SetFillColorAlpha(0.3,kRed);
  hD->SetFillStyle(1001);
  hD->GetYaxis()->SetTitleOffset(1.15);
  hD->GetYaxis()->SetRangeUser(1,1e8);
  hD->Draw("E3");
  hD->GetXaxis()->SetRangeUser( 0, 3);
  hB->SetMarkerSize(1.2);
  hB->SetMarkerStyle(kFullCircle);
  hB->SetFillColorAlpha(0.3,kBlack);
  hB->SetFillStyle(1001);
  hB->Draw("E3same");
  hB->GetXaxis()->SetRangeUser( 0, 3);
  gPad->SetLogy(1);

  hD->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hD->GetYaxis()->SetTitle("d#sigma/dp_{T} (pb/GeV)");

  TLegend* leg = new TLegend(0.68,0.58,0.83,0.75);
  leg->SetTextSize(0.06);
  leg->AddEntry( hD, "D->e","p" );
  leg->AddEntry( hB, "B->e","p" );
  leg->Draw();
  drawLatex(0.65,0.77,"p+p 62.4 GeV",0.06);
  drawLatex(0.6,0.85,"FONLL calculation",0.06);
}
