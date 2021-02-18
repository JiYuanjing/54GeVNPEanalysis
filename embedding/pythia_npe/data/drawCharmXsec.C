#include "../sPhenixStyle.h"
void drawTheorySideLineAs(TGraphAsymmErrors* gr, int color)
{
    int np = gr->GetN();
    double *x = new double[np];
    double *y = new double[np];
    double *yL = new double[np];
    double *yH = new double[np];
    for(int i=0; i<np; i++) {
        double eyH, eyL;
        gr->GetPoint(i, x[i], y[i]);
        eyH = gr->GetErrorYhigh(i);
        eyL = gr->GetErrorYlow(i);
        yL[i] = y[i] - eyL;
        yH[i] = y[i] + eyH;
    }
    TGraph *gL = new TGraph(np, x, yL);
    gL->SetLineColor(color);
    gL->SetLineWidth(2);
    gL->Draw("sameL");
    TGraph *gH = new TGraph(np, x, yH);
    gH->SetLineColor(color);
    gH->SetLineWidth(2);
    gH->Draw("sameL");
}
void drawCharmXsec()
{
 SetsPhenixStyle(); 
  ifstream read;
  read.open("totalXsec.txt");
  char line[100];
  getline(read,line,'\n');
  getline(read,line,'\n');
  cout<<line<<endl;
  double x[10],y[10],errup[10],errlow[10], errsclow[10],errscup[10],errmasslow[10],errmassup[10];
  for (int i=0;i<8;i++)
  {
     read>>x[i]>>y[i]>>errlow[i]>>errup[i]>>errsclow[i]>>errscup[i]>>errmasslow[i]>>errmassup[i]; 
     cout <<x[i]<<" " <<y[i] <<endl;
     errlow[i] = (y[i]-errlow[i])*1e-6;
     errup[i] = (errup[i]-y[i])*1e-6;
     y[i]=y[i]*1e-6;
  }

  TH1F* h = new TH1F("h","h",1,50,15000);
  h->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
  h->GetYaxis()->SetTitle("Total Charm Cross Secton (#mu b)");
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetXaxis()->SetTitleOffset(1.15);
  /* h->GetYaxis()->SetRangeUser(1e6,5e9); */
  h->GetYaxis()->SetRangeUser(1.1,5e3);
  h->Draw("c");

  TGraphAsymmErrors* g = new TGraphAsymmErrors( 8, x, y, 0,0,errlow,errup);
  g->SetMarkerColor(kViolet+7);
  g->SetLineColor(kViolet+6);
  /* g->SetMarkerSize(0); */
  /* g->SetMarkerStyle(kOpenCircle); */
  g->SetFillStyle(3504);
  g->SetFillColor(kViolet+6);

  g->Draw("E3same");
  drawTheorySideLineAs(g, kViolet+8);
  g->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
  g->GetYaxis()->SetTitle("Total Charm Cross Secton (#mu b)");
  g->GetYaxis()->SetTitleOffset(1.15);
  g->GetXaxis()->SetTitleOffset(1.15);
  gPad->SetLogy();
  gPad->SetLogx();
  drawLatex(0.6,0.41,"-1 < #eta < 1",0.06);
  drawLatex(0.6,0.34,"0 < p_{T} < 20 GeV/c",0.06);
  drawLatex(0.6,0.27,"p+p collisions",0.06);
  drawLatex(0.6,0.2,"FONLL calculation",0.06);
  g->SaveAs("TotCharmFonll.root");
}
