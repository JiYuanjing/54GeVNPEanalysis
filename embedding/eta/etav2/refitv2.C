#include "../sPhenixStyle.h"
#include "TCollection.h"
void myFit(double *x, double *p)
{
  if ( x[0]<p[6] ) return p[0]+x[0]*p[1]+x[0]*x[0]*p[2]+x[0]*x[0]*x[0]*p[3]+pow(x[0],4)*p[4];
  if (x[0]>=p[6]) return  p[0]+p[6]*p[1]+p[6]*p[6]*p[2]+pow(p[6],3)*p[3]+pow(p[6],4)*p[4]-p[6]*p[5]+x[0]*p[5];
}
//fit the kaon v2
void refitv2()
{
  SetsPhenixStyle();
  gStyle->SetOptFit(111);
  TCanvas* c = new TCanvas("c","c");
  TPDF* pdf = new TPDF("drawv2.pdf");
  // TFile* file = new TFile("Pionv2.root");
  TFile* file = new TFile("Kaonv2.root");
  pdf->Off();

  //0-10%
  int const nCentbins = 3;
  TString name[nCentbins] = {"0_10","10_40","40_80"};
  TGraphErrors* g200[nCentbins];
  TGraphErrors* g200Low[nCentbins];
  TMultiGraph* g200Combine[nCentbins];
  // TGraphErrors* g62;
  // g62 = (TGraphErrors*)file->Get("pionplus_0_10_62");
  // TF1* fit = new TF1("fit","pol5",0,7);
  TF1* fit[6];
  TF1* fitpol = new TF1("fitpol","pol4",0,18);
  TFile* fout = new TFile("fitKaonv2.root","recreate");
  // c->Divide(3,2);
  for (int i=0;i<nCentbins;i++)
  {
    fit[i] = new TF1(Form("fit_%s",name[i].Data()),myFit,0,18,7);
    fit[i]->SetNpx(10000);
    c->Clear();
    // c->cd(i+1);
    g200Combine[i]=new TMultiGraph(Form("kaonv2Com_%s",name[i].Data()),Form("kaonv2Com_%s",name[i].Data()));

    // g200Combine[i]=new TMultiGraph;
    g200Low[i] = (TGraphErrors*)file->Get(Form("K0s_%s_low",name[i].Data()));
    g200[i] = (TGraphErrors*)file->Get(Form("K0s_%s_high",name[i].Data())); 
    g200[i]->Draw();
    g200[i]->GetYaxis()->SetRangeUser(0,0.45);
    g200Combine[i]->Add(g200[i]);
    g200Low[i]->SetLineColor(kBlue);
    g200Low[i]->SetMarkerColor(kBlue);
    // if (i>0) {
      // g200Low[i]->Draw("psame");
      g200Combine[i]->Add(g200Low[i]);
    // }
    // if (i==0) {
    //   g62->SetMarkerColor(kBlue);
    //   g62->SetLineColor(kBlue);
    //  // g62->Draw("psame"); 
    //  g200Combine[i]->Add(g62);
    // }
    
    g200Combine[i]->Draw("p");
    // g200Combine[i]->GetXaxis()->SetRangeUser(0,1);
    g200[i]->GetXaxis()->SetTitle("p_{T}");
    g200[i]->GetYaxis()->SetTitle("v_{2}");

    g200Combine[i]->Fit(fitpol);
    double par[7];
    par[6]=6;
    fitpol->GetParameters(par);
    fit[i]->SetParameters(par);
    if (i==0)
    {
      fit[i]->FixParameter(6,2.6);
      fit[i]->FixParameter(5,0);
    }
    else if (i==1 )
    {   
      fit[i]->FixParameter(6,2.7);
      fit[i]->SetParLimits(5,-0.007,0.01);
    }
    else if (i==2 )   fit[i]->FixParameter(6,3);
    // else if (i==3)   fit[i]->FixParameter(6,3);
    // else if (i==4)   fit[i]->FixParameter(6,3);
    g200Combine[i]->Fit( fit[i],"B");
    fit[i]->SetLineColor(kMagenta);
    fit[i]->Draw("same");
    drawLatex(0.2,0.8,Form("%s \%",name[i].Data()),0.05);
    addpdf(pdf);
    g200Combine[i]->Write();
    fit[i]->Write();
  } 
 
    // addpdf(pdf);
  pdf->On();
  pdf->Close(); 
   
   
  

}
