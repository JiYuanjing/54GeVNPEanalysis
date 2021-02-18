#include "sPhenixStyle.h"
void projectMom(){
  SetsPhenixStyle();
  gStyle->SetPalette(1);
  // TFile* file = new TFile("PE.root");
  TCanvas* c = new TCanvas("c","c",600,600);
  c->cd();
  TFile* file = new TFile("hePtvsP.root");
  TH2F* hp_pt = (TH2F*)file->Get("hePtvsP");
  hp_pt->Draw("colz");
  gPad->SetGridx(1);
  gPad->SetGridy(1);

  gPad->SetLogz();
  hp_pt->GetYaxis()->SetTitle("p [GeV/c]");
  hp_pt->GetYaxis()->SetRangeUser(0,6);
  hp_pt->GetXaxis()->SetRangeUser(0,6);
  hp_pt->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  c1->cd();
  gPad->SetGridx(1);
  gPad->SetGridy(1);

  TProfile* p_pt = (TProfile*)hp_pt->ProfileY();
  p_pt->GetXaxis()->SetTitle("p [GeV/c]");
  p_pt->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  p_pt->GetXaxis()->SetRangeUser(0,6);
  p_pt->SetMarkerSize(0.6);
  p_pt->Draw();
  p_pt->SaveAs("p_pt_MB_54.root");
}
