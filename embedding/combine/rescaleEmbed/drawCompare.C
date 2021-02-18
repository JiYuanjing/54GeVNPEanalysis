#include "myStyle.h"

void drawCompare()
{
  SetMyStyle();
  TFile* f1 = new TFile("RecoEff_comb.root");
  TH1F* h1 = (TH1F*)f1->Get("hRecoEff_2_8");
  h1->SetDirectory(0);

  // TFile* f2 = new TFile("combineEff.root"); 
  TFile* f2 = new TFile("PhoEff_gamma.root"); 
  // TFile* f2 = new TFile("PhoEff_dirpho.root"); 
  TH1F* h2 = (TH1F*)f2->Get("hRecoEff_2_8");
  h2->SetDirectory(0);
  h1->SetMarkerColor(kRed);
  h2->SetMarkerColor(kBlue);
  h1->Draw();
  h2->Draw("same");

  TLegend* leg = new TLegend(0.65,0.7,0.88,0.88);
  leg->AddEntry(h1,"No weight","pl");
  leg->AddEntry(h2,"Combined","pl");
  
}
