void draw()
{
  TFile* file = new TFile("embeddQa_gamma_Apr20.root");
  TH1F* h1 = (TH1F*)((TH2F*)file->Get("hTagElectron"))->ProjectionX("h1",1,9);
  TH1F* h2 = (TH1F*)((TH2F*)file->Get("hTagElectronPassCut"))->ProjectionX("h2",1,9);
  TH1F* h7 = (TH1F*)file->Get("htest7");
  TH1F* h8 = (TH1F*)file->Get("htest8");

  // h2->Divide(h1);
  // h2->Draw();
  TEfficiency* eff = new TEfficiency(*h2,*h1);
  // eff->Draw();
  TH1F* h3 = (TH1F*)h1->Clone("h3");
  for (int ip=0;ip<h2->GetNbinsX();ip++)
  {
    // h3->SetBinContent(ip+1, eff->GetEfficiency(ip+1));
    // h3->SetBinError(ip+1, eff->GetEfficiencyErrorLow(ip+1));
  }
  h3->SetMarkerStyle(22);
  h3->SetMarkerSize(1.2);
  h3->SetLineColor(kRed);
  h3->SetMarkerColor(kRed);
  h3->Draw("samep");
  // h8->Divide(h7);
  h7->Draw("same");
}
