void norm()
{
  TFile* f = new TFile("embeddQa_gamma_Apr16.root");
  TH2F* h2pi0 = (TH2F*)f->Get("hPi0Pt");
  TH1F* hcent = (TH1F*)f->Get("hEvent");
  TH1F* hpi0cent = (TH1F*)h2pi0->ProjectionY();
  hpi0cent->Divide(hcent);
  hpi0cent->Draw();
  for (int i=0;i<hpi0cent->GetNbinsX();i++)
  {
    cout<<hpi0cent->GetBinContent(i+1)<<endl;;
  }
}
