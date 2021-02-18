void sample()
{
  gRandom->SetSeed(0);
  TFile* f2 = new TFile("fout_pi0_eta_gamma1016.root");
  TF1* fun2 = (TF1*)f2->Get("fGMSp_pi0_8");
  TFile* f = new TFile("fout_pi0_eta_gamma_0926.root");
  TF1* fun = (TF1*)f->Get("fGMSp_pi0_8");
  TH1F* h = new TH1F("h","h",3000,0,6);
  int num=1e7;
  for (int i=0;i<num;i++)
  {
    double x = gRandom->Uniform()*6;
    double y = fun->Eval(x);
    h->Fill(x,y);
    if (i%1000000==0) cout<<i<<endl;
  }

  h->Draw();
  h->Scale(1./h->GetBinWidth(1)/(1.0*num)*6.0);
  fun->Draw("same");
  fun2->SetLineColor(kBlue);
  fun2->Draw("same");

}
