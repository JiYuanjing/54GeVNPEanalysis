void checkTagEdistribution()
{
  TH1F* htest[7]; 

  int color[7] ={kRed,kOrange,kYellow,kGreen+2,kBlue,kCyan,kMagenta};
  
  // TFile* file = new TFile("rescale_combine.root");
  TFile* file = new TFile("rescaleFile/rescale_gamma.root");
  TH2F* h2reco = (TH2F*)file->Get("hTagElectron");
  TH1F* h1reco = (TH1F*)h2reco->ProjectionX();
  TH2F* h2pass = (TH2F*)file->Get("hTagElectronPassCut");
  TH1F* h1pass = (TH1F*)h2pass->ProjectionX();
  
  TLegend* leg = new TLegend(0.1,0.1,0.95,0.95);
  TCanvas* c = new TCanvas("c","c",1200,600);
  gPad->SetLogy();
  h1reco->Draw();
  h1reco->GetXaxis()->SetRangeUser(0,2.3);
  h1pass->SetLineColor(kRed);
  h1pass->Draw("same");
  gPad->SetGridx();


}
