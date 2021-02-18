void calNonFlow()
{
  TFile *file = new TFile("profile_nowt.root");
  TProfile* pM2 = (TProfile*)file->Get("pM2_full");
  TProfile* pH = (TProfile*)file->Get("pHadron_full");
  double ptedge[14]={0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4};
  TH1F* hv2 = new TH1F("hv2","hv2",13,ptedge);;
  for (int i=0;i<hv2->GetNbinsX();i++)
  {
    float tmp=pM2->GetBinContent(i+1);
    tmp = tmp/0.046/242*pH->GetBinContent(i+1);
    /* tmp = tmp/0.046/242.; */
    cout <<tmp << endl;
    hv2->SetBinContent(i+1, tmp);
    hv2->SetBinError(i+1, pM2->GetBinError(i+1)/0.046/242.*pH->GetBinContent(i+1));
  }
  hv2->Draw();
}
