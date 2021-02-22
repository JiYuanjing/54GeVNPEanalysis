void addpdf(TPDF* pdf,TCanvas* c)
{
  c->cd();
  pdf->On();
  pdf->NewPage();
  c->Update();
  pdf->Off();
}
void fitDCAZ()
{
  TFile* f = new TFile("DCA.root");
  TCanvas* c = new TCanvas("c","c");
  c->Divide(4,2);
  TPDF* pdf = new TPDF("test.pdf");
  pdf->Off();
  /* TH2F* hDCAxy = (TH2F*)f->Get("hDcaXY"); */
  TH2F* hDCAxy = (TH2F*)f->Get("hDcaZ");
  TH1F* hXYres = (TH1F*)hDcaXY->ProjectionY("hZres");
  /* int nBins = ; */
  /* double ptedge[]; */
  /* hDCAxy->Rebin(2); */
  TF1* fun = new TF1("fun","gaus",-2,2);
  for (int i =0;i<hDCAxy->GetYaxis()->GetNbins();i++)
  {
     c->cd((i)%8+1);
     TH1* h = hDCAxy->ProjectionX("h",i+1,i+1);
     h->DrawCopy();
     h->GetXaxis()->SetRangeUser(-0.2,0.2);
     h->Fit(fun,"R+");
     hXYres->SetBinContent(i+1,fun->GetParameter(2));
     /* if (((i+1)%8)==0){ */
     /*   addpdf(pdf,c); */
     /*   #<{(| c->Clear(); |)}># */
     /*   #<{(| c->Divide(4,2); |)}># */
     /* } */
  }
  c->cd();
  addpdf(pdf,c);
  
  TF1* fit2 = new TF1("fit2","sqrt([0]+[1]/x/x)",0.2,5);
  hXYres->Draw();
  hXYres->Fit(fit2);
  addpdf(pdf,c);
  pdf->On();
  pdf->Close();
  hXYres->SaveAs("fk2.root");
}
void fitDCA()
{
  TFile* f = new TFile("DCA.root");
  TCanvas* c = new TCanvas("c","c");
  c->Divide(4,2);
  TPDF* pdf = new TPDF("test.pdf");
  pdf->Off();
  TH2F* hDCAxy = (TH2F*)f->Get("hDcaXY");
  TH1F* hXYres = (TH1F*)hDcaXY->ProjectionY("hXYres");
  /* int nBins = ; */
  /* double ptedge[]; */
  /* hDCAxy->Rebin(2); */
  TF1* fun = new TF1("fun","gaus",-2,2);
  for (int i =0;i<hDCAxy->GetYaxis()->GetNbins();i++)
  {
     c->cd((i)%8+1);
     TH1* h = hDCAxy->ProjectionX(Form("p_{T}_%0.1fGeV",hDCAxy->GetYaxis()->GetBinCenter(i+1)),i+1,i+1);
     h->DrawCopy();
     h->GetXaxis()->SetRangeUser(-0.2,0.2);
     h->Fit(fun,"R+");
     h->GetXaxis()->SetRangeUser(-2,2);
     hXYres->SetBinContent(i+1,fun->GetParameter(2));

     if (((i+1)%8)==0){
       addpdf(pdf,c);
       /* c->Clear(); */
       /* c->Divide(4,2); */
     }
  }
  c->cd();
  TF1* fit2 = new TF1("fit2","sqrt([0]+[1]/x/x)",0.2,5);
  hXYres->Draw();
  hXYres->GetXaxis()->SetRangeUser(0.15,3);
  hXYres->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hXYres->GetYaxis()->SetTitle("DCA_{XY} Res (cm)");
  hXYres->Fit(fit2);
  addpdf(pdf,c);
  pdf->On();
  pdf->Close();
  hXYres->SaveAs("fk.root");
}
