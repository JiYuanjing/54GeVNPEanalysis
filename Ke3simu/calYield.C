#include "sPhenixStyle.h"

void scale(TH1D* h, double etarange, double branch, double nevents);
void calspectra(TString histname, TString infilename, TFile* fout, TString name);
void calD2e();
void calKe3Frac();
void checkSpectra();
void projV2();
Double_t fitfun(double *x, double * p)
{
  double xv = x[0];
  double v2;
  if (p[2]!=0)
  {
    v2 = p[0]/(1.+TMath::Exp(-(xv-p[1])/p[2]))-p[0]/(1.+TMath::Exp((p[1])/p[2]));
  }
  else v2 = 0;
  return (Double_t)v2;
}
double myFit(double *x, double *p)
{
  if ( x[0]<p[6] ) return p[0]+x[0]*p[1]+x[0]*x[0]*p[2]+x[0]*x[0]*x[0]*p[3]+pow(x[0],4)*p[4];
  else if (x[0]>=p[6]) return  p[0]+p[6]*p[1]+p[6]*p[6]*p[2]+pow(p[6],3)*p[3]+pow(p[6],4)*p[4]-p[6]*p[5]+x[0]*p[5];
  else return 0;
}

void calYield()
{
  SetsPhenixStyle();
  TFile* fout;
  fout = new TFile( "foutK.root", "recreate");
  calspectra("hKaonSpec","Kplus.root",fout,"Kaon");
  calspectra("hElectronSpecDcaCut","Kplus.root",fout, "KplusDcaCut");
  calspectra("hElectronSpecNoCut","Kplus.root",fout, "KplusNoCut");
  calspectra("hElectronSpecDLCut","KplusTest.root",fout, "KplusDLCut");
  /* calspectra("hElectronSpecDcaCut","K0L.root",fout, "K0LDcaCut"); */
  calspectra("hElectronSpecDcaCut","Klong.root",fout, "K0LDcaCut");
  calspectra("hElectronSpecNoCut","Klong.root",fout, "K0LNoCut");
  /* calspectra("hElectronSpecDcaCut","KL_test.root",fout, "K0LDcaCut"); */
  calspectra("hElectronSpecDcaCut","K0s.root",fout, "K0sDcaCut");
  calspectra("hElectronSpecNoCut","K0s.root",fout, "K0sNoCut");
  fout->Close();

  checkSpectra();
  calD2e();
  calKe3Frac();
  /* projV2(); */
}

void calspectra(TString histname, TString infilename,  TFile* fout, TString name)
{
  TFile* f = TFile::Open(infilename.Data());
  TH3D* h3sp = (TH3D*)f->Get(histname.Data());
  h3sp->Sumw2();
  TH1D* hCent = (TH1D*)f->Get("hCent");
  //project 0-60%
  double etarange;
  TH1D* hsp_0_60 = (TH1D*)h3sp->ProjectionX(Form("h%ssp_0_60",name.Data()), h3sp->GetYaxis()->FindBin(-0.5+1e-6), h3sp->GetYaxis()->FindBin(0.5-1e-6), 
      h3sp->GetZaxis()->FindBin(2+1e-6), h3sp->GetZaxis()->FindBin(8-1e-6)); 
  hsp_0_60->SetDirectory(0);
  hsp_0_60->Rebin(5);
  hsp_0_60->Draw();

  scale(hsp_0_60, 1, 1,hCent->Integral(hCent->GetXaxis()->FindBin(2),  hCent->GetXaxis()->FindBin(8) ) );
  hsp_0_60->Draw();
  /* hsp_0_60->GetBinContent(hsp_0_60->GetXaxis()->FindBin()); */
  int cent[10]={80,70,60,50,40,30,20,10,5,0};
  TH1D* hsp[9]; 
  for (int i=2;i<9;i++){
    hsp[i] = (TH1D*)h3sp->ProjectionX(Form("h%ssp_%d_%d", name.Data(), cent[i+1], cent[i]), h3sp->GetYaxis()->FindBin(-0.5+1e-6), h3sp->GetYaxis()->FindBin(0.5-1e-6), 
                                             h3sp->GetZaxis()->FindBin(i), h3sp->GetZaxis()->FindBin(i)); 
    cout<<i<<" " << Form("hsp_%d_%d", cent[i+1], cent[i])<<" "<< hCent->GetBinContent(hCent->GetXaxis()->FindBin(i+1e-6))<<endl;
    hsp[i]->SetDirectory(0);
    hsp[i]->Rebin(5);
    scale(hsp[i], 1, 1,hCent->GetBinContent(hCent->GetXaxis()->FindBin(i+1e-6)) );
  }
  /* hsp[5]->Draw(); */
  /* TGraphAsymmErrors* g = new TGraphAsymmErrors("data/Ks/20_30.txt", "%lg %lg %lg %lg"); */
  /* g->SetMarkerStyle(kOpenCircle); */
  /* g->SetMarkerColor(kBlack); */
  /* g->SetMarkerSize(1.4); */
  /* g->Draw("psame"); */
  /* cout<<g->Eval(0.4)<<endl;; */

  f->Close();
  /*  */
  fout->cd();
  hsp_0_60->Write();
  for (int i=2;i<9;i++) hsp[i]->Write(); 

}
void scale(TH1D* h, double etarange, double branch, double nevents)
{
  double weight = etarange*h->GetXaxis()->GetBinWidth(1)*2*3.1415*nevents;
  for (int i=1;i<=h->GetNbinsX(); i++)
  {
    h->SetBinContent(i, h->GetBinContent(i)/weight/h->GetBinCenter(i)*branch);
    h->SetBinError(i, h->GetBinError(i)/weight/h->GetBinCenter(i)*branch);
  }
}
void calD2e()
{
  TFile* file = new TFile("fnollD2e.root","recreate"); 
  TGraph* g = new TGraph("D2e.txt","%lg %lg");
  /* TGraph* g = new TGraph("D2e_old.txt","%lg %lg"); */
  for (int i=0;i<g->GetN();i++){
    /* g->SetPointY(i, g->GetPointY(i)/(36*1e9)/2./2/3.14159/g->GetPointX(i)  ); //FONLL is ds/dpt, divide 2pi*dy*pt */
    g->SetPointY(i, g->GetPointY(i)/(36*1e9)); //FONLL is ds/dpt, divide 2pi*dy*pt
  }
  g->SetName("gFonllD2e");
  g->Write();
  file->Close();
}

void calKe3Frac()
{
  TFile* fout = TFile::Open("fnollD2e.root");
  TGraph* gfnoll = (TGraph*)fout->Get("gFonllD2e");
  fout->Close();

  TFile* fKe = TFile::Open("foutK.root");
  TH1F* hKplusCut = (TH1F*)fKe->Get("hKplusDcaCutsp_0_60");
  hKplusCut->SetDirectory(0);

  TH1F* hKplusDLCut = (TH1F*)fKe->Get("hKplusDLCutsp_0_60");
  hKplusDLCut->SetDirectory(0);

  TH1F* hKplus = (TH1F*)fKe->Get("hKplusNoCutsp_0_60");
  hKplus->SetDirectory(0);

  TH1F* hK0LCut = (TH1F*)fKe->Get("hK0LDcaCutsp_0_60");
  hK0LCut->SetDirectory(0);
  TH1F* hK0L = (TH1F*)fKe->Get("hK0LNoCutsp_0_60");
  hK0L->SetDirectory(0);

  TH1F* hK0sCut = (TH1F*)fKe->Get("hK0sDcaCutsp_0_60");
  hK0sCut->SetDirectory(0);
  TH1F* hK0s= (TH1F*)fKe->Get("hK0sNoCutsp_0_60");
  hK0s->SetDirectory(0);
  fKe->Close();

  TH1F* hEffks = (TH1F*)hK0sCut->Clone("hKsEff");
  hEffks->SetDirectory(0);
  hEffks->Divide(hK0s);
  /* hEffks->GetYaxis()->SetTitle("K^{0}_{s}#rightarrowe |gDca|<1.5 / Total"); */
  hEffks->GetYaxis()->SetTitle("K#rightarrowe |gDca|<1.5 / Total");
  hEffks->SetMarkerColor(kGreen+2);
  hEffks->Scale(0.01);
  hEffks->Draw();
  hEffks->GetYaxis()->SetRangeUser(0,0.03);
  TH1F* hEff = (TH1F*)hKplusCut->Clone("hKplusEff");
  hEff->SetDirectory(0);
  hEff->Divide(hKplus);
  hEff->GetYaxis()->SetTitle("K^{+}#rightarrowe |gDca|<1.5 / Total");
  hEff->SetMarkerColor(kRed-3);
  hEff->Draw("same");
  hEff->GetYaxis()->SetRangeUser(0,1e-1);

  TH1F* hEffkpDl = (TH1F*)hKplusCut->Clone("hKplusDLEff");
  hEffkpDl->SetDirectory(0);
  hEffkpDl->Divide(hKplusDLCut);
  hEffkpDl->Scale(0.1);
  hEffkpDl->GetYaxis()->SetTitle("K^{+}#rightarrowe |gDca|<1.5 / Decay in TPC");
  hEffkpDl->SetMarkerColor(kViolet);
  hEffkpDl->SetMarkerStyle(29);
  hEffkpDl->Draw("same");
  hEffkpDl->GetYaxis()->SetRangeUser(0,1e-1);

  /* gPad->SaveAs("plot/EffKp.pdf"); */
  TH1F* hEffkl = (TH1F*)hK0LCut->Clone("hKlongEff");
  hEffkl->SetDirectory(0);
  hEffkl->Divide(hK0L);
  hEffkl->GetYaxis()->SetTitle("K^{0}_{L}#rightarrowe |gDca|<1.5 / Total");
  hEffkl->SetMarkerColor(kBlue+2);
  hEffkl->Draw("same");
  hEffkl->GetYaxis()->SetRangeUser(0,1e-1);
  /* gPad->SaveAs("plot/Effkl.pdf"); */

  /* gPad->SaveAs("plot/Effks.pdf"); */
  TLegend* legEff = new TLegend(0.2,0.7,0.88,0.92);
  legEff->AddEntry(hEff, "K^{+}", "lp");
  legEff->AddEntry(hEffkl, "K^{0}_{L}", "lp");
  legEff->AddEntry(hEffks, "K^{0}_{S} #times 0.01", "lp");
  legEff->AddEntry(hEffkpDl, "K^{+} decay in TPC #times 0.1", "lp");
  legEff->Draw();
  gPad->SaveAs("plot/Eff.pdf");

  hKplusCut->Scale(2.*0.05); //add K+ + K- 
  TH1F* ratio = (TH1F*)hKplusCut->Clone("hRatio");
  // add k0L fraction
  hK0LCut->Scale(0.4); // branch ratio
  ratio->Add(hK0LCut);

  hK0sCut->Scale(7e-4); // Ks
  ratio->Add(hK0sCut);

  for (int ib=1;ib<ratio->GetNbinsX()+1;ib++)
  {
    if ((ratio->GetBinContent(ib)+gfnoll->Eval(ratio->GetBinCenter(ib))*2*335.28)>0)
    { 
      /* hK0LCut->SetBinContent( ib, hK0LCut->GetBinContent(ib)/(ratio->GetBinContent(ib)+gfnoll->Eval(hK0LCut->GetBinCenter(ib))*2*335)); */
      /* hK0sCut->SetBinContent( ib, hK0sCut->GetBinContent(ib)/(ratio->GetBinContent(ib)+gfnoll->Eval(hK0sCut->GetBinCenter(ib))*2*335)); */
      /* hKplusCut->SetBinContent( ib, hKplusCut->GetBinContent(ib)/(ratio->GetBinContent(ib)+gfnoll->Eval(hKplusCut->GetBinCenter(ib))*2*335)); */
      /* ratio->SetBinContent( ib, ratio->GetBinContent(ib)/(ratio->GetBinContent(ib)+gfnoll->Eval(ratio->GetBinCenter(ib))*2*335));  */
      hK0LCut->SetBinContent( ib, hK0LCut->GetBinContent(ib)/(gfnoll->Eval(hK0LCut->GetBinCenter(ib))*2*335.28));
      hK0sCut->SetBinContent( ib, hK0sCut->GetBinContent(ib)/(gfnoll->Eval(hK0sCut->GetBinCenter(ib))*2*335.28));
      hKplusCut->SetBinContent( ib, hKplusCut->GetBinContent(ib)/(gfnoll->Eval(hKplusCut->GetBinCenter(ib))*2*335.28));
      ratio->SetBinContent( ib, ratio->GetBinContent(ib)/(gfnoll->Eval(ratio->GetBinCenter(ib))*2*335.28)); 
    }
  }
  /* ratio->GetYaxis()->SetTitle("Ke3 / (Ke3+D#rightarrowe)"); */
  ratio->GetYaxis()->SetTitle("Ke3 / D#rightarrowe");
  hK0LCut->SetMarkerColor(kRed);
  hK0sCut->SetMarkerColor(kGreen+2);
  hKplusCut->SetMarkerColor(kBlue);
  ratio->Draw();
  hK0LCut->Draw("same");
  hK0sCut->Draw("same");
  hKplusCut->Draw("same");

  TLegend* leg = new TLegend(0.6,0.6,0.88,0.88);
  leg->AddEntry( hK0LCut,"K^{0}_{L}","lep" );
  leg->AddEntry( hKplusCut,"K^{+}+K^{-}","lep" );
  leg->AddEntry( hK0sCut,"K^{0}_{s}","lep" );
  leg->Draw();
  gPad->SaveAs("plot/ratio.pdf");

  TH1F* htest = (TH1F*)hK0LCut->Clone("htest");
  htest->Divide(hKplusCut);
  htest->Draw();
  gPad->SaveAs("test.pdf");
}

void checkSpectra()
{
  TFile* fKe = TFile::Open("foutK.root");
  TH1F* hKaon= (TH1F*)fKe->Get("hKaonsp_0_60");
  hKaon->SetDirectory(0);
  fKe->Close();

  TFile* f = new TFile("data/Ks/Ksnew.root");
  TGraphAsymmErrors* g = (TGraphAsymmErrors*)f->Get("0_60"); 
  g->SetMarkerSize(1.5);
  g->SetMarkerStyle(kFullCircle);
  g->SetMarkerColor(kRed);
  gPad->SetLogy();

  hKaon->Draw();
  hKaon->GetYaxis()->SetTitle("d^{2}N/2#pip_{T}dydp_{T}");
  g->Draw("same p");
  drawLatex(0.25, 0.25, "Kaon 0-60%", 0.045);
  gPad->SaveAs("plot/checkspectra.pdf");
  gPad->SetLogy(0);

}

void projV2()
{
  TFile* f = TFile::Open("Kplus.root");
  TProfile3D* pKplus3D = (TProfile3D*)f->Get("pElectronV2DcaCut");
  /* TProfile3D* pKplus3D = (TProfile3D*)f->Get("pElectronV2NoCut"); */
  /* TProfile3D* pKplus3D = (TProfile3D*)f->Get("pKaonV2"); */
  pKplus3D->GetYaxis()->SetRange(pKplus3D->GetYaxis()->FindBin(-0.8+1e-6), pKplus3D->GetYaxis()->FindBin(0.8-1e-6));
  TProfile2D* pKplus2D = (TProfile2D*)pKplus3D->Project3DProfile("zx");
  TProfile* pKplus = (TProfile*)pKplus2D->ProfileX("pKplus", 3, 9);
  pKplus->SetDirectory(0);
  f->Close();

  f = TFile::Open("Klong.root");
  TProfile3D* pK0L3D = (TProfile3D*)f->Get("pElectronV2DcaCut");
  TProfile2D* pK0L2D = (TProfile2D*)pK0L3D->Project3DProfile("zx");
  TProfile* pK0L = (TProfile*)pK0L2D->ProfileX("pK0L", pK0L3D->GetZaxis()->FindBin(2+1e-6), pK0L3D->GetZaxis()->FindBin(8-1e-6));
  pK0L->SetDirectory(0);

  f = TFile::Open("K0s.root");
  TProfile3D* pK0s3D = (TProfile3D*)f->Get("pElectronV2DcaCut");
  TProfile2D* pK0s2D = (TProfile2D*)pK0s3D->Project3DProfile("zx");
  TProfile* pK0s = (TProfile*)pK0s2D->ProfileX("pK0s", pK0s3D->GetZaxis()->FindBin(2+1e-6), pK0s3D->GetZaxis()->FindBin(8-1e-6));
  pK0s->SetDirectory(0);
  f->Close();

  pKplus->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  pKplus->GetYaxis()->SetTitle("K#rightarrowe v_{2}");
  /* pKplus->GetXaxis()->SetRangeUser(0,3); */

  pKplus->Rebin(5);
  pK0L->Rebin(5);
  pK0s->Rebin(5);

  pKplus->SetMarkerColor(kRed);
  pKplus->SetLineColor(kRed);
  pK0s->SetMarkerColor(kGreen+2);
  pK0s->SetLineColor(kGreen+2);
  pK0L->SetMarkerColor(kBlue+3);
  pK0L->SetLineColor(kBlue+3);

  pKplus->Draw();
  pK0L->Draw("same");
  pK0s->Draw("same");

  TProfile* pKaonV2 = (TProfile*)pK0L->Clone("pKaonV2");
  pKaonV2->Add(pKplus, pK0L, 0.05*2, 0.4);
  pKaonV2->Add(pK0s, 7e-4);
  pKaonV2->SetMarkerStyle(kFullSquare);
  pKaonV2->Draw("same");

  TF1* fun = new TF1("fitfun", fitfun, 0, 4, 3);
  /* TF1* fun = new TF1("fitfun", myFit, 0, 4, 7); */
  /* fun->FixParameter(6,2.); */
  fun->SetParameters(1,1,1);
  pKaonV2->Fit(fun);

  TLegend* legv2 = new TLegend(0.7,0.7,0.88,0.92);
  legv2->AddEntry(pKplus, "K^{+}#rightarrowe", "lp");
  legv2->AddEntry(pK0L, "K^{0}_{L}#rightarrowe", "lp");
  legv2->AddEntry(pK0s, "K^{0}_{S}#rightarrowe", "lp");
  legv2->AddEntry(pKaonV2, "Total", "lp");
  legv2->Draw();

  gPad->SaveAs("plot/Kaonv2.pdf");
  TFile *fv2 = new TFile("fKaonV2.root","recreate");
  pKaonV2->Write();
  pKplus->Write();
  pK0L->Write();
  pK0s->Write();
  fun->Write();
  fv2->Close();

  //////////////////////////////fit result/////////////////////////////////////
  //   NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE  //
  //   1  p0           2.18526e-01   2.17567e-03   2.56879e-06   2.60905e+00 //
  //   2  p1           3.00555e-01   4.59579e-03   4.07522e-05   4.48280e-01 //
  //   3  p2           4.54157e-01   5.13039e-03   6.08919e-06  -8.85434e-01 //
  /////////////////////////////////////////////////////////////////////////////
}
