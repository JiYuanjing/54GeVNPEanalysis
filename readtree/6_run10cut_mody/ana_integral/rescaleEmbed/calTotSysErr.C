#include "sPhenixStyle.h"
void calTotSysErr()
{
  SetsPhenixStyle();
  // int const nItems = 3;
  // double sys[nItems] = {0.42,0.18,1.5};
    // for (int i=0;i<nItems;i++)
  // {
  //   totalsys=sys[i]*sys[i];
  // }
  //totalsys=sqrt(totalsys);
  int const nbins = 10;
  double ptedge[nbins+1]={0.2,0.4,0.6,0.7,0.85,1,1.2,1.6,2.0,2.4,2.8};
  double totalsys[nbins],x[nbins];
  // double partPtCutSys = 0.0084;
  double partPtCutSys = 0.00276299;  //this value is after full gamma sample

  TFile* fCombineSys = new TFile("output/fCombineSys_30.root");
  TGraph* gCombineSys = (TGraph*)fCombineSys->Get("gCombineSys");
  fCombineSys->Close();

  TFile* fsys = TFile::Open("sysFromEmbedd.root");
  TGraph* gInvMassSys = (TGraph*)fsys->Get("gInvMassSys");
  TGraph* gPairDca = (TGraph*)fsys->Get("gPairDca");
  TGraph* gPartEnHitsSys = (TGraph*)fsys->Get("gPartEnHitsSys");
  
  for (int i=0;i<nbins;i++)
  { 
    x[i]=0.5*(ptedge[i]+ptedge[i+1]);
    totalsys[i] = 0;
    totalsys[i] += gPartEnHitsSys->Eval(x[i])*gPartEnHitsSys->Eval(x[i]);
    totalsys[i] += gPairDca->Eval(x[i])*gPairDca->Eval(x[i]);
    totalsys[i] += gInvMassSys->Eval(x[i])*gInvMassSys->Eval(x[i]);
    totalsys[i] += partPtCutSys*partPtCutSys;
    totalsys[i] += gCombineSys->Eval(x[i])*gCombineSys->Eval(x[i]); //not consider this part
    totalsys[i]=sqrt(totalsys[i]);
    cout <<x[i]<<" "<<totalsys[i]<< endl;
  }


  TGraph* gtot = new TGraph(nbins,x,totalsys);
  gtot->GetXaxis()->SetTitle("Tagged e p_{T} (GeV/c)");
  gtot->GetYaxis()->SetTitle("Sys. Error of Reco. Eff.");
  gtot->Draw("PAL");
  gtot->GetYaxis()->SetRangeUser(0.0,0.048);
  gInvMassSys->SetMarkerColor(kBlue);
  gInvMassSys->SetLineColor(kBlue);
  gInvMassSys->Draw("pLsame");
  gPairDca->SetMarkerColor(kRed);
  gPairDca->SetLineColor(kRed);
  gPairDca->Draw("pLsame");
  gPartEnHitsSys->SetMarkerColor(kGreen+2);
  gPartEnHitsSys->SetLineColor(kGreen+2);
  gPartEnHitsSys->Draw("pLsame");
  TLine* line = new TLine(0.1, partPtCutSys, 2.8, partPtCutSys);
  line->SetLineColor(kMagenta);
  line->SetLineWidth(2);
  line->Draw("same");
  gCombineSys->SetMarkerColor(kGray);
  gCombineSys->SetLineColor(kGray);
  gCombineSys->Draw("pLsame");

  TLegend* leg = new TLegend(0.2,0.65,0.5,0.9);
  leg->AddEntry(gtot,"Total sys. error","pe" );
  leg->AddEntry(gCombineSys, "Conversion probability","lep");
  leg->AddEntry(gInvMassSys,"Inv. mass","pe" );
  leg->AddEntry(gPairDca,"Pair DCA","pe" );
  leg->AddEntry(gPartEnHitsSys, "nHitsFit","pe" );
  leg->AddEntry(line, "Part. e p_{T} cut","l");
  leg->Draw("same");
  gtot->SaveAs("gTotSysErr.root");
  gPad->SaveAs("totsyserr.png");
  gPad->SaveAs("totsyserr.pdf");
  /* TGraph* gtot = new TGraph(nbins,x,totalsys); */
  /* gtot->GetXaxis()->SetTitle("Tag e p_{T}"); */
  /* gtot->Draw("PA");  */
  /* gtot->SaveAs("gTotSysErr.root"); */
 
}
