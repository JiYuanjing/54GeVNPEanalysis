#include "PID.h"
#include "../sPhenixStyle.h"
void fitfun(double *x, double * p)
{
  /* if (1!=-1*exp(x[0]/p[2])) return p[1]/(1.0+exp((x[0]-p[0])/p[2])); */
  /* return p[0]+p[1]*TMath::Log(x+p[2]); */
  /* return p[0]*x[0]+p[1]*x[0]*x[0]+p[2]*pow(x[0],3)+p[3]*pow(x[0],4); */
  double xv = x[0];
  double v2;
  if (p[2]!=0) 
  {
    /* v2 = p[0]*xv/(1.+TMath::Exp(-(xv-p[1])/p[2]))-p[3]*xv; */
    /* v2 = p[0]/(1.+TMath::Exp(-(xv-p[1])/p[2]))-p[3]; */
    v2 = p[0]/(1.+TMath::Exp(-(xv-p[1])/p[2]))-p[0]/(1.+TMath::Exp((p[1])/p[2]));
  }
  else v2 = 0;
  return v2;

}
void fitMt()
{
  SetsPhenixStyle();
  int centbin = 0;
  TString centbinname[4] = {"0-80%","0-10%","10-40%","40-80%"};
  //pi
  if (centbin==0){ 
    TGraphErrors* gPiPlus54 = new TGraphErrors( 14, MeanPT_pionPlus_0, val_pion_54_0,0,err_pion_54_0);
    gPiPlus54->SetName("gPiPlus54");
    TGraphErrors* gPiMinus54 = new TGraphErrors( 14, MeanPT_pionMinus_0, val_antipion_54_0,0,err_antipion_54_0);
    gPiMinus54->SetName("gPiMinus54");
    //K
    TGraphErrors* gKaonPlus54 = new TGraphErrors(14,MeanPT_KPlus_0,val_kaon_54_0,0,err_kaon_54_0);
    gKaonPlus54->SetName("gKaonPlus54");
    TGraphErrors* gKaonMinus54 = new TGraphErrors(14,MeanPT_KMinus_0,val_antikaon_54_0,0,err_kaon_54_0);
    gKaonMinus54->SetName("gKaonMinus54");
    //phi
    TGraphErrors* gPhi54= new TGraphErrors(11, MeanPT_phi_54_0, val_phi_54_0, 0, err_phi_54_0);
    gPhi54->SetName("gPhi54");
  }
  else if (centbin==1)
  {
    TGraphErrors* gPiPlus54 = new TGraphErrors( 14, MeanPT_pionPlus_1, val_pion_54_1,0,err_pion_54_1);
    gPiPlus54->SetName("gPiPlus54");
    TGraphErrors* gPiMinus54 = new TGraphErrors( 14, MeanPT_pionMinus_1, val_antipion_54_1,0,err_antipion_54_1);
    gPiMinus54->SetName("gPiMinus54");
    //K
    TGraphErrors* gKaonPlus54 = new TGraphErrors(14,MeanPT_KPlus_1,val_kaon_54_1,0,err_kaon_54_1);
    gKaonPlus54->SetName("gKaonPlus54");
    TGraphErrors* gKaonMinus54 = new TGraphErrors(14,MeanPT_KMinus_1,val_antikaon_54_1,0,err_kaon_54_1);
    gKaonMinus54->SetName("gKaonMinus54");
    //phi
    TGraphErrors* gPhi54= new TGraphErrors(11, MeanPT_phi_54_1, val_phi_54_1, 0, err_phi_54_1);
    gPhi54->SetName("gPhi54");

  }
  else if (centbin==2)
  {
    TGraphErrors* gPiPlus54 = new TGraphErrors( 14, MeanPT_pionPlus_2, val_pion_54_2,0,err_pion_54_2);
    gPiPlus54->SetName("gPiPlus54");
    TGraphErrors* gPiMinus54 = new TGraphErrors( 14, MeanPT_pionMinus_2, val_antipion_54_2,0,err_antipion_54_2);
    gPiMinus54->SetName("gPiMinus54");
    //K
    TGraphErrors* gKaonPlus54 = new TGraphErrors(14,MeanPT_KPlus_2,val_kaon_54_2,0,err_kaon_54_2);
    gKaonPlus54->SetName("gKaonPlus54");
    TGraphErrors* gKaonMinus54 = new TGraphErrors(14,MeanPT_KMinus_2,val_antikaon_54_2,0,err_kaon_54_2);
    gKaonMinus54->SetName("gKaonMinus54");
    //phi
    TGraphErrors* gPhi54= new TGraphErrors(11, MeanPT_phi_54_2, val_phi_54_2, 0, err_phi_54_2);
    gPhi54->SetName("gPhi54");

  }
  else if (centbin==3)
  {
    TGraphErrors* gPiPlus54 = new TGraphErrors( 14, MeanPT_pionPlus_3, val_pion_54_3,0,err_pion_54_3);
    gPiPlus54->SetName("gPiPlus54");
    TGraphErrors* gPiMinus54 = new TGraphErrors( 14, MeanPT_pionMinus_3, val_antipion_54_3,0,err_antipion_54_3);
    gPiMinus54->SetName("gPiMinus54");
    //K
    TGraphErrors* gKaonPlus54 = new TGraphErrors(14,MeanPT_KPlus_3,val_kaon_54_3,0,err_kaon_54_3);
    gKaonPlus54->SetName("gKaonPlus54");
    TGraphErrors* gKaonMinus54 = new TGraphErrors(14,MeanPT_KMinus_3,val_antikaon_54_3,0,err_kaon_54_3);
    gKaonMinus54->SetName("gKaonMinus54");
    //phi
    TGraphErrors* gPhi54= new TGraphErrors(11, MeanPT_phi_54_3, val_phi_54_3, 0, err_phi_54_3);
    gPhi54->SetName("gPhi54");
  }
  else { cout<<"please input correct cent bin 0-3" <<endl;return;}
  TString name[3]={"#pi^{+}","K^{+}","#phi"};
  int color[3]={kRed, kCyan+2, kMagenta};
  int style[3]={27,28,25};
  TGraphErrors* gPiPlus54Mt = (TGraphErrors*)ptToNCQ(gPiPlus54,2.,0.13957,"gPiPlus54Mt",style[0],color[0],2.4);
  TGraphErrors* gPiMinus54Mt =(TGraphErrors*)ptToNCQ(gPiMinus54,2.,0.13957,"gPiMinus54Mt",style[0],color[0],2.4);
  TGraphErrors* gKaonPlus54Mt= (TGraphErrors*)ptToNCQ(gKaonPlus54,2.,0.493677,"gKaonPlus54Mt",style[1],color[1],2.4);
  TGraphErrors* gPhi54Mt= (TGraphErrors*)ptToNCQ(gPhi54,2.,1.019461,"gPhi54Mt",style[2],color[2],1.8);
  TMultiGraph* gmult = new TMultiGraph("mult","mult");
  gmult->Add(gPhi54Mt);
  gmult->Add(gPiPlus54Mt);
  gmult->Add(gKaonPlus54Mt);

  gPhi54->SetMarkerColor(color[2]);
  gPhi54->SetMarkerStyle(style[2]);
  gPhi54->SetMarkerSize(1.8);
  gPiPlus54->SetMarkerColor(color[0]);
  gPiPlus54->SetMarkerStyle(style[0]);
  gPiPlus54->SetMarkerSize(2.4);
  gKaonPlus54->SetMarkerColor(color[1]);
  gKaonPlus54->SetMarkerStyle(style[1]);
  gKaonPlus54->SetMarkerSize(2.4);

  TH1D* h = new TH1D("h","h",1,0,1.5);
  h->Draw();
  gmult->Draw("samep");
  gmult->GetXaxis()->SetTitle("(m_{T}-m_{0})/n_{q} (GeV/c^{2})");
  gmult->GetYaxis()->SetTitle("v_{2}/n_{q}");

  TF1* f2 = new TF1("fun2", "pol4", 0,2.5);
  TF1* f = new TF1("fun", fitfun, 0,2.5,4);
  f->SetParameters(1,0.5,4,1);
  f->SetLineStyle(2);
  /* gmult->Fit("pol3"); */
  gmult->Fit(f);
  TLegend* leg = new TLegend(0.7,0.2,0.88,0.6);
  leg->SetTextSize(0.055);
  leg->AddEntry(gPiPlus54,"#pi^{+}","p");
  leg->AddEntry(gKaonPlus54,"K^{+}","p");
  leg->AddEntry(gPhi54,"#phi","p");
  leg->AddEntry(f,"fit","l");
  leg->Draw();
  drawLatex(0.2,0.78,centbinname[centbin].Data(),0.055);
  drawLatex(0.2,0.85,"Au+Au 54.4 GeV",0.055);
  drawLatex(0.3,0.2,"data: STAR preliminary",0.045,kGray+2);
  /* TLatex* lat; */
  /* lat->DrawLatexNDC(0.2,0.88,centbinname[centbin].Data()); */
  f->SaveAs(Form("fitPIDV2_%d.root",centbin));
}

TGraphErrors* ptToNCQ(TGraphErrors* g, double nq,double m0,TString name,int style,int color,double size)
{
  double x[50],y[50],err[50],sys[50];
  for (int ip=0;ip<g->GetN();ip++)
  {
    g->GetPoint(ip,x[ip],y[ip]);
    y[ip]=y[ip]/nq;
    x[ip]=(sqrt(x[ip]*x[ip]+m0*m0)-m0)/nq;
    err[ip] = g->GetErrorY(ip)/nq;
  }

  TGraphErrors* g = new TGraphErrors(ip,x,y,0,err);
  g->SetName(name.Data());
  g->SetMarkerStyle(style);
  g->SetMarkerColor(color);
  g->SetMarkerSize(size);
  return g;
}
