//change the pi/K/p as the histogram fit
#include "rootlogon.h"
#include <string>
#include "Common.h"
TPDF* pdf;
TH1F* hesamp[400];
TH1F* hpsamp[400];
TH1F* hpisamp[400];
TGraph* gpisamp[400];
TH1F* hKsamp[400];
TH1F* htofsamp[400];
int ipar=0;
void addpdf(TPDF* pdf)
{
  pdf->On();
  pdf->NewPage();
  gPad->Update();
  pdf->Off(); 
}
void drawtitle(TPDF* pdf,TCanvas* c,string s){
  c->cd();
  c->Draw();
  // setPad(0.1,0.1,0.05,0.12);
  TLatex t;
  t.SetTextSize(0.05);
  t.DrawText(0.2,0.5,s.c_str());
  TLatex la;
  la.SetTextSize(0.035);
  la.DrawText(0.1,0.3,(new TDatime())->AsSQLString());
  la.DrawText(0.1,0.2,"by Yuanjing");
  pdf->On();
  pdf->NewPage();
  gPad->Update();
  pdf->Off();
  c->Clear();
}
void projectionAndFit(TH2F* h, float lowpt, float highpt, float &mean,float &sigma,
    float range1,float range2,float meanL,float meanH,float sigmaL,float sigmaH, string p,int centL,int centH)
{
  int lbin = h->GetXaxis()->FindBin(lowpt);
  int hbin = h->GetXaxis()->FindBin(highpt);
  TH1F* hx = (TH1F*)h->ProjectionY("hx",lbin,hbin);

  if (lowpt>1.5) 
  {
    hx->Rebin();
  }
  if (p=="#pi"){ 
    // hx->Rebin();
    if (lowpt>2) 
    {
      hx->Rebin();
    }
  }

  // hx->Smooth();
  gPad->SetLogy(1);
  hx->DrawCopy();
  TF1* fit  = new TF1("fit","gaus", range1,range2);
  hx->Fit(fit,"RN");
  fit->SetLineColor(kRed);
  mean = fit->GetParameter(1); 
  if (mean>meanH || mean <meanL) {
    fit->SetParLimits(2, meanL, meanH);
    hx->Fit(fit,"BRN");
  }
  sigma = fit->GetParameter(2); 
  if (sigma>sigmaH|| sigma<sigmaL) {
    fit->SetParLimits(2, sigmaL, sigmaH);
    hx->Fit(fit,"BRN");
  }
  fit->DrawCopy("same"); 
  mean = fit->GetParameter(1); 
  sigma = fit->GetParameter(2); 
  TLatex lat;
  lat.SetTextSize(0.035);
  lat.DrawLatexNDC(0.2,0.85,Form("%.2f<p_{T}<%.2f GeV %s",lowpt,highpt,p.c_str()));
  if (p.find("proton")!=std::string::npos) {
    float total = hx->Integral("width");
    hx->Scale(1./total); 
    hpsamp[ipar] = (TH1F*)hx->Clone(Form("h%s_%.1f_%.1f",p.c_str() ,lowpt,highpt));
    hpsamp[ipar]->SetDirectory(0);
    cout<< Form("h%s_%.3f_%.3f",p.c_str() ,lowpt,highpt)<<endl;
    // ipar++;
  }
  if (p.find("pi")!=std::string::npos) {
    float total = hx->Integral("width");
    hx->Scale(1./total); 
    hpisamp[ipar] = (TH1F*)hx->Clone(Form("h%s_%.1f_%.1f",p.c_str() ,lowpt,highpt));
    hpisamp[ipar]->SetDirectory(0);
    cout<< Form("h%s_%.3f_%.3f",p.c_str() ,lowpt,highpt)<<endl;
    // ipar++;
    gpisamp[ipar] = new TGraph(hpisamp[ipar]);
  }

  if (p.find("K")!=std::string::npos) {
    float total = hx->Integral("width");
    hx->Scale(1./total); 
    hKsamp[ipar] = (TH1F*)hx->Clone(Form("h%s_%.1f_%.1f",p.c_str() ,lowpt,highpt));
    hKsamp[ipar]->SetDirectory(0);
    cout<< Form("h%s_%.3f_%.3f",p.c_str() ,lowpt,highpt)<<endl;
    hKsamp[ipar]->Smooth();
    // ipar++;
  }
}
void projAndFitForMg(TH2F* h, float lowpt, float highpt, float &mean,float &sigma,
    float range1,float range2,float meanL,float meanH,float sigmaL,float sigmaH, string p,int centL,int centH)
{
  int lbin = h->GetXaxis()->FindBin(lowpt);
  int hbin = h->GetXaxis()->FindBin(highpt);

  TH1F* hx = (TH1F*)h->ProjectionY("hx",lbin,hbin);
  if (lowpt>2.2) hx->Rebin();
  //hx->SetDirectory(0);
  gPad->SetLogy(1);
  hx->DrawCopy();
  hx->GetXaxis()->SetRangeUser(2,6);
  int meanbin =  hx->GetMaximumBin();
  TF1* fit  = new TF1("fit","gaus",-5 ,15);
  fit->SetParLimits(2, sigmaL, sigmaH);
  fit->SetParLimits(1, meanL, meanH);
  hx->GetXaxis()->SetRangeUser(hx->GetBinCenter(meanbin)-1.5,hx->GetBinCenter(meanbin)+1.5);
  TF1* fitall  = new TF1("fitall","gaus(0)+gaus(3)+gaus(6)",-5, 15);
  hx->Fit(fit,"BR");
  fit->SetLineColor(kRed);
  hx->GetXaxis()->SetRangeUser(-10,35);
  fit->DrawCopy("same"); 
  mean = fit->GetParameter(1);
  sigma = fit->GetParameter(2); 
  TLatex lat;
  lat.SetTextSize(0.035);
  lat.DrawLatexNDC(0.2,0.85,Form("%.2f<p_{T}<%.2f GeV %s",lowpt,highpt,p.c_str()));
  // addpdf(pdf);
}
void funtofsamp(double* x,double *p)
{
  int idx = p[7]; //p[8] is no use now
  double e_part = p[0]*TMath::Gaus(x[0],p[1],p[2],1);
  // double p_part = p[3]*TMath::Gaus(x[0],p[4],p[5],1);
  double p_part = p[3]*hpsamp[idx]->GetBinContent(hpsamp[idx]->FindBin(x[0]));
  double pi_part=0;
  
  if (x[0]<10&&x[0]>-10)  
  {
    // if (hpisamp[idx]->GetBinContent(hpisamp[idx]->FindBin(x[0]))>0) pi_part= p[6]*hpisamp[idx]->GetBinContent(hpisamp[idx]->FindBin(x[0]));
    if (hpisamp[idx]->GetBinContent(hpisamp[idx]->FindBin(x[0]))>0) 
      pi_part= p[6]*gpisamp[idx]->Eval(x[0]);
  }
  // double K_part = p[9]*TMath::Gaus(x[0],p[10],p[11],1);
  double K_part = p[9]*hKsamp[idx]->GetBinContent(hKsamp[idx]->FindBin(x[0]));
  double Mg_part = p[12]*TMath::Gaus(x[0],p[13],p[14],1);
  return e_part+pi_part+K_part+p_part+Mg_part;
}
void funpisamp(double* x,double *p )
{
  int idx = p[1]; //p[8] is no use now
  double pi_part=0;
  if (x[0]<10&&x[0]>-10)  
  {
    // if (hpisamp[idx]->GetBinContent(hpisamp[idx]->FindBin(x[0]))>0) pi_part= p[6]*hpisamp[idx]->GetBinContent(hpisamp[idx]->FindBin(x[0]));
    if (hpisamp[idx]->GetBinContent(hpisamp[idx]->FindBin(x[0]))>0) 
      pi_part= p[0]*gpisamp[idx]->Eval(x[0]);
  }

  // if (x[0]<10&&x[0]>-10)  pi_part = p[0]*hpisamp[idx]->GetBinContent(hpisamp[idx]->FindBin(x[0]));
  return pi_part;
}
void funpsamp(double* x,double *p )
{
  int idx = p[1]; //p[8] is no use now
  double p_part=0;
  if (x[0]<10&&x[0]>-10)  p_part = p[0]*hpsamp[idx]->GetBinContent(hpsamp[idx]->FindBin(x[0]));
  return p_part;
}
void funKsamp(double* x,double *p )
{
  int idx = p[1]; //p[8] is no use now
  double K_part=0;
  if (x[0]<10&&x[0]>-10)  K_part = p[0]*hKsamp[idx]->GetBinContent(hKsamp[idx]->FindBin(x[0]));
  return K_part;
}
void fitNsigE(int centL=2, int centH=8){
  myStyle();
  TCanvas* c = new TCanvas("c","c",900,600);
  pdf = new TPDF(Form("NsigE_%d_%d.pdf",centL,centH));
  pdf->Off();
  drawtitle(pdf,c,"Photonic electron plots");
  TFile* f = TFile::Open("../0630/PE.root");
  TH2F* hels = (TH2F*)f->Get("hnSigE_e_ls"); 
  TH2F* he = (TH2F*)f->Get("hnSigE_e"); 
  he->Add(hels,-1);
  he->SetDirectory(0);
  TH2F* hp = (TH2F*)f->Get("hnSigE_p");
  hp->SetDirectory(0);
  TH2F* hk = (TH2F*)f->Get("hnSigE_k");
  hk->SetDirectory(0);
  TH2F* hpimg = (TH2F*)f->Get("hnSigE_piMg");
  hpimg->SetDirectory(0);

  TH2F* hpi = (TH2F*)f->Get("hnSigE_pi");
  hpi->SetDirectory(0);
  f->Close();

  f = TFile::Open("../PE.root");
  TH3F* htof = (TH3F*)f->Get("hnSigE_tof");
  htof->SetDirectory(0);
  f->Close();

  // f = TFile::Open("../0630/Ks.root");
  // f = TFile::Open("../pion/Ks.root");
  // TH3F* hpi3 = (TH3F*)f->Get("hnSigE_pi_tof");
  // TH3F* hpi3 = (TH3F*)f->Get("hnSigE_pi");
  // // hpi->RebinY(4);
  // TH3F* hpils3 = (TH3F*)f->Get("hnSigE_pi_ls");
  // // hpils->RebinY(3);
  // hpi3->Add(hpils3,-1);
  // TH2F* hpi = (TH2F*)hpi3->Project3D("yx");
  // hpi->SetDirectory(0);
  // f->Close();

  float pt[300]; //0.2-5
  int i=0;
  pt[0]=0.2;
  while (pt[i]<0.5){
    pt[i+1]=pt[i]+0.05;
    i++;
  }

  while (pt[i]<0.85){
    pt[i+1]=pt[i]+0.05;
    i++;
  }
  // cout<<pt[i]<<" "<<i<<endl;
  while  (pt[i]<1.5){
    pt[i+1]=pt[i]+0.05;
    i++;
  }
  cout<<pt[i]<<" "<<i<<endl;
  //
  while  (pt[i]<3.2){
    pt[i+1]=pt[i]+0.05;
    i++;
  }
  cout<<pt[i]<<" "<<i<<endl;

  int  bin =i;
  TH1F* hpurity = new TH1F("hpurity", "hpurity;electron p_{T}(GeV);purity",bin,pt);
  // TH1F* hpideff = new TH1F("hpideff", "hpideff;electron p_{T}(GeV);eff",bin,pt);
  TH1F* hkratio = new TH1F("hkratio", "hkratio;p(GeV);purity",bin,pt);
  TH1F* hpiratio = new TH1F("hpiratio", "hpiratio;p(GeV);purity",bin,pt);
  TH1F* hpratio = new TH1F("hpratio", "hpratio;p(GeV);purity",bin,pt);
  TH1F* hmgpiratio = new TH1F("hmgpiratio", "hmgpiratio;p(GeV);eff",bin,pt);
  TH1F* hpideff = new TH1F("hpideff", "hpideff;electron p(GeV);eff",bin,pt);
  TH1F* htemp = (TH1F*)hpi->ProjectionY("htemppi",1,10);
  htemp->SetLineColor(kGreen);
  htemp->SetLineStyle(7);

  c->Divide(3,2);
  int color[5]={kRed,kMagenta+1,kGreen+2,kBlue,kOrange+1};
  int line[5]={2,3,4,5,6};
  string pfitname[5]={"e","p","pi","k","pimg"};
  string plegname[5]={"e","proton","#pi","K","merged #pi"};
  TH1F* hyield[5];
  TH1F* hmean[5];
  TH1F* hsigma[5];
  for (int ip=0;ip<5;ip++){
    hyield[ip] = new TH1F(Form("hyield_%s",pfitname[ip].c_str()),Form("fitting const. %s;%s p_{T};Const.",plegname[ip].c_str(),plegname[ip].c_str()),bin,pt);
    hsigma[ip] = new TH1F(Form("hsigma_%s",pfitname[ip].c_str()),Form("sigma of nSigE for %s;%s p_{T};n#sigma e",plegname[ip].c_str(),plegname[ip].c_str()),bin,pt);
    hmean[ip] = new TH1F(Form("hmean_%s",pfitname[ip].c_str()),Form("mean of nSigE for  %s;%s p_{T};mean",plegname[ip].c_str(),plegname[ip].c_str()),bin,pt);
  }
  // TF1* ftot = new TF1("ftot","gausn(0)+gausn(3)+gausn(6)+gausn(9)+gausn(12)",-10,15);
  TF1* ftot = new TF1("ftot",funtofsamp,-10,15,15);
  TF1* fpion = new TF1("fpion",funpisamp,-10,10,2);
  TF1* fK = new TF1("fkaon",funKsamp,-10,10,2);
  TF1* fp = new TF1("fproton",funpsamp,-10,10,2);
  TF1* fpartical[5]; 
  for (int ip=0;ip<5;ip++){
    fpartical[ip] = new TF1(Form("f%s",pfitname[ip].c_str()),"gausn(0)",-10,15);
    fpartical[ip]->SetLineColor(color[ip]);
    fpartical[ip]->SetLineStyle(line[ip]);
  }
  TF1* fpimg = new TF1("fpimg","gausn(0)",-10,15);
  TLatex lax;
  TH1F* hpad = new TH1F("hpad","hpad",120,-10,15);
  TF1* frg_p = new TF1("frg_p","[0]+[1]/x/x+[2]/x+[3]*x",-10,35); 
  double frg_p_par[]={-17.8,-1.335,16.87,2.195};
  frg_p->SetParameters(frg_p_par);
  double frg_p_rg[2]={6,-7};
  TF1* frg_k = new TF1("frg_k","[0]+[1]/x/x+[2]/x+[3]*x",-10,35); 
  double frg_k_par[]={-12.5,0.289,5.469,1.946};
  frg_k->SetParameters(frg_k_par);
  double frg_k_rg[2]={6.5,-7.5};

  for (int j=0;j<bin;j++){
    cout<<pt[j]<<" "<<pt[j+1]<<endl;
    float mean[5]={0}, sigma[5]={0};
    c->cd(1);
    projectionAndFit(he,pt[j],pt[j+1],mean[0],sigma[0],-2,2,-0.6,0.1,0.9,1.1,"e",centL,centH);
    c->cd(2);
    projectionAndFit(hp,pt[j],pt[j+1],mean[1],sigma[1],frg_p->Eval(pt[j+1])-7,frg_p->Eval(pt[j])+6.5, -6, 30,0.5,2.5 ,"proton",centL,centH);
    c->cd(3);
    projectionAndFit(hpi,pt[j],pt[j+1],mean[2],sigma[2],-10,10, -6, 0, 0.2, 1.2,"#pi",centL,centH);
    c->cd(4);
    projectionAndFit(hk, pt[j],pt[j+1],mean[3],sigma[3],frg_k->Eval(pt[j+1])-6.5,frg_k->Eval(pt[j])+7.5,-6,15,0.5,3,"K",centL,centH);
    ipar++;
    c->cd(5);
    double meanpimg =  3;
    if (pt[j]>2) meanpimg = 4.8;
    projAndFitForMg(hpimg,pt[j],pt[j+1],mean[4],sigma[4],0,10,meanpimg,5.2,1.1,1.4,"merged #pi",centL,centH);
    c->cd(6);
    double par[15];
    int lbin = htof->GetXaxis()->FindBin(pt[j]);
    int hbin = htof->GetXaxis()->FindBin(pt[j+1]);
    int centbinL= htof->GetZaxis()->FindBin(centL);
    int centbinH= htof->GetZaxis()->FindBin(centH);

    TH1F* htotsamp = (TH1F*)htof->ProjectionY("htotx",lbin,hbin,centbinL,centbinH); 
    // htotsamp->Rebin(5);
    double constL[5]={ 100,1e2,1e2,1e2,1e1};
    // double constL[5]={ 100,1e4,1e4,1e2,1e2};
    double constH[5]={ 5e9,5e9,1e9,5e9,1e9};
    // double constH[5]={ 5e6,1e9,1e8,5e7,1e5};
    if (centH==8 && centL==6)
    {
    }
    if (centH==5 && centL==2)
    {
    }
    if (centH==8 && centL==0)
    {
    }
    if (centH==8 && centL==2)
    {
    }
    for (int i=0;i<5;i++){
      // if (i!=2) ftot->FixParameter(i*3+1, mean[i]);
      ftot->FixParameter(i*3+1, mean[i]);
      ftot->FixParameter(i*3+2, sigma[i]);
      fpartical[i]->FixParameter(1,mean[i]);
      fpartical[i]->FixParameter(2,sigma[i]);
      // ftot->SetParLimits(i*3, 0,htotsamp->Integral()/htotsamp->GetBinWidth(1));
      ftot->SetParLimits(i*3, constL[i], constH[i]);
      hmean[i]->SetBinContent(j+1, mean[i]);
      hsigma[i]->SetBinContent(j+2, sigma[i]);
    }
    ftot->FixParameter(7,j); //for pion
    //e
    htotsamp->GetXaxis()->SetRangeUser(-2,2); 
    htotsamp->Fit(fpartical[0],"RN");
    double e_const  = fpartical[0]->GetParameter(0);
    ftot->SetParameter(0,e_const);
    htotsamp->GetXaxis()->SetRangeUser(-10,20); 
    //p 
    // htotsamp->GetXaxis()->SetRangeUser( frg_p->Eval(pt[j+1])-7,frg_p->Eval(pt[j])+6.5); 
    // htotsamp->Fit(fpartical[1],"R");
    // double pi_const  = fpartical[1]->GetParameter(0);
    // ftot->SetParameter(3,pi_const);
    // htotsamp->GetXaxis()->SetRangeUser(-10,20); 
    //pi
    htotsamp->GetXaxis()->SetRangeUser(-10,-2); 
    htotsamp->Fit(fpartical[2],"RN");
    double pi_const  = fpartical[2]->GetParameter(0);
    ftot->SetParameter(6,pi_const);
    htotsamp->GetXaxis()->SetRangeUser(-10,20); 
    // //K
    // htotsamp->GetXaxis()->SetRangeUser( frg_k->Eval(pt[j+1])-6.5,frg_k->Eval(pt[j])+7.5); 
    // htotsamp->Fit(fpartical[3],"R");
    // double k_const  = fpartical[3]->GetParameter(0);
    // ftot->SetParameter(9,k_const);
    // htotsamp->GetXaxis()->SetRangeUser(-10,20); 
    // //merged  pi
    // htotsamp->GetXaxis()->SetRangeUser(3,6); 
    // htotsamp->Fit(fpartical[4],"R");
    // double pimg_const = fpartical[4]->GetParameter(0);
    // ftot->SetParameter(12,pimg_const);
    // htotsamp->GetXaxis()->SetRangeUser(-10,20); 
    for (int ip=0;ip<5;ip++){
      ftot->SetParLimits(i*3, constL[ip], constH[ip]); 
      // ftot->SetParameter(i*3, 0.5*(constL[ip]+constH[ip])); 
    }
    gPad->SetLogy();
    htotsamp->Draw();
    htotsamp->GetXaxis()->SetRangeUser(-10,5);
    htotsamp->Fit(ftot,"RB");
    htotsamp->GetXaxis()->SetRangeUser(-10,15);
    htotsamp->Fit(ftot,"RB");
    ftot->GetParameters(par);

    TLegend* leg = new TLegend(0.7,0.7,0.88,0.88);
    for (int ip=0;ip<5;ip++){
      fpartical[ip]->SetParameters(par+ip*3);  
      if (ip!=2)  fpartical[ip]->Draw("same");
      else if (ip==2)
      {
        htemp->Rebin(hpisamp[ipar-1]->GetBinWidth(1)/htemp->GetBinWidth(1));
        htemp->Add(htemp,hpisamp[ipar-1],0,1); 
        htemp->Scale(par[ip*3]);
        htemp->Draw("same");
        fpion->SetParameter(0,par[ip*3]);
        fpion->SetParameter(1,ipar-1);
      }
      // if (ip!=2) leg->AddEntry(fpartical[ip],plegname[ip].c_str(),"l");
      leg->AddEntry(fpartical[ip],plegname[ip].c_str(),"l");
      hyield[ip]->SetBinContent(j+1, par[ip*3] );
      hyield[ip]->SetBinError(j+1, ftot->GetParError(ip*3) );
    }

    leg->Draw();
    c->cd();
    addpdf(pdf);
    double meanpt = 0.5*(pt[j]+pt[j+1]);
    // double lowsigma = (meanpt-0.2)*1.5625-2.34;
    // double lowsigma = meanpt*2-2.8;
    double lowsigma = 0;
    // double lowsigma = -1;
    double lowsigma = meanpt*3.5-2.8;
    if(meanpt>0.8) lowsigma = 0;
    // float total = htotsamp->Integral(htotsamp->FindBin(lowsigma),htotsamp->FindBin(2));
    float total = ftot->Integral(lowsigma,2)/htotsamp->GetBinWidth(1);
    // float total = ftot->Integral(lowsigma,3)/htotsamp->GetBinWidth(1);
    float electron = fpartical[0]->Integral(lowsigma,2)/htotsamp->GetBinWidth(1);
    // float electron = fpartical[0]->Integral(lowsigma,1.66)/htotsamp->GetBinWidth(1);
    float proton = fpartical[1]->Integral(lowsigma,2)/htotsamp->GetBinWidth(1);
    // float pi = fpartical[2]->Integral(lowsigma,2)/htotsamp->GetBinWidth(1);
    float pi = fpion->Integral(lowsigma,2)/htotsamp->GetBinWidth(1);
    float kaon = fpartical[3]->Integral(lowsigma,2)/htotsamp->GetBinWidth(1);
    float mgpion = fpartical[4]->Integral(lowsigma,2)/htotsamp->GetBinWidth(1);

    float purity = electron*1.0/total;
    float eff = fpartical[0]->Integral(-3,3)/fpartical[0]->Integral(-6,6);
    hpideff->SetBinContent(j+1,eff);
    hpurity->SetBinContent(j+1,purity);
    hpratio->SetBinContent(j+1,proton*1./total);
    hpiratio->SetBinContent(j+1,pi*1./total);
    hkratio->SetBinContent(j+1,kaon*1./total);
    hmgpiratio->SetBinContent(j+1,mgpion*1./total);
  }
  gPad->SetLogy(0);
  hpurity->GetYaxis()->SetRangeUser(0.2,1.1);
  hpurity->Draw();
  hpurity->SetDirectory(0);
  addpdf(pdf);
  hpideff->Draw();
  hpideff->SetDirectory(0);
  addpdf(pdf);
  hpurity->SetLineColor(color[0]);
  hpratio->SetLineColor(color[1]);
  hpiratio->SetLineColor(color[2]);
  hkratio->SetLineColor(color[3]);
  hmgpiratio->SetLineColor(color[4]);
  hpurity->Draw();
  hpurity->GetYaxis()->SetRangeUser(0,1.1);
  hpratio->Draw("same");
  hpiratio->Draw("same");
  hkratio->Draw("same");
  hmgpiratio->Draw("same");
  addpdf(pdf);

  TCanvas* cp = new TCanvas("cp","cp");
  cp->Divide(2,2);
  cp->cd();
  TLegend* lp = new TLegend(0.2,0.2,0.8,0.8);
  for (int ip=0;ip<5;ip++){
    cp->cd();
    cp->cd(1);
    gPad->SetLogy();
    hyield[ip]->SetLineColor(color[ip]);
    hyield[ip]->SetMarkerColor(color[ip]);
    hyield[ip]->SetMarkerSize(0.6);
    hyield[ip]->GetYaxis()->SetRangeUser(10,1e9);
    hyield[ip]->Draw("same");
    hyield[ip]->SetDirectory(0);
    // gPad->Update();
    cp->cd(2);
    hsigma[ip]->SetLineColor(color[ip]);
    hsigma[ip]->GetYaxis()->SetRangeUser(0.5,2.5);
    hsigma[ip]->Draw("same");
    hsigma[ip]->SetDirectory(0);
    // gPad->Update();
    cp->cd(3);
    hmean[ip]->SetLineColor(color[ip]);
    hmean[ip]->GetYaxis()->SetRangeUser(-10,15);
    hmean[ip]->Draw("same");
    hmean[ip]->SetDirectory(0);
    lp->AddEntry(hyield[ip],plegname[ip].c_str(),"l");
    // gPad->Update();
  }

  cp->cd(4);
  gPad->Draw();
  lp->Draw();
  cp->cd();
  pdf->On();
  pdf->NewPage();
  cp->Update();
  pdf->Off(); 

  pdf->On();
  pdf->Close();
  TFile* fNsig = new TFile(Form("Nsigma_%d_%d.root",centL,centH),"recreate");
  hpurity->Write();
  hpiratio->Write();
  hkratio->Write();
  hpratio->Write();
  hmgpiratio->Write();

  hpideff->Write();
  for (int ip=0;ip<5;ip++){
    hyield[ip]->Write();
    hsigma[ip]->Write();
    hmean[ip]->Write();
  }
}
