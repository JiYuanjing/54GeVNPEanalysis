#include "sPhenixStyle.h"
double funK0s(double *x, double *par)
{
	float m_pi0 = 0.135;
	float m_kaon = 0.493677;

	float pt = sqrt(x[0]*x[0] + m_kaon*m_kaon - m_pi0*m_pi0);
  return par[0]*(pow(exp(-par[1]*pt - par[2]*pt*pt) + pt/par[3], -par[4]));
}
void fit()
{
  TCanvas* c = new TCanvas("c","c");
  TPDF* pdf = new TPDF("test.pdf");
  TGraphAsymmErrors* g[7];
  int name[8]={0,5,10,20,30,40,60,80};
  int color[9]={kBlack,kRed,kViolet-1,kBlue,kGreen+2,kGray+2,kAzure-1,kPink-1,kOrange+9};
  for (int i=1;i<8;i++)
  {
    g[i-1] = new TGraphAsymmErrors(Form("%d_%d.txt",name[i-1],name[i]), "%lg %lg %lg %lg");
    g[i-1]->SetName(Form("%d_%d",name[i-1],name[i]));
    g[i-1]->SetMarkerColor(color[i-1]);
    g[i-1]->SetMarkerSize(1.5);
    g[i-1]->SetMarkerStyle(kOpenCircle);
    g[i-1]->SetLineColor(color[i-1]);
    cout << Form("%d_%d.txt",name[i-1],name[i]) <<endl;
  float Npart[9]={ 346.5,293.9,229.8,164.1,114.3,76.3,47.9,27.8,15.3};
  }

  TGraphAsymmErrors* gcomb = (TGraphAsymmErrors*)g[1]->Clone("0_60");
  double ratio[7] = {0.5,0.5,1,1,1,2,2};
  float Npart[9]={ 346.5,293.9,229.8,164.1,114.3,76.3,47.9,27.8,15.3};
  float Npart2[7]={ 346.5,293.9,229.8,164.1,114.3,(76.3+47.9)*0.5,(27.8+15.3)*0.5};
  double sum=0;
  for (int i=0;i<7;i++) sum+=ratio[i]; 
  for (int ip=0;ip<gcomb->GetN();ip++)
  {
    double tmp = 0;
    double tmperrh = 0,tmperrl =0;
    for (int i=0;i<6;i++){
      tmperrh+=pow(g[i]->GetErrorYhigh(ip),2);
      tmperrl+=pow(g[i]->GetErrorYlow(ip),2);
      tmp+=g[i]->GetPointY(ip);
    }
    cout<<tmperrh <<endl;
    tmperrl = sqrt(tmperrl)/sum;
    tmperrh = sqrt(tmperrh)/sum;
    tmp/=sum;
    gcomb->SetPointEYhigh(ip, tmperrh);
    gcomb->SetPointEYlow(ip, tmperrl);
    gcomb->SetPointY(ip, tmp);
  } 

  TF1* fun[8];
  double para1[7][5];
  double para2[7][5];
  /* float fitParams62GeV[5] = {6923.29,-0.169028,0.308199,0.82446,11.5055}; */
  double fitParams62GeV[5] = {6.71676e+01, 2.94912e-01,5.24846e-02,1.61808e+00,1.41987e+01};
  fun[7] = new TF1("fun_0_60", funK0s, 0,15,5);
  for (int ip=0;ip<5;ip++) fun[7]->SetParameter(ip,fitParams62GeV[ip]);
  
  gcomb->Draw("AP");
  gcomb->SetMarkerColor(kBlue);
  gcomb->Fit(fun[7]);
  gPad->SetLogy();
  fun[7]->GetParameters(fitParams62GeV);
  fun[7]->Draw("same");
  addpdf(pdf);

  for (int i=0; i<7; i++)
  {
    g[i]->Draw("pA");
    /* fun[i] = new TF1(Form("fun_%d_%d",name[i-1],name[i]),"[0]*2.*3.1415*pow(TMath::Exp(-1*[1]*x-[2]*x*x)+x/[3], -1*[4])",0,15); */
    fun[i] = new TF1(Form("fun_%d_%d",name[i],name[i+1]), funK0s, 0,5,5);
    for (int ip=0;ip<5;ip++)fun[i]->SetParameter(ip,fitParams62GeV[ip]);
    fun[i]->SetParameter(0, Npart2[i]/158.7*fitParams62GeV[0]);
     
    /* for (int ip=0;ip<5;ip++) fun[i]->SetParameter(ip,para[i][ip]); */
    g[i]->Fit(fun[i],"R+","",1.2,5);
    fun[i]->GetParameters(para2[i]);

    /* g[i]->Fit(fun[i],"R+","",0,1.75); */
    /* fun[i]->GetParameters(para1[i]); */
    gPad->SetLogy();
    drawLatex(0.15,0.15,Form("%d_%d",name[i],name[i+1]));
    addpdf(pdf);
  }

  for (int i=0; i<7; i++)
  {
    g[i]->Draw("pA");
    /* fun[i] = new TF1(Form("fun_%d_%d",name[i-1],name[i]), funK0s, 0,5,5); */
    /* for (int ip=0;ip<5;ip++)fun[i]->SetParameter(ip,fitParams62GeV[ip]); */
    /* fun[i]->SetParameter(0, Npart[i]/158.7*fitParams62GeV[0]); */
     
    /* g[i]->Fit(fun[i],"R+","",1.5,5); */
    /* fun[i]->GetParameters(para2[i]); */

    g[i]->Fit(fun[i],"R+","",0,1.75);
    fun[i]->GetParameters(para1[i]);
    gPad->SetLogy();
    drawLatex(0.15,0.15,Form("%d_%d",name[i],name[i+1]));
    addpdf(pdf);
  }

  cout<< "0_60 ";
  for (int ip=0;ip<5;ip++) cout << fitParams62GeV[ip]<<" ";
  cout << endl;
  cout <<"0-1.75: " <<endl;
  for (int i=6;i>-1;i--)
  {
    cout << Form("%d_%d",name[i],name[i+1])<<" : ";
    for (int ip=0;ip<5;ip++)
      cout << para1[i][ip] <<", "; 
    cout<<endl;
  }
  cout <<"1-5GeV: " <<endl;
  for (int i=6;i>-1;i--)
  {
    cout << Form("%d_%d",name[i],name[i+1])<<" : {";
    for (int ip=0;ip<5;ip++)
      cout << para2[i][ip] <<", "; 
    cout<<"}"<<endl;
  }

  pdf->On();
  pdf->Close();

  /* TH1F* h = new TH1F("h","h",200,0,5); */
  /* for (int i=0;i<200;i++) */
  /* { */
  /*    h->SetBinContent(i+1, g[0]->Eval(h->GetBinCenter(i+1))); */
  /* } */
  /* h->Draw(); */
  /* g[0]->Draw("same"); */

  TFile*fnew = new TFile("Ksnew.root","recreate");
  for (int i=0;i<7;i++){
    g[i]->Write();
  }
  gcomb->Write();

}
