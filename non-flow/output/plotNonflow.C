#include "myStyle.h"
void plotNonflow()
{
  myStyle();
  TFile* full = new TFile("nonflow_fullrange.root");
  TH1F* hfull = (TH1F*)full->Get("hM2_EtaGap_dot0")->Clone("hnonflow_away");
  hfull->SetDirectory(0);
  TFile* half = new TFile("nonflow_halfrange.root");
  TH1F* hhalf = (TH1F*)half->Get("hM2_EtaGap_dot0")->Clone("hnonflow_near");
  hhalf->SetDirectory(0);

  hfull->SetMarkerStyle(20);
  hhalf->SetMarkerStyle(25);
  hfull->Draw();
  hfull->GetYaxis()->SetTitleOffset(1.2);
  hhalf->Draw("same");
  int const  np = hfull->GetNbinsX();
  double x[100],y[100],err[100];
  double xup[100],yup[100], lowerr[100];
  for (int ip=0;ip<np;ip++)
  {
    x[ip]=hfull->GetBinCenter(ip+1);
    y[ip]=0.5*(hfull->GetBinContent(ip+1)+hhalf->GetBinContent(ip+1));
    err[ip]=(hfull->GetBinError(ip+1)*hfull->GetBinError(ip+1)+hhalf->GetBinError(ip+1)*hhalf->GetBinError(ip+1))*0.5*0.5;
    err[ip]=err[ip]+pow(hfull->GetBinContent(ip+1)-y[ip],2);
    err[ip]=sqrt(err[ip]);
    yup[ip]=hfull->GetBinContent(ip+1);
    lowerr[ip]=yup[ip];
  }

  TGraphErrors* gnonflow = new TGraphErrors(np,x,y,0,err);
  gnonflow->SetFillStyle(3004);
  gnonflow->SetFillColor(kRed);
  gnonflow->Draw("same3"); 
  gnonflow->SetName("nonflow");

  TGraph* gnonflowfull = new TGraph(np-1,x+1,yup+1);
  TGraphAsymmErrors* gnonflowfullas = new TGraphAsymmErrors(np-1,x+1,yup+1,0,0,lowerr+1,0);
  gnonflowfullas->SetName("nonflowfullAs");
  gnonflowfull->SetFillStyle(3004);
  gnonflowfull->SetFillColor(kRed);
  /* gnonflowfull->Draw("same3");  */
  gnonflowfull->SetName("nonflowfull");

  /* gnonflow->Draw("same l");  */
  //nonflow for 200GeV
  double nonflow200_pt[6] = {1.11351,1.35135,1.60541,1.85946,2.47027,2.82162};
  double nonflow200_val[6] = {0.0277577,0.0369801,0.0396926,0.0497288,0.0556962,0.0660036 };
  double nonflow200_errL[6] = {0.00135624,0.00271248,0.00379747,0.00786618,0.00352622,0.00678119};
  double nonflow200_errH[6] = {0.00135624,0.00271248,0.00379747,0.00786618,0.00352622,0.00678119};
  TGraphAsymmErrors* g200nonflow = new TGraphAsymmErrors( 6, nonflow200_pt, nonflow200_val,0,0,nonflow200_errL,nonflow200_errH);
  g200nonflow->SetName("g200nonflow");
  g200nonflow->SetFillColor(kGray);
  g200nonflow->SetFillStyle(1001);
  g200nonflow->Draw("same3");
  TLegend* leg = new TLegend(0.7,0.7,0.88,0.88);
  leg->SetTextSize(0.055);
  leg->SetTextFont(42);
  leg->AddEntry(hfull,"full #phi range","p");
  leg->AddEntry(hhalf,"near side","p");
  leg->AddEntry(g200nonflow,"200 GeV","f");
  leg->Draw();

  TLatex tex;
  tex.DrawLatexNDC(0.2,0.68, "v_{2}^{AA}=0.046"); 
  tex.DrawLatexNDC(0.2,0.78, "<N_{h}^{AA}>=122");
  tex.DrawLatexNDC(0.2,0.88,"0.2<p_{T}^{h}<2 GeV"); 
  TFile* f = new TFile("nonflow54.root","recreate");
  g200nonflow->Write();
  gnonflowfull->Write();
  gnonflowfullas->Write();

}
