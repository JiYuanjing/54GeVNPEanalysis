#include "sPhenixStyle.h"
void drawrecenter()
{
  SetsPhenixStyle();
  // TFile* file = TFile::Open("getshift.root");
  TFile* file = TFile::Open("incEv2_0.root");
  TH2F* hEventPlaneCent2 = (TH2F*)file->Get("hEventPlaneCent_M_Sh"); 
  TH2F* hEventPlaneCent1 = (TH2F*)file->Get("hEventPlaneCent_M_Re"); 
  TH2F* hEventPlaneCent0 = (TH2F*)file->Get("hEventPlaneCent_M"); 
  // TH2F* hEventPlaneCent2 = (TH2F*)file->Get("hEventPlaneCent_P_Sh"); 
  // TH2F* hEventPlaneCent1 = (TH2F*)file->Get("hEventPlaneCent_P_Re"); 
  // TH2F* hEventPlaneCent0 = (TH2F*)file->Get("hEventPlaneCent_P"); 

  TH1F* mbep1 = (TH1F*)hEventPlaneCent1->ProjectionY("hh1",3,9);
  TH1F* mbep2 = (TH1F*)hEventPlaneCent2->ProjectionY("hh2",3,9);
  TH1F* mbep0 = (TH1F*)hEventPlaneCent0->ProjectionY("hh0",3,9);
  mbep2->SetLineColor(kRed);
  mbep0->SetLineColor(kBlue);
  double  mean = mbep2->GetEntries()*1.0/(1.0*mbep2->GetNbinsX());
  TLine* line = new TLine(0,mean,3.14,mean);
  line->SetLineStyle(2);
  mbep0->Draw();
  mbep0->GetXaxis()->SetTitle("#Phi");
  mbep0->GetYaxis()->SetTitle("Counts");
  mbep2->Draw("same");
  mbep1->Draw("same");
  line->Draw("same");
  TLegend*  l = new TLegend(0.55,0.7,0.88,0.88);
  l->SetTextSize(0.055);
  l->AddEntry(mbep0,"raw","l" );
  l->AddEntry(mbep1,"after recenter","l" );
  l->AddEntry(mbep2,"after shift","l" );
  l->Draw();
  TLatex lat;
  lat.DrawLatexNDC(0.2,0.28,"Au+Au 54.4 GeV"); 
  lat.DrawLatexNDC(0.25,0.23,"0-60%  #eta<0"); 
}
