
#include "sPhenixStyle.h"
Double_t  myFit(double *x, double *p)
{
  if ( x[0]<p[6] ) return p[0]+x[0]*p[1]+x[0]*x[0]*p[2]+x[0]*x[0]*x[0]*p[3]+pow(x[0],4)*p[4];
  else if (x[0]>=p[6]) return  p[0]+p[6]*p[1]+p[6]*p[6]*p[2]+pow(p[6],3)*p[3]+pow(p[6],4)*p[4]-p[6]*p[5]+x[0]*p[5];
  else return 0;
}
Double_t fitfun(double *x, double * p)
{
  double xv = x[0];
  double v2;
  if (p[2]!=0)
  {
    v2 = p[0]/(1.+TMath::Exp(-(xv-p[1])/p[2]))-p[0]/(1.+TMath::Exp((p[1])/p[2]));
  }
  else v2 = 0;
  return v2;

}
void fit54v2()
{
  SetsPhenixStyle();
  double meanpt_kaplus_0060_54gev[14] ={0.297735, 0.500618, 0.69614, 0.894909, 1.09392, 1.29262, 1.49072, 1.68779, 1.88403, 2.09178, 2.31548, 2.51433, 2.71073, 2.90832};
  double meanpt_kaminus_0060_54gev[14]={0.296627, 0.500496, 0.696005, 0.894732, 1.09367, 1.29227, 1.49017, 1.68695, 1.88355, 2.09818, 2.31643, 2.51349, 2.71013, 2.90793};
 double val_kaplus_0060_54gev[14] ={0.0125859,0.0280836,0.0474321,0.0652961,0.0802486,0.0928426,0.102274,0.108638,0.1152,0.12043,0.126605,0.132931,0.140953,0.143233,};
 double err_kaplus_0060_54gev[14] ={9.20666e-05,5.67282e-05,5.62234e-05,6.57851e-05,8.36465e-05,0.00010862,0.000147658,0.000211829,0.00030532,0.000456551,0.000677507,0.000958456,0.0012659,0.001623,};
 double val_kaminus_0060_54gev[14] ={0.012131,0.0273535,0.0465513,0.0645764,0.0798733,0.0925996,0.10249,0.109127,0.114397,0.12037,0.125258,0.130451,0.133229,0.138455,};
 double err_kaminus_0060_54gev[14] ={0.00010137,6.27377e-05,6.24903e-05,7.37491e-05,9.43174e-05,0.000123005,0.00016962,0.000250868,0.000378593,0.000557424,0.000825085,0.00118238,0.00161504,0.00215492,};
  TGraphErrors* gKm = new TGraphErrors(14, meanpt_kaminus_0060_54gev, val_kaminus_0060_54gev, 0, err_kaminus_0060_54gev);
  TGraphErrors* gKp = new TGraphErrors(14, meanpt_kaplus_0060_54gev, val_kaplus_0060_54gev, 0, err_kaplus_0060_54gev);

  TMultiGraph* gcom=new TMultiGraph("gcom","gcom");
  TF1* fit = new TF1("fit", myFit, 0, 5, 7);
  /* TF1* fit = new TF1("fit", fitfun, 0, 3, 3); */
  /* fit->FixParameter(6,2.5); */
  fit->SetParameters(-0.01,0.05,0.07,-0.05,0.0,0.0027,3);
  fit->FixParameter(6,2.9);
  /* TF1* fit = new TF1("fit", fitfun, 0, 5, 3); */
  /* fit->SetParameters(0.1,-0.1,1); */
  gKm->SetMarkerStyle(kOpenCircle);
  /* gKs->SetMarkerStyle(kOpenCircle); */
  gKp->SetMarkerStyle(kOpenCircle);
  gKm->SetMarkerColor(kBlue);
  /* gKs->SetMarkerColor(2); */
  gKp->SetMarkerColor(kGreen+2);

  /* gcom->Add(gKs); */
  gcom->Add(gKm);
  gcom->Add(gKp);
  gKm->SetMarkerSize(1.5);
  /* gKs->SetMarkerSize(1.5); */
  gKp->SetMarkerSize(1.5);
  gcom->Draw("AP");
  gcom->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gcom->GetYaxis()->SetTitle("Kaon v_{2}");
  /* gcom->Fit(fit); */
  TGraphErrors* gnew = (TGraphErrors*)gKp->Clone("gCombine");
  for (int i=0;i<gnew->GetN();i++)
  {
      gnew->SetPointY( i, (gKp->GetPointY(i)+gKm->GetPointY(i))*0.5);
  }
  gnew->SetMarkerStyle(kFullCircle);
  gnew->SetMarkerColor(kRed);
  TLegend* leg = new TLegend(0.7,0.2,0.88,0.48);
  leg->AddEntry(gnew,"Average","pe");
  leg->AddEntry(gKm,"K^{-}","pe");
  leg->AddEntry(gKp,"K^{+}","ep");
  leg->Draw();
  gnew->Draw("p c same");
  gPad->SaveAs("Kaon54V2.pdf");
  /* fit->SaveAs("Kaon54V2.root"); */
  TFile* fnew = new TFile("Kaon54V2.root","recreate");
  gnew->Write();
  gcom->Write();
  fnew->Close();

}
