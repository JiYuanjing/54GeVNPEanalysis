#include "sPhenixStyle.h"
void myFit(double *x, double *p)
{
  if ( x[0]<p[6] ) return p[0]+x[0]*p[1]+x[0]*x[0]*p[2]+x[0]*x[0]*x[0]*p[3]+pow(x[0],4)*p[4];
  if (x[0]>=p[6]) return  p[0]+p[6]*p[1]+p[6]*p[6]*p[2]+pow(p[6],3)*p[3]+pow(p[6],4)*p[4]-p[6]*p[5]+x[0]*p[5];
}
void fitfun(double *x, double * p)
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
void fitv2()
{
  SetsPhenixStyle();
  Double_t Ks_pt_bin_center[16] = {0.3055,0.5055,0.7055,0.9055,1.1055,1.2955,1.5055,1.7055,1.9045,2.1055,2.3055,2.5055,2.7055,3.0055,3.4055,3.9755};
  Double_t Ks_v2_values[16] = {0.0126035,0.0283912,0.0477735,0.0655328,0.0810368,0.0942137,0.106462,0.115122,0.124632,0.130506,0.132923,0.137253,0.137584,0.13788,0.133286,0.119046};
  Double_t Ks_v2_stat_error[16] = {0.000905458,0.000526192,0.000436614,0.000453511,0.000533842,0.000670305,0.000874881,0.00116362,0.00155222,0.00206399,0.00273512,0.00361251,0.00475974,0.00494931,0.00833496,0.011402};
  Double_t Ks_v2_syst_low_error[16] = {6.10066e-05,9.64478e-06,8.13875e-06,6.17024e-05,1.83154e-05,0.000186825,0.000135778,0.000112175,0.000177097,0.000254422,0.00118586,8.49751e-05,0.000872825,0.00119956,0.00157152,0.00590387};
  Double_t Ks_v2_syst_high_error[16] = {5.94386e-05,9.78101e-06,8.23048e-06,6.05527e-05,1.78633e-05,0.000181043,0.000141012,0.000107432,0.000186081,0.000268904,0.00126215,7.93461e-05,0.000804378,0.00120266,0.0015422,0.00514164};

  Double_t Kp_pt_bin_center[15] = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
  Double_t Kp_v2_values[15] = {0.0127996,0.0281504,0.0479773,0.0667651,0.0833767,0.0954043,0.106832,0.116863,0.122605,0.128074,0.132259,0.137232,0.13796,0.137573,0.145755};
  Double_t Kp_v2_stat_error[15] = {0.00033299,0.000204027,0.000201565,0.000232314,0.000285152,0.000365635,0.000481724,0.000640665,0.000868916,0.001205,0.00168448,0.0023502,0.00332723,0.00374709,0.00529857};
  Double_t Kp_v2_syst_low_error[15] = {2.29511e-05,2.53346e-05,3.87823e-05,2.24753e-05,1.7626e-05,1.76556e-05,2.49181e-05,0.000108266,9.6143e-05,0.000189906,0.000497042,0.000317001,0.000361833,0.00103781,0.0017297};
  Double_t Kp_v2_syst_high_error[15] = {3.82051e-05,3.81784e-05,5.82484e-05,5.40282e-05,6.28527e-05,0.000160831,0.000181208,0.000117551,0.000119296,0.000291899,0.0005016,0.000373801,0.000509016,0.000919941,0.0022869};

  Double_t Km_pt_bin_center[15] = {0.315,0.495,0.705,0.885,1.095,1.275,1.485,1.695,1.875,2.085,2.295,2.475,2.685,2.955,3.345};
  Double_t Km_v2_values[15] = {0.0119524,0.0275439,0.0478142,0.0668999,0.0828742,0.0968429,0.107873,0.116527,0.125162,0.127785,0.132031,0.131966,0.133757,0.136496,0.151896};
  Double_t Km_v2_stat_error[15] = {0.000364338,0.00022339,0.000221358,0.000256186,0.000316043,0.000406005,0.000540067,0.000728892,0.00100546,0.00142133,0.00182739,0.00284911,0.00384224,0.00460709,0.00645661};
  Double_t Km_v2_syst_low_error[15] = {2.36593e-05,3.69839e-05,6.01284e-05,1.78588e-05,1.52751e-05,1.85946e-05,9.33403e-05,8.34216e-05,9.25011e-05,0.000456588,0.000584889,0.000917781,0.00170289,0.00223047,0.000581813};
  Double_t Km_v2_syst_high_error[15] = {3.74267e-05,5.26462e-05,8.28083e-05,5.74421e-05,0.000150701,0.000177623,9.62286e-05,0.000118425,0.000135989,0.000459601,0.000698035,0.000912625,0.00180228,0.00239182,0.00071318};

  TGraphErrors* gKm = new TGraphErrors(15,Km_pt_bin_center, Km_v2_values, 0, Km_v2_stat_error );
  TGraphErrors* gKp = new TGraphErrors(15,Kp_pt_bin_center, Kp_v2_values, 0, Kp_v2_stat_error );
  TGraphErrors* gKs = new TGraphErrors(16,Ks_pt_bin_center, Ks_v2_values, 0, Ks_v2_stat_error );

  TMultiGraph* gcom=new TMultiGraph("gcom","gcom");
  TF1* fit = new TF1("fit", myFit, 0, 5, 7);
  fit->FixParameter(6,2.5);
  fit->SetParameters(-0.01,0.05,0.07,-0.05,0.00965,0.0027,2.5);
  /* TF1* fit = new TF1("fit", fitfun, 0, 5, 3); */
  /* fit->SetParameters(0.1,-0.1,1); */
  gKm->SetMarkerStyle(kOpenCircle);
  gKs->SetMarkerStyle(kOpenCircle);
  gKp->SetMarkerStyle(kOpenCircle);
  gKm->SetMarkerColor(1);
  gKs->SetMarkerColor(2);
  gKp->SetMarkerColor(3);

  gcom->Add(gKs);
  gcom->Add(gKm);
  gcom->Add(gKp);
  gKm->SetMarkerSize(1.5);
  gKs->SetMarkerSize(1.5);
  gKp->SetMarkerSize(1.5);
 
  gcom->Draw("AP");
  gcom->Fit(fit);
  fit->SaveAs("KaonV2.root");
}
