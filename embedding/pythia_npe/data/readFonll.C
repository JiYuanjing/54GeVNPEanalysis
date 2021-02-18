
TGraph* g;
double myfun(double *x, double *par)
{
  if (x[0]<10) return g->Eval(x[0])*par[0];
  else return 0;
}
void readFonll()
{
  ifstream read;
  read.open("fonll62D0.txt");
  char line[1000];
  /* getline( read, line,'\n' ); */
  /* cout<< line<<endl; */
  /* getline( read, line,'\n' ); */
  /* getline( read, line,'\n' ); */
  /* cout<< line<<endl; */

  double pt[500], yield[500];
  for (int i=0;i<500;i++)
  {
    read>>pt[i]>>yield[i];
    cout<<pt[i]<<endl;
  }
  
  g = new TGraph(500,pt,yield);
  TH1F* h = new TH1F("h","h",1,0,10);
  /* h->Fill(1); */
  h->Draw("c");
  h->GetYaxis()->SetRangeUser(1e-5,1e9.);
  gPad->SetLogy();
  /* TF1* fpispectra = new TF1("fpi0spectra","[0]*([2]-1)*([2]-2)/([2]*[1]*([2]*[1]+1.864500*([2]-2)))*TMath::Power(([2]*[1]+TMath::Sqrt(x[0]*x[0]+1.864500*1.864500)-1.864500)/([2]*[1]),-[2])*x[0]",0,10); */
  /* TF1* fpispectra = new TF1("fpi0spectra","[0]*([2]-1)*([2]-2)/([2]*[1]*([2]*[1]+1.864500*([2]-2)))*TMath::Power(([2]*[1]+TMath::Sqrt(x[0]*x[0]+1.864500*1.864500)-1.864500)/([2]*[1]),-[2])*x[0]",0,10); */
  TF1* fpispectra = new TF1("fpi0spectra","[0]*([2]-1)*([2]-2)/(([2]-3)*([2]-3))/[1]/[1]*TMath::Power(1+2*x[0]/[1]/([2]-3)/2,-1*[2])*x[0]",1,10);
  g->Draw("psame");
  /* for () { */
    /* fpi0spectra->SetParameters(3e7,1.3,30); */
    /* g->Fit(fpispectra); */
  /* } */
  TF1* fun = new TF1( "fonll62D0" ,myfun,0,10 , 1 );
  fun->SetParameter(0,1);
  fun->SetNpx(2000);
  fun->Draw("same");
  /* fun->SaveAs("test.root"); */
  TFile* f = new TFile("fonll62.root","recreate");
  double Ncoll[9]={ 10,25, 55, 105,186.45,310.8 ,494 ,691.6,861};
  TF1* fit[9];
  for (int ic=0;ic<9;ic++)
  {
    fit[ic] = new TF1(Form("FONLL62D0_Ncoll%d",ic), myfun,0,10 , 1 );
    fit[ic]->SetParameter(0,Ncoll[ic]);

    fit[ic]->SetNpx(2000);
    fit[ic]->Write();
  }

}
