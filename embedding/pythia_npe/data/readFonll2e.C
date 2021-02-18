TGraph* g;
double myfun(double *x, double *par)
{
  if (x[0]<10) return g->Eval(x[0])*par[0];
  else return 0;
}
void readFonll2e(int mode=2)
{
  TString name="";
  if (mode==1) name.Append("D");
  else (mode==2) name.Append("B");
  ifstream read;
  read.open(Form("fonll%s2e.txt",name.Data()));
  char line[1000];
  /* getline( read, line,'\n' ); */
  /* cout<< line<<endl; */
  /* getline( read, line,'\n' ); */
  /* getline( read, line,'\n' ); */
  /* cout<< line<<endl; */

  TH1F* h = new TH1F(Form("h%s2e",name.Data()),Form("h%s2e",name.Data()),400,0,4);
  double pt[500], yield[500], error[500];
  int remove=0;
  for (int i=0;i<400;i++)
  {
    read>>pt[i]>>yield[i]>>error[i];
    /* yield[i]=yield[i]*1e-6; */
    /* error[i]=error[i]*1e-6; */
    cout<<pt[i]<< " "<<error[i]<<endl;
    if (yield[i]<0) remove++;
    continue;
    h->SetBinContent(i+1, yield[i]);
    h->SetBinError(i+1, error[i]);
  }
 
  TGraphErrors* gerr = new TGraphErrors(400-remove, pt+remove,yield+remove,0,error+remove);  
  gerr->Draw("pA");
  gerr->SetName(Form("gErrFonll%s2e",name.Data()));
  g = new TGraph(400-remove,pt+remove,yield+remove);
  g->SetName(Form("gFonll%s2e",name.Data()));
  /* TH1F* h = new TH1F("h","h",1,0,4); */
  /* h->Fill(1); */
  /* h->Draw("p"); */
  h->GetYaxis()->SetRangeUser(1e-2,1e8);
  gPad->SetLogy();
  /* TF1* fpispectra = new TF1("fpi0spectra","[0]*([2]-1)*([2]-2)/([2]*[1]*([2]*[1]+1.864500*([2]-2)))*TMath::Power(([2]*[1]+TMath::Sqrt(x[0]*x[0]+1.864500*1.864500)-1.864500)/([2]*[1]),-[2])*x[0]",0,10); */
  /* TF1* fpispectra = new TF1("fpi0spectra","[0]*([2]-1)*([2]-2)/([2]*[1]*([2]*[1]+1.864500*([2]-2)))*TMath::Power(([2]*[1]+TMath::Sqrt(x[0]*x[0]+1.864500*1.864500)-1.864500)/([2]*[1]),-[2])*x[0]",0,10); */
  /* TF1* fpispectra = new TF1("fpi0spectra","[0]*([2]-1)*([2]-2)/(([2]-3)*([2]-3))/[1]/[1]*TMath::Power(1+2*x[0]/[1]/([2]-3)/2,-1*[2])*x[0]",1,10); */
  /* g->Draw("psame"); */
  /* TF1* fun = new TF1( "fonll62B2e" ,myfun,0,4 , 1 ); */
  TF1* fun = new TF1( Form("fonll62%s2e",name.Data()) ,myfun,0,4 , 1 );
  fun->SetParameter(0,1);
  fun->SetNpx(2000);
  /* fun->Draw("same"); */
  /* fun->SaveAs("test.root"); */
  TFile* f = new TFile(Form("fonll%s2e.root",name.Data()),"recreate");
  fun->Write();
  g->Write();
  h->Write();
  gerr->Write();
}
