#include "Energy.h"
#include "sPhenixStyle.h"
#include "AnaCuts.h"
TH1F* ProjectionAndFit(TString name,int centL,int centH,TString name,TPDF* pdf);
TH1F* getPhev2(TH3F* phe, TH3F* phels, int centL,int centH, double* ptedge,int const nBins,char* hname);
double getHFev2Err(double v1,double v1err,double n1,double n1err,double v2,double v2err,double n2,double n2err,double eff,double eff2err);
double getHFev2ErrWithP(double v1,double v1err,double n1,double n1err,double v2,double v2err,double n2,double n2err,double eff,double eff2err, double effstaterr, double purity,double vk,double kr,double vp,double pr,double vpi,double pir,double vmgpi,double mgpir,double& sys,double& stat,double pt,double hfv2);

double funv2(double*x, double *par);

void SetTH1(TH1* h, int color)
{
   h->SetLineColor(color);
   h->SetMarkerColor(color);
   h->Sumw2();
}
TH1F* hNpheError;
TH1F* hNpheError1;
TH1F* hNpheError2;
TH1F* hNpheError3;
TH1F* hNincError;
TH1F* hphV2error;
TH1F* hincV2error;
TH1F* hS2Bstaterror;
int errorbin = 0;

void HFev2()
{
  getHFev2("plots.pdf","Nsigma_2_8_purity.root","def","out.root");
  getHFev2("plots_syseff.pdf","Nsigma_2_8_purity.root","def","out_syseff.root","fitfun2");
  getHFev2("plots_sys.pdf","Nsigma_2_8_hist.root","sys","out_sys.root");
}

void getHFev2(TString pdfname="plots.pdf",TString purityfilename="Nsigma_2_8_purity.root",TString puritytype="def",TString outfilename="out.root",TString efftype="fitfun2")
{
  int centL = 3,centH=9;
  SetsPhenixStyle(); 
  TFile* file = TFile::Open("incEv2_1012.root"); //fix nHits>20 issue..previous is >=
  TCanvas* c = new TCanvas("c","c");
  TPDF* pdf = new TPDF(pdfname); 
  pdf->Off();
  //book hists
  TFile* fpurity = TFile::Open(purityfilename);
  // TFile* fpurity = TFile::Open("Nsigma_2_8_purity.root");
  TGraph* hpurity = (TGraph*)fpurity->Get(Form("gpurity_pt%s",puritytype.Data())); 
  TGraph* hp = (TGraph*)fpurity->Get(Form("gpratio_pt%s",puritytype.Data())); 
  TGraph* hpi = (TGraph*)fpurity->Get(Form("gkratio_pt%s",puritytype.Data())); 
  TGraph* hk = (TGraph*)fpurity->Get(Form("gpiratio_pt%s",puritytype.Data())); 
  TGraph* hmgpi = (TGraph*)fpurity->Get(Form("gmgpiratio_pt%s",puritytype.Data())); 
  fpurity->Close();
  TFile* fPIDv2 = TFile::Open("chargeparticle/prev2.root");
  TGraphErrors* gKs = (TGraphErrors*)fPIDv2->Get("ks_0_80_62");
  TGraphErrors* gPi = (TGraphErrors*)fPIDv2->Get("pionplus_0_80_62");
  TGraphErrors* gP = (TGraphErrors*)fPIDv2->Get("proton_0_80_62");

  //for sys err of reco eff
  // TFile* fsysReco = TFile::Open("embed/gTotSysErr_Mar06.root");
  /* TFile* fsysReco = TFile::Open("embed/gTotSysErr_scale30.root"); */
  TFile* fsysReco = TFile::Open("embed/gTotSysErr_Jul2.root");
  TGraph* gTotSysErr = (TGraph*)fsysReco->Get("Graph")->Clone("gTotSysErr");

  TFile* freco = TFile::Open("embed/fitRecoEff.root");
  /* TF1* fitreco = (TF1*)freco->Get("fitfun2"); */
  TF1* fitreco = (TF1*)freco->Get(efftype.Data());

  TH1F* hPheEff = ProjectionAndFit("rescaleEmbed/rescale_combine.root", centL-1,centH-1 ,"RecoEff",pdf );
  
  hNpheError = new TH1F("hNpheError","hNpheError",nbin,ptedge);
  hNpheError1 = new TH1F("hNpheError_stat","hNpheError",nbin,ptedge);
  hNpheError3 = new TH1F("hNpheError_stat2","hNpheError",nbin,ptedge);
  hNpheError2 = new TH1F("hNpheError_sys","hNpheError",nbin,ptedge);
  hNincError = new TH1F("hNincError","hNincError",nbin,ptedge);
  hphV2error = new TH1F("hphV2error","hphV2error",nbin,ptedge);
  hincV2error = new TH1F("hincV2error","hincV2error",nbin,ptedge);
  hS2Berror= new TH1F("hS2Berror","hS2Berror",nbin,ptedge);
  hS2Bstaterror= new TH1F("hS2Bstaterror","hS2Berror",nbin,ptedge);
  
  //hists
  TH3F* hphotols = (TH3F*)file->Get("hphoto_LS_hitcut"); 
  // TH3F* hphotols = (TH3F*)file->Get("hphoto_LS"); 
  // TH3F* hphoto = (TH3F*)file->Get("hphoto"); 
  hphotols->Sumw2();
  TH3F* hphoto = (TH3F*)file->Get("hphoto_hitcut"); 
  // TH3F* hphoto = (TH3F*)file->Get("hphoto"); 
  hphoto->Sumw2();
  TH3F* hphotoul = (TH3F*)file->Get("hphoto_hitcut")->Clone("hphotoul"); 
  // TH3F* hphotoul = (TH3F*)file->Get("hphoto")->Clone("hphotoul"); 
  hphoto->Add(hphotols,-1);
  TH1F* centcorr = (TH1F*)file->Get("hcentwg");
  double nEvents = centcorr->Integral(centL,centH);
  cout<<nEvents<<endl;
  TH1F* hPhe = new TH1F("hPhectra", "hPhectra;electron p_{T}(GeV);Counts",nbin,ptedge);
  TH1F* hPheUL = new TH1F("hPheUL", "hPhectra;electron p_{T}(GeV);Counts",nbin,ptedge);
  TH1F* hPheLS = new TH1F("hPheLS", "hPhectra;electron p_{T}(GeV);Counts",nbin,ptedge);
   
  for(int j=0;j<nbin;j++){
    int lbin = hphoto->GetYaxis()->FindBin(ptedge[j]+1e-6);
    int hbin = hphoto->GetYaxis()->FindBin(ptedge[j+1]-1e-6);
    TH1F* hpx = (TH1F*)hphoto->ProjectionX("hpx",lbin,hbin,centL,centH);
    TH1F* hpx_ul = (TH1F*)hphotoul->ProjectionX("hpx_ul",lbin,hbin,centL,centH);
    TH1F* hpxls = (TH1F*)hphotols->ProjectionX("hpx_ls",lbin,hbin,centL,centH);
    SetTH1(hpx_ls,kBlue);
    SetTH1(hpx,kRed);

    hpx_ul->Draw();
    hpx_ul->GetXaxis()->SetTitle("M_{ee} [GeV/c^{2}]");
    hpx_ul->GetYaxis()->SetTitle(Form("Counts/(%0.f MeV/c^{2})",hpx_ul->GetBinWidth(1)*1000));
    hpx_ls->Draw("same");
    hpx->Draw("same");
    TLegend* lpx = new TLegend(0.2,0.7,0.55,0.88);
    lpx->AddEntry(hpx_ul,"Unlike-sign (UL)","l");
    lpx->AddEntry(hpx_ls,"Like-sign (LS)","l");
    lpx->AddEntry(hpx,"UL-LS","l");
    lpx->Draw();
    TLatex lat;
    lat.SetTextSize(0.05);
    lat.DrawLatexNDC(0.55,0.78,"0-60%");
    lat.DrawLatexNDC(0.47,0.3,Form("Tagged e: %.2f<p_{T}<%.2f GeV/c",ptedge[j],ptedge[j+1]));
    lat.DrawLatexNDC(0.55,0.85,"Au+Au #sqrt{s_{NN}} = 54.4 GeV");
    lat.SetTextColor(kGray+2);
    lat.DrawLatexNDC(0.47,0.2,"STAR Preliminary");
    drawLine(0.1,0,0.1,hpx->GetMaximum()*0.5,2,9,kBlack);
    int highbin=hpx->GetXaxis()->FindBin(0.1-1e-6);
    double integralerr;
    hPhe->SetBinContent(j+1,hpx->IntegralAndError(1,highbin,integralerr));
    hPhe->SetBinError(j+1,integralerr);
    hPheUL->SetBinContent(j+1,hpx_ul->Integral(1,highbin));
    hPheLS->SetBinContent(j+1,hpx_ls->Integral(1,highbin));
    // cout <<" integralerr: "<<integralerr<< endl;
    addpdf(pdf);
  }
  // TH2F* hincEptCent= (TH2F*)file->Get("hIncEvsPtvsCent");
  // TH3F* hincEptCent= (TH3F*)file->Get("hIncEv2vsPtvsCent_hitcut");
  TH2F* hincEptCent= (TH2F*)file->Get("hIncEPtvsCent_hitcut");
  TH1F* hincE = (TH1F*)hincEptCent->ProjectionX("hincEpt", centL,centH);
  TH1F* hEsp4eff = (TH1F*)hincE->Clone("hEsp4eff"); // use for reweight the purity 
  hincE = (TH1F*)hincE->Rebin(nbin,"hincE",ptedge);
  hincE->Draw();
  TH1F* hHFe_cor = (TH1F*)hincE->Clone("hHFe_cor");
  hPhe->Draw("same");
  SetTH1(hincE,kRed);
  hincE->GetYaxis()->SetTitle("N/N_{events}");
  hincE->GetYaxis()->SetRangeUser(1e3,3e9);
  TLegend* lince = new TLegend(0.65,0.65,0.85,0.85);
  lince->AddEntry(hincE,"pass cut inclusive e","l");
  lince->AddEntry(hPhe,"reconstructed pho.c e","l");
  lince->Draw();
  gPad->SetLogy(1);
  addpdf(pdf);
  gPad->SetLogy(0);
  TH1F* hratio = (TH1F*)hPhe->Clone("ratio");
  hratio->Divide(hincE);
  hratio->Draw();
  addpdf(pdf);

  // inclusive electron
  // TProfile2D* pE2D = (TProfile2D*)file->Get("pIncEv2");
  TProfile2D* pE2D = (TProfile2D*)file->Get("pIncEv2_hitcut");
  TProfile2D* pPhoE2D = (TProfile2D*)file->Get("pTagEv2");
  TProfile2D* pPhoE2Dls = (TProfile2D*)file->Get("pTagEv2_LS");
  TProfile*  pEv2;
  TProfile*  pPhoEls;
  TProfile*  pPhoEul;
  pEv2 = (TProfile*)pE2D->ProfileX(Form("pEv2_%d_%d",centL,centH),centL,centH); 
  pEv2->SetDirectory(0);
  pEv2 = (TProfile*)pEv2->Rebin(nbin,"pEv2",ptedge);
  pPhoEul = (TProfile*)pPhoE2D->ProfileX(Form("pPhoEul_%d_%d",centL,centH),centL,centH);
  pPhoEul = (TProfile*)pPhoEul->Rebin(nbin,"pPhoEul",ptedge);
  pPhoEls = (TProfile*)pPhoE2Dls->ProfileX(Form("pPhoEls_%d_%d",centL,centH),centL,centH);
  pPhoEls = (TProfile*)pPhoEls->Rebin(nbin,"pPhoEls",ptedge);
  
  //photonic electron
  cout << "photonic electron v2 from measurement" <<endl; 
  TH1F* hphoE = (TH1F*)pPhoEul->ProjectionX("hphoE");
  for (int i=1;i<=hphoE->GetNbinsX();i++)
  {
     // double Nls = pPhoEls->GetBinEntries(i);
     double Nls = hPheLS->GetBinContent(i);
     double v2ls = pPhoEls->GetBinContent(i);
     double v2lserr = pPhoEls->GetBinError(i);
     // double Nul = pPhoEul->GetBinEntries(i);
     double Nul = hPheUL->GetBinContent(i);
     double v2ul = pPhoEul->GetBinContent(i);   
     double v2ulerr = pPhoEul->GetBinError(i);   
     double v2PhE=0;
     if (Nul!=0 && (Nul-Nls)!=0) 
       v2PhE = (Nul*v2ul-Nls*v2ls)*1./(double)(Nul-Nls);
     // cout<<Nls << "  "<< Nul << " " <<v2ls <<" "<<v2ul <<endl;
     hphoE->SetBinContent(i,v2PhE);
      
     double err = v2ulerr*v2ulerr*Nul*Nul+v2lserr*v2lserr*Nls*Nls;
     err = sqrt(err);
     err/=(1.*(Nul-Nls)); 
     hphoE->SetBinError(i,err);
  }

  TFile* fprevious = TFile::Open("phoEv2.root");
  TGraphErrors* gHFe200 = fprevious->Get("HFe200"); 
  // TGraphErrors* gPhoE62 = fprevious->Get("phoE62"); 
  TGraphErrors* gIncE62 = fprevious->Get("incE62"); 
  // TGraphErrors* gPhoE39 = fprevious->Get("phoE39"); 
  // TGraphErrors* gIncE39 = fprevious->Get("incE39"); 
  gIncE62->SetMarkerColor(kMagenta); 
  gIncE62->SetMarkerStyle(20); 
  // gIncE62->Draw("psame");
  // gPhoE62->SetMarkerColor(kOrange+5); 
  // gPhoE62->SetMarkerStyle(20); 
  // gPhoE62->Draw("psame");
  fprevious->Close();
  //
  TFile* fPhev2 = new TFile("PhoEff_comb.root");
  TF1* fitPhev2 = (TF1*)fPhev2->Get("fPhoE_2_8");
  TH1F* hphev2Mc = (TH1F*)fPhev2->Get("hPhoEv2_2_8");
  hphev2Mc->SetDirectory(0);

  fPhev2->Close(); 
  double par62[6];
  fitPhev2->GetParameters(par62);
  TF1* phe62v2 = new TF1("phe62v2","0.963*pol5(0)",0,5);

  phe62v2->SetParameters(par62); 
  phe62v2->Draw();
  hphev2Mc->Draw("same");
  addpdf(pdf);
  TH1F* heratio = (TH1F*)hEsp4eff->Clone("heratio");
  TH1F* hpratio = (TH1F*)hEsp4eff->Clone("hpratio");
  TH1F* hkratio = (TH1F*)hEsp4eff->Clone("hkratio");
  TH1F* hmgpiratio = (TH1F*)hEsp4eff->Clone("hmgpiratio");
  TH1F* hpiratio = (TH1F*)hEsp4eff->Clone("hpiratio");
  TH1F* hreco= (TH1F*)hEsp4eff->Clone("hreco");
  for (int i=1;i<=hEsp4eff->GetNbinsX();i++)
  {
    double pt = hEsp4eff->GetBinCenter(i);
    heratio->SetBinContent(i,hpurity->Eval(pt)*heratio->GetBinContent(i));
    hpratio->SetBinContent(i, hp->Eval(pt)*hpratio->GetBinContent(i));
    hpiratio->SetBinContent(i, hpi->Eval(pt)*hpiratio->GetBinContent(i));
    hkratio->SetBinContent(i, hk->Eval(pt)*hkratio->GetBinContent(i));
    hmgpiratio->SetBinContent(i,hmgpi->Eval(pt)*hmgpiratio->GetBinContent(i));
    hreco->SetBinContent(i,hreco->GetBinContent(i)*fitreco->Eval(hEsp4eff->GetBinCenter(i)));
  }
  TH1F* hmgpiratio = (TH1F*)hmgpiratio->Rebin(nbin,"hmgpiratio",ptedge);
  TH1F* hkratio = (TH1F*)hkratio->Rebin(nbin,"hkratio",ptedge);
  TH1F* hpiratio = (TH1F*)hpiratio->Rebin(nbin,"hpiratio",ptedge);
  TH1F* hpratio = (TH1F*)hpratio->Rebin(nbin,"hpratio",ptedge);
  TH1F* heratio = (TH1F*)heratio->Rebin(nbin,"heratio",ptedge);
  TH1F* hreco = (TH1F*)hreco->Rebin(nbin,"hreco",ptedge);
  TH1F* hEsp4eff = (TH1F*)hEsp4eff->Rebin(nbin,"hEsp4eff",ptedge);
  hmgpiratio->Divide(hEsp4eff);
  hkratio->Divide(hEsp4eff);
  hpiratio->Divide(hEsp4eff);
  hpratio->Divide(hEsp4eff);
  heratio->Divide(hEsp4eff);
  hreco->Divide(hEsp4eff);
  addpdf(pdf);

  hreco->Draw();
  fitreco->Draw("same");
  addpdf(pdf);

  heratio->Draw();
  heratio->GetYaxis()->SetRangeUser(0,1.1);
  hmgpiratio->Draw("same");
  hpiratio->Draw("same");
  hpratio->Draw("same");
  hkratio->Draw("same");
  addpdf(pdf);
  TH1F* hPheCor = (TH1F*)hPhe->Clone("hPheCor");  
  TH1F* hHFv2 = new TH1F("hHFv2","hHFv2",nbin,ptedge);
  hHFv2->SetDirectory(0);
  TH1F* hHFv2sys = new TH1F("hHFv2sys","hHFv2sys",nbin,ptedge);
  hHFv2sys->SetDirectory(0);
  TH1F* hHFv2stat = new TH1F("hHFv2stat","hHFv2stat",nbin,ptedge);
  hHFv2stat->SetDirectory(0);

  hratio = (TH1F*)hincE->Clone("hratioCor");
  for (int ib=0;ib<nbin;ib++) 
  {
    double nphe =  hPhe->GetBinContent(ib+1);
    double pt = hPhe->GetBinCenter(ib+1);
    // nphe/=pheeff->Eval(0.5*(ptedge[ib]+ptedge[ib+1]));
    // double eff =pheeff->Eval(hPhe->GetBinCenter(ib+1));
    /* double eff =pheeff->Eval(ptedge[ib]); */
    double eff = hreco->GetBinContent(ib+1);
    /* double eff = hPheEff->GetBinContent(ib+1); // directly from histogram */
    double effstaterr = hPheEff->GetBinError(ib+1)/eff; // directly from histogram
    // double eff = fitreco->Eval(pt);
    nphe/=eff;
    hPheCor->SetBinContent(ib+1,nphe);
    double nince =  hincE->GetBinContent(ib+1);
    // double phev2 = phe62v2->Eval(hPhe->GetBinCenter(ib+1));
    double autoscale=0.9633;
    /* double autoscale=0.958; */
    double phev2 = hphev2Mc->GetBinContent(ib+1)*autoscale;
    // if (pt>1)  phev2 = hphoE->GetBinContent(ib+1);
    // double phev2err = 0.028;
    // if (pt>1) phev2err = 0.04;
    /* double phev2err = 0.04; */
    double phev2err = 0.03;
    
    double incev2 = pEv2->GetBinContent(ib+1);
    double piv2 = gPi->Eval(hPhe->GetBinCenter(ib+1)); 
    double pv2 = gP->Eval(hPhe->GetBinCenter(ib+1)); 
    double kv2 = gKs->Eval(hPhe->GetBinCenter(ib+1)); 
    double eratio = heratio->GetBinContent(ib+1); 
    double piratio = hpiratio->GetBinContent(ib+1); 
    double pratio = hpratio->GetBinContent(ib+1); 
    double kratio = hkratio->GetBinContent(ib+1); 
    double mgpiratio = hmgpiratio->GetBinContent(ib+1); 

    double hfev2 = (nince*incev2-phev2*nphe);
    hfev2=hfev2-nince*piratio*piv2-nince*kratio*kv2-nince*pratio*pv2-nince*mgpiratio*piv2;
    hfev2/=1.*(eratio*nince-nphe);
   // HF electron yield before single electron efficiency correction
    hHFe_cor->SetBinContent(ib+1,nince*eratio-nphe);
    hratio->SetBinContent(ib+1,nince*eratio);
    // double hfev2err_wo_purity = getHFev2Err(incev2, pEv2->GetBinError(ib+1), nince,0,phev2,phev2err,hPhe->GetBinContent(ib+1), hPhe->GetBinError(ib+1),eff,gTotSysErr->Eval(pt)); 
    double hfev2syserr,hfev2staterr;
    double hfev2err = getHFev2ErrWithP(incev2, pEv2->GetBinError(ib+1), nince,hincE->GetBinError(ib+1), phev2,phev2err,hPhe->GetBinContent(ib+1), hPhe->GetBinError(ib+1),eff,sqrt(pow(gTotSysErr->Eval(pt),2)), effstaterr, eratio,kv2,kratio,pv2,pratio,piv2,piratio,piv2,mgpiratio,hfev2syserr,hfev2staterr,pt,hfev2); 
    /* double hfev2err = getHFev2ErrWithP(incev2, pEv2->GetBinError(ib+1), nince,hincE->GetBinError(ib+1), phev2,phev2err,hPhe->GetBinContent(ib+1), hPhe->GetBinError(ib+1),eff,sqrt(pow(0.035,2)), effstaterr, eratio,kv2,kratio,pv2,pratio,piv2,piratio,piv2,mgpiratio,hfev2syserr,hfev2staterr,pt);  */

    hHFv2->SetBinContent(ib+1,hfev2);
    hHFv2->SetBinError(ib+1,hfev2err);
    hHFv2sys->SetBinContent(ib+1,hfev2);
    hHFv2sys->SetBinError(ib+1,hfev2syserr);
    hHFv2stat->SetBinContent(ib+1,hfev2);
    hHFv2stat->SetBinError(ib+1,hfev2staterr);

  }

  gPad->SetLogy(1);
  SetTH1(hincE,kBlack);
  SetTH1(hPhe,kRed);
  SetTH1(hPheCor,kBlue);

  hincE->Scale(1./nEvents);
  hPheCor->Scale(1./nEvents);
  hPhe->Scale(1./nEvents);
  hratio->Scale(1./nEvents);
  hincE->GetYaxis()->SetRangeUser(1e-6,1);
  hincE->Draw();
  SetTH1(hratio,kGreen+2);
  hratio->Draw("same");
  hPheCor->Draw("same");
  hPhe->Draw("same");
  TLegend* lsp = new TLegend(0.6,0.6,0.85,0.9);
  lsp->AddEntry(hincE,"inc. e","l");
  lsp->AddEntry(hPheCor,"corrected pho. e","l");
  lsp->AddEntry(hPhe,"pho. e","l");
  lsp->Draw();

  addpdf(pdf);
  
  hratio->Add(hPheCor,-1); 
  hratio->Divide(hPheCor); 
  hratio->GetYaxis()->SetTitle("HF e/pho e");
  hratio->Draw();
  addpdf(pdf);

  gPad->SetLogy(0); 
  pEv2->SetMarkerStyle(21);
  pEv2->SetMarkerColor(kRed);
  pEv2->GetYaxis()->SetRangeUser(-0.1,0.2);
  pEv2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  pEv2->GetYaxis()->SetTitle("v_{2}");
  pEv2->Draw();

  hphoE->SetMarkerStyle(21);
  hphoE->SetMarkerColor(kGreen+2);
  hphoE->Draw("psame");
  hphoE->SetDirectory(0);
  phe62v2->Draw("same"); 
  hHFv2->SetMarkerStyle(21);
  hHFv2->SetMarkerColor(kBlue);
  hHFv2->Draw("psame");   
  gHFe200->SetMarkerColor(kMagenta); 
  gHFe200->SetMarkerStyle(20); 
  // gHFe200->Draw("psame");
  //gIncE62->Draw("psame");
  TLegend* l = new TLegend(0.5,0.2,0.8,0.5);
  l->AddEntry(hphoE,"54GeV pho. e v2","pl");
  l->AddEntry(pEv2,"54GeV inc. e v2","pl");
  l->AddEntry(phe62v2,"pho. e v2 embed","l");
  l->AddEntry(hHFv2,"54GeV NPE v2","pl");
  // l->AddEntry(gHFe200,"200GeV NPE v2","pl");
  l->Draw();
  addpdf(pdf);
  
  // hNpheError->SetLineColor(kBlue);
  hNpheError1->SetLineColor(kBlue+1);
  hNpheError3->SetLineColor(kBlack);
  hNpheError2->SetLineColor(kCyan+2);
  hNincError->SetLineColor(kMagenta);
  hphV2error->SetLineColor(kGreen+2);
  hincV2error->SetLineColor(kRed);
  hNpheError1->Draw();
  /* cout << hNpheError1->GetBinContent(hNpheError1->FindBin()); */
  hNpheError1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hNpheError2->Draw("same");
  hNpheError3->Draw("same");
  hNpheError1->GetYaxis()->SetRangeUser(1e-8,2);
  // hNincError->Draw("same");
  hincV2error->Draw("same");
  hphV2error->Draw("same");
  TLegend* lerr = new TLegend(0.45,0.2,0.9,0.5);
  // lerr->AddEntry(hNincError,"Num of inclu. e","l");
  // lerr->AddEntry(hNpheError,"Num of Pho. e","l");
  lerr->AddEntry(hNpheError1,"N^{pho} stat (data)","l");
  lerr->AddEntry(hNpheError3,"N^{pho} stat (embed)","l");
  lerr->AddEntry(hNpheError2,"N^{pho} sys","l");
  lerr->AddEntry(hincV2error,"e^{inc} v_{2}","l");
  lerr->AddEntry(hphV2error,"e^{pho} v_{2}","l");
  gPad->SetLogy();
  lerr->Draw();
  addpdf(pdf);

  TFile* fout = new TFile(outfilename,"recreate");   
  hHFv2->Write();
  hHFv2sys->Write();
  hHFv2stat->Write();
  pEv2->Write();
  hphoE->Write(); 
  hPheCor->Write();
  hPhe->Write();
  hHFe_cor->Write();  
  hratio->Write();
  hS2Berror->Write();
  hS2Bstaterror->Write();
  heratio->Write();
  pdf->On();
  pdf->Close();
  fout->Close();
}

double funv2(double*x, double *par)
{
  return par[0]*(1+par[1]*cos(2*x[0]));
}
TH1F* getIncev2(TH3F* phe,  int centL,int centH, double* ptedge,int const nBins,char* hname,TPDF* pdf)
{
  TH1F* hphev2 = new TH1F(hname,hname,nBins,ptedge);
  TF1* fitfun = new TF1("fitfun","[0]*(1+[1]*cos(2*x))", 0,5); 
  for (int ib=0;ib<nBins;ib++)
  {
    TH1F* h = (TH1F*)phe->ProjectionX(Form("ince%d_%d",ptedge[ib],ptedge[ib+1]),phe->GetYaxis()->FindBin(ptedge[ib]+1e-6),phe->GetYaxis()->FindBin(ptedge[ib+1]-1e-6),centL,centH);
    // TH1F* hls = (TH1F*)phels->ProjectionX(Form("%d_%d_ls",ptedge[ib],ptedge[ib+1]),phe->GetYaxis()->FindBin(ptedge[ib]),phe->GetYaxis()->FindBin(ptedge[ib+1]),centL,centH);
    // h->Add(hls,-1); 
    h->Rebin(20);
    if (ptedge[ib]>2.5) h->Rebin();
    h->Fit(fitfun);
    h->Draw("p");
    drawLatex(0.4,0.8,Form("%0.f<p_{T}<%0.f inc. e",ptedge[ib],ptedge[ib+1]),0.05);
    addpdf(pdf);
    // gPad->Update();
    // gSystem->Sleep(500);
    hphev2->SetBinContent( ib+1, fitfun->GetParameter(1));
    hphev2->SetBinError( ib+1, fitfun->GetParError(1));
  } 
  return hphev2;
}
TH1F* getPhev2(TH3F* phe, TH3F* phels, int centL,int centH, double* ptedge,int const nBins,char* hname,TPDF* pdf)
{
  TH1F* hphev2 = new TH1F(hname,hname,nBins,ptedge);
  TF1* fitfun = new TF1("fitfun","[0]*(1+[1]*cos(2*x))", 0,5); 
  for (int ib=0;ib<nBins;ib++)
  {
    TH1F* h = (TH1F*)phe->ProjectionX(Form("phe%d_%d",ptedge[ib],ptedge[ib+1]),phe->GetYaxis()->FindBin(ptedge[ib]+1e-6),phe->GetYaxis()->FindBin(ptedge[ib+1]-1e-6),centL,centH);
    TH1F* hls = (TH1F*)phels->ProjectionX(Form("%d_%d_ls",ptedge[ib],ptedge[ib+1]),phe->GetYaxis()->FindBin(ptedge[ib]+1e-6),phe->GetYaxis()->FindBin(ptedge[ib+1]-1e-6),centL,centH);
    h->Add(hls,-1); 
    h->Rebin(30);
    if (ptedge[ib]>2.5) h->Rebin();
    h->Fit(fitfun);
    h->Draw("p");
    drawLatex(0.4,0.8,Form("%0.f<p_{T}<%0.f pho. e",ptedge[ib],ptedge[ib+1]),0.05);
    addpdf(pdf);
    // gPad->Update();
    // gSystem->Sleep(500);
    hphev2->SetBinContent( ib+1, fitfun->GetParameter(1));
    hphev2->SetBinError( ib+1, fitfun->GetParError(1));
  } 
  return hphev2;
}
double getHFev2ErrWithP(double v1,double v1err,double n1,double n1err,double v2,double v2err,double n2,double n2err,double eff,double eff2err, double effstaterr, double purity, double vk,double kr,double vp,double pr,double vpi,double pir,double vmgpi,double mgpir,double& sys,double& stat,double pt, double hfv2)
{
  //
  // eff2err and v2err use relative error (photonic electron reco eff and v2)
  double sumerr = 0;
  double hadronv2 = vk*kr+vp*pr+vpi*pir+vmgpi*mgpir;
  double nperr2_1 = pow(n2/eff,2)*(n2err*n2err/(1.*n2*n2)); //stat
  // double nperr2_2 = pow(n2/eff,2)*(eff2err*eff2err);  //stat
  double nperr2_2 = pow(n2/eff,2)*(0*0);  //stat
  double nperr2_3 = pow(n2/eff,2)*(eff2err*eff2err);  //sys
  double nperr2_4 = pow(n2/eff,2)*(effstaterr*effstaterr);  //sys
  // double nperr2_3 = pow(n2/eff,2)*(0.1*0.1);  //sys
  double nperr2 = nperr2_1+nperr2_2;
  // double nperr2 = pow(n2/eff,2)*(n2err*n2err/(1.*n2*n2)+eff2err*eff2err);
  nperr2_1 = nperr2_1*pow(( n1*(v1-hadronv2)-purity*n1*v2 )/(purity*n1-n2/eff)/(purity*n1-n2/eff),2); //num of pho e  //stat
  nperr2_2 = nperr2_2*pow(( n1*(v1-hadronv2)-purity*n1*v2 )/(purity*n1-n2/eff)/(purity*n1-n2/eff),2); //num of pho e // sys
  nperr2_3 = nperr2_3*pow(( n1*(v1-hadronv2)-purity*n1*v2 )/(purity*n1-n2/eff)/(purity*n1-n2/eff),2); //num of pho e // sys
  nperr2_4 = nperr2_4*pow(( n1*(v1-hadronv2)-purity*n1*v2 )/(purity*n1-n2/eff)/(purity*n1-n2/eff),2); //num of pho e // sys
  // nperr2 = nperr2*pow(( n1*(v1-hadronv2)-purity*n1*v2 )/(purity*n1-n2/eff)/(purity*n1-n2/eff),2); //num of pho e
  nperr2 = nperr2_1+nperr2_2+nperr2_3+nperr2_4; //num of pho e
  double vpherr = pow(v2err*v2*(n2/eff)/(purity*n1-n2/eff),2);  // v2 phe
  double vincerr = pow(n1/(purity*n1-n2/eff)*v1err,2);  // v2 inc e
  double nincerr = pow(n1err, 2)*pow((purity*n2/eff*v2-(v1-hadronv2)*n2/eff)/(purity*n1-n2/eff)/(purity*n1-n2/eff),2);
  sumerr = nperr2+vpherr+vincerr+nincerr;
  //v1err stat n2err stat 
  // if (pt<1)
  // {
  sys = sqrt(nperr2_3 + vpherr);
  stat = sqrt(nperr2_1+nperr2_2+nperr2_4 + vincerr + nincerr);
  // }
  // else if (pt>1)
  // {
  //    sys = sqrt(nperr2_2 );
  //    stat = sqrt(nperr2_1 + vpherr + vincerr + nincerr);
  // }

  /* cout<< "the error for each component: "<<endl; */
  /* cout<<"no. of PE: "<< nperr2/sumerr<<endl; */
  /* cout<<"PE v2: "<<vpherr/sumerr <<endl; */
  /* cout<<"no. inc E: " <<nincerr/sumerr<<endl; */
  /* cout<<"v2 inc E: "<< vincerr/sumerr <<endl; */
  // cout<< pow(n2/eff,2)*(n2err*n2err/(1.*n2*n2)+eff2err*eff2err)<< " "<<n1 <<endl;
  cout<<endl;
   
  errorbin++;
  /* hNpheError->SetBinContent(errorbin,nperr2/sumerr); */
  /* hNpheError1->SetBinContent(errorbin,nperr2_1/sumerr); */
  /* hNpheError3->SetBinContent(errorbin,nperr2_4/sumerr); */
  /* hNpheError2->SetBinContent(errorbin,nperr2_3/sumerr); */
  /* hNincError->SetBinContent(errorbin,nincerr/sumerr); */
  /* hphV2error->SetBinContent(errorbin,vpherr/sumerr); */
  /* hincV2error->SetBinContent(errorbin,vincerr/sumerr); */
  hNpheError->SetBinContent(errorbin,sqrt(nperr2)/hfv2);
  hNpheError1->SetBinContent(errorbin,sqrt(nperr2_1)/hfv2);
  cout << "N^{pho} stat (data) " << pt << " "<< sqrt(nperr2_1)/hfv2<< endl;
  hNpheError3->SetBinContent(errorbin,sqrt(nperr2_4)/hfv2);
  cout << "N^{pho} stat (embed) " << pt << " "<< sqrt(nperr2_4)/hfv2<< endl;
  hNpheError2->SetBinContent(errorbin,sqrt(nperr2_3)/hfv2);
  cout << "N^{pho} sys " << pt << " "<< sqrt(nperr2_3)/hfv2<< endl;
  /* hNincError->SetBinContent(errorbin,sqrt(nincerr)/hfv2); */
  /* cout << "N^{incE} sys" << pt << " "<< sqrt(nincerr)/hfv2<< endl; */
  hphV2error->SetBinContent(errorbin,sqrt(vpherr)/hfv2);
  cout << "e^{pho} v_{2} " << pt << " "<< sqrt(vpherr)/hfv2<< endl;
  hincV2error->SetBinContent(errorbin,sqrt(vincerr)/hfv2);
  cout << "e^{inc} v_{2} " << pt << " "<< sqrt(vincerr)/hfv2<< endl;
  
  cout << "e^{HF} v_{2} " << pt << " sys "<< sys/hfv2 << " stat "<<stat/hfv2 << endl;

  hS2Berror->SetBinContent(errorbin, fabs( purity*n1/n2*(eff2err*eff)));
  hS2Bstaterror->SetBinContent(errorbin, fabs( purity*n1/n2*(effstaterr*eff)));
  return sqrt(sumerr);
}

TH1F* ProjectionAndFit(TString filename,int centL,int centH,TString name,TPDF* pdf)
{
  TFile* file = new TFile(filename.Data());
  TH2F* hMc = (TH2F*)file->Get("hTagElectron");
  TH2F* hRc = (TH2F*)file->Get("hTagElectronPassCut");
  int centHbin = hMc->GetYaxis()->FindBin(centH);
  int centLbin = hMc->GetYaxis()->FindBin(centL);
  TH1F* hMcX =(TH1F*)hMc->ProjectionX("hMcX",centLbin,centHbin);
  TH1F* hRcX =(TH1F*)hRc->ProjectionX("hRcX",centLbin,centHbin);
  hRcX = (TH1F*)hRcX->Rebin(nbin,"hRcX",ptedge);
  hMcX = (TH1F*)hMcX->Rebin(nbin,"hMcX" ,ptedge);
  TH1F* hRecoEff = (TH1F*)hRcX->Clone(Form("h%s_%d_%d",name.Data(),centL,centH));
  hRecoEff->SetDirectory(0);
  /* hRecoEff->Divide(hMcX); */
  hRecoEff->GetYaxis()->SetRangeUser(0,0.9);
  TEfficiency*  geff= new TEfficiency(*hRcX,*hMcX);
  // geff->SetStatisticOption(TEfficiency::kBUniform);
  geff->SetStatisticOption(TEfficiency::kBJeffrey);
  geff->SetName(Form("Eff%s_%d_%d",name.Data(),centL,centH));
  for (int ip=0;ip<hRecoEff->GetNbinsX();ip++)
  {
    if (hRecoEff->GetBinCenter(ip+1)<0.2) continue;
    hRecoEff->SetBinContent(ip+1,geff->GetEfficiency(ip+1));
    hRecoEff->SetBinError(ip+1, geff->GetEfficiencyErrorLow(ip+1));
  }
  hRecoEff->Draw();
  // TF1* fit = new TF1(Form("f%s_%d_%d",name.Data(),centL,centH),"pol3",0,4);
  // fit->SetLineColor(kRed);
  hRecoEff->GetXaxis()->SetRangeUser(0.2,4);
  // hRecoEff->Fit(fit,"R");

  drawLatex(0.2,0.8,Form("%s %d-%d \%",name.Data(),CentNum[centH+1],CentNum[centL]), 0.035);
  addpdf(pdf);
  file->Close();
  return hRecoEff;
}
