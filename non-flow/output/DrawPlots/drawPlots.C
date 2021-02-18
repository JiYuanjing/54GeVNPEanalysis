// min-bias Au+Au hadron
float RefFlowMultAuAu = 281.7;
float RefFlowMultErrAuAu = 7.58e-02;

float RefFlowAuAu = 0.0481507;
float RefFlowErrAuAu = 7.5517e-06;

// pp
const int nEtaGaps=1;
const int nPtBins=7;
float etagap[nEtaGaps]={0.05};
float ptBin[nPtBins+1]={1.2, 1.5, 2.0, 2.5, 3.5, 4.5, 5.5, 8.5};
float Nhadron_ref[nEtaGaps+1] = {0., 0.};
float NhadronErr_ref[nEtaGaps+1] = {0., 0.}; // zyj
//==========================================================================================

void calcM2part(TH1F* dNdCos2DPhi, float& sumCosine, float& sumErrCosine, double& nHadPP, double& nHadErrPP) {
  for(Int_t ibin=0; ibin< dNdCos2DPhi->GetNbinsX(); ibin++)
  {
    sumCosine += dNdCos2DPhi->GetBinCenter(ibin+1) * dNdCos2DPhi->GetBinContent(ibin+1);
    sumErrCosine += dNdCos2DPhi->GetBinCenter(ibin+1) * dNdCos2DPhi->GetBinCenter(ibin+1) * dNdCos2DPhi->GetBinError(ibin+1) * dNdCos2DPhi->GetBinError(ibin+1);
  }
  nHadPP = dNdCos2DPhi->IntegralAndError(1, dNdCos2DPhi->GetNbinsX(), nHadErrPP);
  sumErrCosine = TMath::Sqrt(sumErrCosine);
  //cout << "nHadron_PP +- nHadronErr_PP: " << nHadPP << "+-" << nHadErrPP << endl;
}
//==========================================================================================

void drawPlots(const bool isFullRange=true) {
  gROOT->ProcessLine(".x ./myStyle.C");
  TString inDir;
  if(isFullRange)
    inDir = "fullRange";
  else
    inDir = "nearSide";
  
  TFile *fin_s200GeV  = new TFile(Form("%s/NpeCHcorr_TPC_allpT_lowpT.root", inDir.Data()));
  TFile *fin_s200GeV_fullRange  = new TFile("fullRange/NpeCHcorr_TPC_allpT_lowpT.root");
  
  TString triggerType="pp@200 GeV, Pythia8 Simu.";
  TString trigTitle="Pythia8 Simu";
  
  TH3F::SetDefaultSumw2();
  TH2F::SetDefaultSumw2();
  TH1F::SetDefaultSumw2();
  
  TString paveLabel_system = "PYTHIA8 p+p@200 GeV, STAR HF tuning";
  TString paveLabel_event  = "gg(qqbar)->ccbar->e, 100M c-rich events";
  
  TPaveLabel *pl_system = new TPaveLabel(0.45, 0.78, 0.85, 0.88, paveLabel_system, "brNDC");
  pl_system->SetBorderSize(0);
  pl_system->SetShadowColor(0);
  pl_system->SetTextSize(0.2);
  pl_system->SetTextColor(kBlack);
  pl_system->SetFillStyle(0);
  
  TPaveLabel *pl_event = new TPaveLabel(0.45, 0.68, 0.85, 0.78, paveLabel_event, "brNDC");
  pl_event->SetBorderSize(0);
  pl_event->SetShadowColor(0);
  pl_event->SetTextSize(0.2);
  pl_event->SetTextColor(kBlack);
  pl_event->SetFillStyle(0);
  //==========================================================================================
  //event setting
  //TCanvas *cPtHat = new TCanvas("cPtHat", "", 1200, 600);
  TCanvas *cPtHat = new TCanvas("cPtHat", "");
  cPtHat->Divide(2,1);
  
  cPtHat->cd(1)->SetLogy();
  TH1F *hPTHat = (TH1F*)fin_s200GeV->Get("hPTHat");
  hPTHat->Draw();
  hPTHat->SetStats(kFALSE);
  hPTHat->SetMarkerStyle(2);
  hPTHat->GetXaxis()->SetRangeUser(0,50);
  hPTHat->SetTitle("unweighted #hat{p_{T}}");
  hPTHat->GetXaxis()->SetTitle("#hat{p_{T}} [GeV/c]");
  hPTHat->GetYaxis()->SetTitle("Counts");
  
  TLegend *leg_PtHat1 = new TLegend(0.48,0.58,0.88,0.68);
  leg_PtHat1->SetFillColor(kWhite);
  leg_PtHat1->SetBorderSize(0);
  leg_PtHat1->SetShadowColor(0);
  leg_PtHat1->SetFillStyle(0);
  leg_PtHat1->AddEntry(hPTHat,"#hat{p_{T}}");
  leg_PtHat1->Draw();
  
  pl_system->Draw();
  pl_event->Draw();
  
  cPtHat->cd(2)->SetLogy();
  TH1F *hStatistics = (TH1F*)fin_s200GeV->Get("hStatistics");
  hStatistics->Draw("text00 H");
  hStatistics->SetStats(kFALSE);
  hStatistics->SetTitle("Pythia event statistics");
  hStatistics->GetXaxis()->SetBinLabel(1,"Total cross section (b#bar{b}) [mb]");
  hStatistics->GetXaxis()->SetBinLabel(2,"Number of accepted events");
  
  cPtHat->Print(Form("pythiaEventInfo_%s.pdf", inDir.Data()));
  //==========================================================================================
  //NPE kinematics
  TCanvas *cNpeKine = new TCanvas("cNpeKine", "");
  //TCanvas *cNpeKine = new TCanvas("cNpeKine", "", 2400, 1200);
  cNpeKine->Divide(4,2);
  
  cNpeKine->cd(1)->SetLogy();
  TH1I *hElecMult_gen = (TH1I*)fin_s200GeV->Get("hElecMult_gen");
  hElecMult_gen->Draw();
  hElecMult_gen->SetStats(kFALSE);
  hElecMult_gen->SetMarkerStyle(2);
  hElecMult_gen->GetXaxis()->SetRangeUser(0,100); //zyj
  hElecMult_gen->GetXaxis()->SetTitle("c->e Multiplicity");
  hElecMult_gen->GetYaxis()->SetTitle("Counts");
  
  TLegend *leg_NpeKine1 = new TLegend(0.58,0.58,0.88,0.68);
  leg_NpeKine1->SetFillColor(kWhite);
  leg_NpeKine1->SetBorderSize(0);
  leg_NpeKine1->SetShadowColor(0);
  leg_NpeKine1->SetFillStyle(0);
  leg_NpeKine1->SetHeader("|#eta|<1.0, p_{T}>0."); // ??
  leg_NpeKine1->AddEntry(hElecMult_gen,"c->e");
  leg_NpeKine1->Draw("same");
  
  pl_system->Draw("same");
  pl_event->Draw("same");
  
  cNpeKine->cd(5)->SetLogy();
  TH1I *hElecMult = (TH1I*)fin_s200GeV->Get("hElecMult");
  hElecMult->Draw();
  hElecMult->SetStats(kFALSE);
  hElecMult->SetMarkerStyle(2);
  hElecMult->GetXaxis()->SetRangeUser(0,100);
  hElecMult->GetXaxis()->SetTitle("c->e Multiplicity");
  hElecMult->GetYaxis()->SetTitle("Counts");
  
  TLegend *leg_NpeKine5 = new TLegend(0.58,0.58,0.88,0.68);
  leg_NpeKine5->SetFillColor(kWhite);
  leg_NpeKine5->SetBorderSize(0);
  leg_NpeKine5->SetShadowColor(0);
  leg_NpeKine5->SetFillStyle(0);
  leg_NpeKine5->SetHeader("|#eta|<0.7, p_{T}>1.2 GeV/c");
  leg_NpeKine5->AddEntry(hElecMult, "c->e");
  leg_NpeKine5->Draw();
  
  pl_system->Draw();
  pl_event->Draw();
  
  cNpeKine->cd(2)->SetLogy();
  TH1F *hEtaE_gen = (TH1F*)fin_s200GeV->Get("hEtaE_gen");
  hEtaE_gen->Draw();
  hEtaE_gen->SetStats(kFALSE);
  hEtaE_gen->SetMarkerStyle(2);
  hEtaE_gen->GetXaxis()->SetRangeUser(-1.05,1.05);
  //hEtaE_gen->SetTitle("c->e #eta, |#eta|<1.0, p_{T}>0.2 GeV/c");
  hEtaE_gen->GetXaxis()->SetTitle("#eta^{c->e}");
  hEtaE_gen->GetYaxis()->SetTitle("dN/d#eta");
  
  TLegend *leg_NpeKine2 = new TLegend(0.58,0.58,0.88,0.68);
  leg_NpeKine2->SetFillColor(kWhite);
  leg_NpeKine2->SetBorderSize(0);
  leg_NpeKine2->SetShadowColor(0);
  leg_NpeKine2->SetFillStyle(0);
  leg_NpeKine2->SetHeader("|#eta|<1.0, p_{T}>0. GeV/c");
  leg_NpeKine2->AddEntry(hEtaE_gen, "c->e");
  leg_NpeKine2->Draw();
  
  pl_system->Draw();
  pl_event->Draw();
  
  cNpeKine->cd(6)->SetLogy();
  TH1F *hEtaE = (TH1F*)fin_s200GeV->Get("hEtaE");
  hEtaE->Draw();
  hEtaE->SetStats(kFALSE);
  hEtaE->SetMarkerStyle(2);
  hEtaE->GetXaxis()->SetRangeUser(-1.05,1.05);
  //hEtaE->SetTitle("c->e #eta, |#eta|<1.0, p_{T}>0.2 GeV/c");
  hEtaE->GetXaxis()->SetTitle("#eta^{c->e}");
  hEtaE->GetYaxis()->SetTitle("dN/d#eta");
  
  TLegend *leg_NpeKine6 = new TLegend(0.58,0.58,0.88,0.68);
  leg_NpeKine6->SetFillColor(kWhite);
  leg_NpeKine6->SetBorderSize(0);
  leg_NpeKine6->SetShadowColor(0);
  leg_NpeKine6->SetFillStyle(0);
  leg_NpeKine6->SetHeader("|#eta|<0.7, p_{T}>1.2 GeV/c");
  leg_NpeKine6->AddEntry(hEtaE, "c->e");
  leg_NpeKine6->Draw();
  
  pl_system->Draw();
  pl_event->Draw();
  
  cNpeKine->cd(3)->SetLogy();
  TH1F *hPhiE_gen = (TH1F*)fin_s200GeV->Get("hPhiE_gen");
  hPhiE_gen->Draw();
  hPhiE_gen->SetStats(kFALSE);
  hPhiE_gen->SetMarkerStyle(2);
  hPhiE_gen->GetYaxis()->SetRangeUser(1e5,3e9); // zyj dev
  hPhiE_gen->GetXaxis()->SetTitle("#phi^{c->e}");
  hPhiE_gen->GetYaxis()->SetTitle("dN/d#phi");
  
  TLegend *leg_NpeKine3 = new TLegend(0.58,0.58,0.88,0.68);
  leg_NpeKine3->SetFillColor(kWhite);
  leg_NpeKine3->SetBorderSize(0);
  leg_NpeKine3->SetShadowColor(0);
  leg_NpeKine3->SetFillStyle(0);
  leg_NpeKine3->SetHeader("|#eta|<1.0, p_{T}>0. GeV/c");
  leg_NpeKine3->AddEntry(hPhiE_gen, "c->e");
  leg_NpeKine3->Draw();
  
  pl_system->Draw();
  pl_event->Draw();
  
  cNpeKine->cd(7)->SetLogy();
  TH1F *hPhiE = (TH1F*)fin_s200GeV->Get("hPhiE");
  hPhiE->Draw();
  hPhiE->SetStats(kFALSE);
  hPhiE->SetMarkerStyle(2);
  hPhiE->GetYaxis()->SetRangeUser(1e3,3e9);
  hPhiE->GetXaxis()->SetTitle("#phi^{c->e}");
  hPhiE->GetYaxis()->SetTitle("dN/d#phi");
  
  TLegend *leg_NpeKine7 = new TLegend(0.58,0.58,0.88,0.68);
  leg_NpeKine7->SetFillColor(kWhite);
  leg_NpeKine7->SetBorderSize(0);
  leg_NpeKine7->SetShadowColor(0);
  leg_NpeKine7->SetFillStyle(0);
  leg_NpeKine7->SetHeader("|#eta|<0.7, p_{T}>1.2 GeV/c");
  leg_NpeKine7->AddEntry(hPhiE, "c->e");
  leg_NpeKine7->Draw();
  
  pl_system->Draw();
  pl_event->Draw();
  
  cNpeKine->cd(4)->SetLogy();
  TH1F *hPTTrkE_gen = (TH1F*)fin_s200GeV->Get("hPTTrkE_gen");
  hPTTrkE_gen->Draw();
  hPTTrkE_gen->SetStats(kFALSE);
  hPTTrkE_gen->SetMarkerStyle(2);
  hPTTrkE_gen->GetXaxis()->SetRangeUser(0,20);
  hPTTrkE_gen->GetXaxis()->SetTitle("p^{c->e}_{T} [GeV/c]");
  hPTTrkE_gen->GetYaxis()->SetTitle("dN/dp_{T}");
  
  TLegend *leg_NpeKine4 = new TLegend(0.58,0.58,0.88,0.68);
  leg_NpeKine4->SetFillColor(kWhite);
  leg_NpeKine4->SetBorderSize(0);
  leg_NpeKine4->SetShadowColor(0);
  leg_NpeKine4->SetFillStyle(0);
  leg_NpeKine4->SetHeader("|#eta|<1.0, p_{T}>0. GeV/c");
  leg_NpeKine4->AddEntry(hPTTrkE_gen, "c->e");
  leg_NpeKine4->Draw();
  
  pl_system->Draw();
  pl_event->Draw();
  
  cNpeKine->cd(8)->SetLogy();
  TH1F *hPTTrkE = (TH1F*)fin_s200GeV->Get("hPTTrkE");
  hPTTrkE->Draw();
  hPTTrkE->SetStats(kFALSE);
  hPTTrkE->SetMarkerStyle(2);
  hPTTrkE->GetXaxis()->SetRangeUser(0,20);
  hPTTrkE->GetXaxis()->SetTitle("p^{c->e}_{T} [GeV/c]");
  hPTTrkE->GetYaxis()->SetTitle("dN/dp_{T}");
  
  TLegend *leg_NpeKine8 = new TLegend(0.58,0.58,0.88,0.68);
  leg_NpeKine8->SetFillColor(kWhite);
  leg_NpeKine8->SetBorderSize(0);
  leg_NpeKine8->SetShadowColor(0);
  leg_NpeKine8->SetFillStyle(0);
  leg_NpeKine8->SetHeader("|#eta|<0.7, p_{T}>1.2 GeV/c");
  leg_NpeKine8->AddEntry(hPTTrkE, "c->e");
  leg_NpeKine8->Draw();
  
  pl_system->Draw();
  pl_event->Draw();
  
  cNpeKine->cd();
  cNpeKine->Print(Form("npeKinematics_%s.pdf", inDir.Data()));
  //=================================================================================

  // zyj
  double ng = hElecMult_gen->GetEntries();
  double ne = hPTHat->GetEntries();
  cout<<"number of events with HF->e is: "<<ne<<endl;
  cout<<"number of generate events is: "<<ng<<endl;
  //=================================================================================
  // pp hadron dist. as a function of pt
  TCanvas *cPtHadron = new TCanvas("cPtHadron", "");
  //TCanvas *cPtHadron = new TCanvas("cPtHadron", "", 1200, 600);

  cPtHadron->Divide(2);
  TH1F *hHadron_Pt_PP[2];
  hHadron_Pt_PP[0]= (TH1F*)fin_s200GeV->Get("hHad_noEtaGap_pTBin7");
  hHadron_Pt_PP[1]= (TH1F*)fin_s200GeV->Get("hHad_EtaGap_dot0_pTBin7");
  cPtHadron->cd(1);
  hHadron_Pt_PP[0]->Draw();
  cPtHadron->cd(2);
  hHadron_Pt_PP[1]->Draw();
  cPtHadron->SaveAs("associate_hadron_pp.pdf");
  
  cout<<"number of HF->e in noEtagap is: "<<hHadron_Pt_PP[0]->GetMean()<<", error is: "<<hHadron_Pt_PP[0]->GetMeanError()<<endl;
  cout<<"number of HF->e in Etagap1 is:  "<<hHadron_Pt_PP[1]->GetMean()<<", error is: "<<hHadron_Pt_PP[1]->GetMeanError()<<endl;
  
  for(auto i=0; i<nEtaGaps+1; i++){
    Nhadron_ref[i] = hHadron_Pt_PP[i]->GetMean();
    NhadronErr_ref[i] = hHadron_Pt_PP[i]->GetMeanError();
  }
  //=================================================================================
  //deltaPhi dist.
  TH1F *hNumTrigNPE_allTracks = (TH1F*)fin_s200GeV_fullRange->Get("hNumTrigNPE_allTracks");
  TH1F *hNumTrigNPE_deltaEta1 = (TH1F*)fin_s200GeV_fullRange->Get("hNumTrigNPE_deltaEta1");
  
  //TCanvas *cDPhi = new TCanvas("cDPhi", "", 1800, 1200);
  TCanvas *cDPhi = new TCanvas("cDPhi", "");

  cDPhi->Divide(2,2);
  
  //double numofNPE_PtBin0 = hPTTrkE->Integral(21, 25);
  //cout << "Number of trigger NPEs at PtBin0: " << numofNPE_PtBin0 << " vs " <<  hNumTrigNPE_allTracks->GetBinContent(1) << endl;
  TString ptRange[nPtBins+1]={"1.2 < p^{NPE}_{T} < 1.5GeV/c", "1.5 < p^{NPE}_{T} < 2.0 GeV/c", "2.0 < p^{NPE}_{T} < 2.5 GeV/c", "2.5 < p^{NPE}_{T} < 3.5 GeV/c", "3.5 < p^{NPE}_{T} < 4.5 GeV/c", "4.5 < p^{NPE}_{T} < 5.5 GeV/c", "5.5 < p^{NPE}_{T} < 8.5 GeV/c", "2.0 < p^{NPE}_{T} < 8.5 GeV/c"};
  TH1F *hDPhi_allTracks_PtBin[nPtBins+1];
  TH1F *hDPhi_deltaEta01_PtBin[nPtBins+1];
  TLegend *leg_dPhi[nPtBins+1];
  //for(auto i=0; i<nPtBins+1; i++){
  for(auto i=0; i<4; i++){
    cDPhi->cd(i+1);
    hDPhi_allTracks_PtBin[i] = (TH1F*)fin_s200GeV->Get(Form("hDPhi_allTracks_PtBin%d", i));
    hDPhi_allTracks_PtBin[i]->SetStats(kFALSE);
    hDPhi_allTracks_PtBin[i]->Scale(1.0/hNumTrigNPE_allTracks->GetBinContent(i+1));
    //hDPhi_allTracks_PtBin[i]->Draw("E1");
    hDPhi_allTracks_PtBin[i]->SetMarkerStyle(24);
    hDPhi_allTracks_PtBin[i]->SetMarkerSize(0.5);
    hDPhi_allTracks_PtBin[i]->SetMarkerColor(kBlack);
    hDPhi_allTracks_PtBin[i]->SetTitle("");
    hDPhi_allTracks_PtBin[i]->GetYaxis()->SetRangeUser(0., 0.12);
    hDPhi_allTracks_PtBin[i]->GetXaxis()->SetTitle("#Delta#phi");
    hDPhi_allTracks_PtBin[i]->GetYaxis()->SetTitle("1/N_{NPE} dN/d#Delta#phi");
    hDPhi_allTracks_PtBin[i]->SetTitleOffset(1.5,"Y");
    
    hDPhi_deltaEta01_PtBin[i]= (TH1F*)fin_s200GeV->Get(Form("hDPhi_deltaEta01_PtBin%d", i));
    hDPhi_deltaEta01_PtBin[i]->SetStats(kFALSE);
    hDPhi_deltaEta01_PtBin[i]->Scale(1.0/hNumTrigNPE_deltaEta1->GetBinContent(i+1));
    hDPhi_deltaEta01_PtBin[i]->Draw("E1");
    hDPhi_deltaEta01_PtBin[i]->SetMarkerStyle(25);
    hDPhi_deltaEta01_PtBin[i]->SetMarkerSize(0.5);
    hDPhi_deltaEta01_PtBin[i]->SetMarkerColor(kRed);
    hDPhi_deltaEta01_PtBin[i]->SetTitle("");
    hDPhi_deltaEta01_PtBin[i]->GetYaxis()->SetRangeUser(0., 0.12);
    hDPhi_deltaEta01_PtBin[i]->GetXaxis()->SetTitle("#Delta#phi");
    hDPhi_deltaEta01_PtBin[i]->GetYaxis()->SetTitle("1/N_{NPE} dN/d#Delta#phi");
    hDPhi_deltaEta01_PtBin[i]->SetTitleOffset(1.5,"Y");
    
    
    leg_dPhi[i] = new TLegend(0.28,0.58,0.78,0.88);
    leg_dPhi[i]->SetFillColor(kWhite);
    leg_dPhi[i]->SetBorderSize(0);
    leg_dPhi[i]->SetShadowColor(0);
    leg_dPhi[i]->SetFillStyle(0);
    leg_dPhi[i]->SetHeader(Form("%s, |#eta|<0.7", ptRange[i].Data()));
    //leg_dPhi[i]->AddEntry(hDPhi_allTracks_PtBin[i],"w/o eta gap, all hadrons");
    leg_dPhi[i]->AddEntry(hDPhi_deltaEta01_PtBin[i],"w/ eta gap, #Delta#eta>0.05");
    leg_dPhi[i]->Draw();
  }
  
  cDPhi->Print(Form("dPhi_%s.pdf", inDir.Data()));
  //=================================================================================
  //Cos2DPhi
  TH1F *hCos2DPhi_allTracks_PtBin[nPtBins+1];
  TH1F *hCos2DPhi_deltaEta01_PtBin[nPtBins+1];
  for(auto i=0; i<nPtBins+1; i++){
    hCos2DPhi_allTracks_PtBin[i] = (TH1F*)fin_s200GeV->Get(Form("hCos2DPhi_allTracks_PtBin%d", i));
    hCos2DPhi_allTracks_PtBin[i]->Scale(1.0/hNumTrigNPE_allTracks->GetBinContent(1));
    
    hCos2DPhi_deltaEta01_PtBin[i] = (TH1F*)fin_s200GeV->Get(Form("hCos2DPhi_deltaEta01_PtBin%d", i));
    hCos2DPhi_deltaEta01_PtBin[i]->Scale(1.0/hNumTrigNPE_deltaEta1->GetBinContent(1));
  }
  
  //calculate M2 and non-flow
  float sumCos2DPhi[nPtBins+1][nEtaGaps+1];
  float sumErrCos2DPhi[nPtBins+1][nEtaGaps+1];
  double nHadronPP[nPtBins+1][nEtaGaps+1];
  double nHadronErrPP[nPtBins+1][nEtaGaps+1];
  for(int iPtBin=0; iPtBin<nPtBins+1; iPtBin++) {
    for(int iGap=0; iGap<nEtaGaps+1; iGap++) {
      sumCos2DPhi[iPtBin][iGap] = 0.;
      sumErrCos2DPhi[iPtBin][iGap] = 0.;
      nHadronPP[iPtBin][iGap] = 0.;
      nHadronErrPP[iPtBin][iGap] = 0.;
    }
  }
  for(int i=0; i<nPtBins+1; i++) {
      calcM2part(hCos2DPhi_allTracks_PtBin[i], sumCos2DPhi[i][0], sumErrCos2DPhi[i][0], nHadronPP[i][0], nHadronErrPP[i][0]);
      calcM2part(hCos2DPhi_deltaEta01_PtBin[i], sumCos2DPhi[i][1], sumErrCos2DPhi[i][1], nHadronPP[i][1], nHadronErrPP[i][1]);
  }
  //=================================================================================
  //draw M2 and non-flow
  TH1F *hNHadron_refAA = new TH1F("hNHadron_refAA","",nPtBins, ptBin);
  hNHadron_refAA->Sumw2();
  TH1F *hV2_refAA = new TH1F("hV2_refAA","",nPtBins, ptBin);
  hV2_refAA->Sumw2();
  TH1F *hNHadron_refPP[nEtaGaps+1];
  for(auto i=0; i<nEtaGaps+1; i++){
    hNHadron_refPP[i] = new TH1F(Form("hNHadron_refPP%d", i),"",nPtBins, ptBin);
    hNHadron_refPP[i]->Sumw2();
  }

  for(int iPtBin=0; iPtBin<nPtBins; iPtBin++) {
    hNHadron_refAA->SetBinContent(iPtBin+1, RefFlowMultAuAu);
    hNHadron_refAA->SetBinError(iPtBin+1, RefFlowMultErrAuAu);
    
    hV2_refAA->SetBinContent(iPtBin+1, RefFlowAuAu);
    hV2_refAA->SetBinError(iPtBin+1, RefFlowErrAuAu);
    
    hNHadron_refPP[0]->SetBinContent(iPtBin+1, Nhadron_ref[0]);
    hNHadron_refPP[0]->SetBinError(iPtBin+1, NhadronErr_ref[0]);
    
    hNHadron_refPP[1]->SetBinContent(iPtBin+1, Nhadron_ref[1]);
    hNHadron_refPP[1]->SetBinError(iPtBin+1, NhadronErr_ref[1]);
  }
  
  TH1F *hNHadron_PP[nEtaGaps+1];
  TH1F *hM2[nEtaGaps+1];
  TH1F *hV2[nEtaGaps+1];
  TString hName[5] = {"noEtaGap","EtaGap_dot0","EtaGap_dot1","EtaGap_dot2","EtaGap_dot4"};
  int mMarkerColor[5] = {1,2,3,4,6};
  
  for(int iGap=0; iGap<nEtaGaps+1; iGap++) {
    hNHadron_PP[iGap] = new TH1F(Form("hNHadron_PP_%s",hName[iGap].Data()),"",nPtBins, ptBin);
    hNHadron_PP[iGap]->Sumw2();
    
    hM2[iGap] = new TH1F(Form("hM2_%s",hName[iGap].Data()),"",nPtBins, ptBin);
    hM2[iGap]->Sumw2();
    
    hV2[iGap] = new TH1F(Form("hV2_%s",hName[iGap].Data()),"",nPtBins, ptBin);
    hV2[iGap]->Sumw2();
    
    //for(int iPtBin=0; iPtBin<nPtBins; iPtBin++) {
    for(int iPtBin=0; iPtBin<4; iPtBin++) {
      hNHadron_PP[iGap]->SetBinContent(iPtBin+1, nHadronPP[iPtBin][iGap]);
      hNHadron_PP[iGap]->SetBinError(iPtBin+1, nHadronErrPP[iPtBin][iGap]);

      hM2[iGap]->SetBinContent(iPtBin+1, sumCos2DPhi[iPtBin][iGap]);
      hM2[iGap]->SetBinError(iPtBin+1, sumErrCos2DPhi[iPtBin][iGap]);
    }
  }

  //float downLimit_M2 = 1e-4, upLimit_M2 = 0.05;
  //float downLimit_V2 = 2e-4, upLimit_V2 = 0.05;
  //float downLimit_M2 = 1e-2, upLimit_M2 = 0.3;
  //float downLimit_V2 = 2e-3, upLimit_V2 = 0.3;
  float downLimit_M2 = 0, upLimit_M2 = 0.25;
  float downLimit_V2 = 0, upLimit_V2 = 0.25;
  if(!isFullRange) {
    downLimit_M2 = 5e-4;
    upLimit_M2 = 1.1;
    downLimit_V2 = 5e-4;
    upLimit_V2 = 1.1;
  }
  
  TCanvas *caM2 = new TCanvas("caM2", "", 1800, 600);
  //TCanvas *caM2 = new TCanvas("caM2", "");
  caM2->Divide(2,1);
  
  caM2->cd(1);
  //caM2->cd(1)->SetLogy();
  hM2[0]->Divide(hNHadron_PP[0]);
  hM2[0]->SetStats(kFALSE);
  hM2[0]->GetYaxis()->SetRangeUser(downLimit_M2, upLimit_M2);
  hM2[0]->GetXaxis()->SetRangeUser(1., 5.5);
  //hM2[0]->Draw("E1");
  hM2[0]->SetMarkerStyle(24);
  hM2[0]->SetMarkerSize(1.0);
  hM2[0]->SetMarkerColor(kBlack);

  hM2[0]->SetTitle("");
  hM2[0]->GetXaxis()->SetTitle("p^{NPE}_{T} [GeV/c]");
  hM2[0]->GetYaxis()->SetTitle("M^{pp}_{2} (= #Sigma cos[2(#phi_{NPE} - #phi_{h})]/<N^{pp}_{NPE}> / <N^{pp}_{h}>)");
  hM2[0]->SetTitleOffset(1.5,"Y");

  hM2[1]->Divide(hNHadron_PP[1]);
  hM2[1]->GetYaxis()->SetRangeUser(downLimit_M2, upLimit_M2);
  hM2[1]->GetXaxis()->SetRangeUser(1, 5.5);
  hM2[1]->SetStats(kFALSE);
  hM2[1]->Draw("E1");
  hM2[1]->SetMarkerStyle(25);
  hM2[1]->SetMarkerSize(1.0);
  hM2[1]->SetMarkerColor(kRed);
  hM2[1]->SetTitle("");
  hM2[1]->GetXaxis()->SetTitle("p^{NPE}_{T} [GeV/c]");
  hM2[1]->GetYaxis()->SetTitle("M^{pp}_{2} (= #Sigma cos[2(#phi_{NPE} - #phi_{h})]/<N^{pp}_{NPE}> / <N^{pp}_{h}>)");
  hM2[1]->SetTitleOffset(1.5,"Y");
  
  //TLegend *leg_M2 = new TLegend(0.15,0.55,0.55,0.85);
  TLegend *leg_M2 = new TLegend(0.45,0.15,0.95,0.35);
  leg_M2->SetFillColor(kWhite);
  leg_M2->SetBorderSize(0);
  leg_M2->SetShadowColor(0);
  leg_M2->SetFillStyle(0);
  leg_M2->SetHeader("0.15 < p^{h}_{T} < 2.0 GeV/c, |#eta_{h}|<1");
  //leg_M2->AddEntry(hM2[0],"w/o eta gap");
  leg_M2->AddEntry(hM2[1],"w/ eta gap, #Delta#eta>0.05");
  leg_M2->Draw();
  
  for(int iGap=0; iGap<nEtaGaps+1; iGap++) {
    hV2[iGap]=(TH1F*)hM2[iGap]->Clone();
  }
  
  caM2->cd(2);
  //caM2->cd(2)->SetLogy();
  hV2[0]->Multiply(hNHadron_refPP[0]);
  hV2[0]->Divide(hNHadron_refAA);
  hV2[0]->Divide(hV2_refAA);
  hV2[0]->SetStats(kFALSE);
  hV2[0]->GetYaxis()->SetRangeUser(downLimit_V2, upLimit_V2);
  hV2[0]->GetXaxis()->SetRangeUser(1, 5.5);
  //hV2[0]->Draw("E1");
  hV2[0]->SetMarkerStyle(24);
  hV2[0]->SetMarkerSize(1.0);
  hV2[0]->SetMarkerColor(kBlack);

  hV2[0]->SetTitle("");
  hV2[0]->GetXaxis()->SetTitle("p^{NPE}_{T} [GeV/c]");
  hV2[0]->GetYaxis()->SetTitle("#nu^{e,AA}_{2} (= M^{pp}_{2} <N^{pp}_{h}> / v^{h, AA}_{2} / <N^{AA}_{h}>)");
  hV2[0]->SetTitleOffset(1.5,"Y");
  
  hV2[1]->Multiply(hNHadron_refPP[1]);
  hV2[1]->Divide(hNHadron_refAA);
  hV2[1]->Divide(hV2_refAA);
  hV2[1]->GetYaxis()->SetRangeUser(downLimit_V2, upLimit_V2);
  hV2[1]->GetXaxis()->SetRangeUser(1, 5.5);
  hV2[1]->SetStats(kFALSE);
  hV2[1]->Draw("E1");
  hV2[1]->SetMarkerStyle(25);
  hV2[1]->SetMarkerSize(1.0);
  hV2[1]->SetMarkerColor(kRed);
  hV2[1]->SetTitle("");
  hV2[1]->GetXaxis()->SetTitle("p^{NPE}_{T} [GeV/c]");
  hV2[1]->GetYaxis()->SetTitle("#nu^{e,AA}_{2} (= M^{pp}_{2} <N^{pp}_{h}> / v^{h, AA}_{2} / <N^{AA}_{h}>)");
  hV2[1]->SetTitleOffset(1.5,"Y");
  
  //TLegend *leg_V2 = new TLegend(0.15,0.55,0.55,0.85);
  TLegend *leg_V2 = new TLegend(0.45,0.15,0.95,0.35);
  leg_V2->SetFillColor(kWhite);
  leg_V2->SetBorderSize(0);
  leg_V2->SetShadowColor(0);
  leg_V2->SetFillStyle(0);
  leg_V2->SetHeader("0.15 < p^{h}_{T} < 2.0 GeV/c, |#eta_{h}|<1");
  //leg_V2->AddEntry(hV2[0],"w/o eta gap");
  leg_V2->AddEntry(hV2[1],"w/ eta gap, #Delta#eta>0.05");
  leg_V2->Draw();
  
  caM2->cd();
  caM2->Print(Form("m2v2_%s.pdf", inDir.Data()));
  // print out the result
  for(auto i=0; i<nPtBins+1; i++){
    cout<<ptBin[i]<<"   "<<hV2[1]->GetBinContent(i+1)<<"   "<<hV2[1]->GetBinError(i+1)<<endl;
  }
  
}
void drawtext(TH1F *name)
{
  Int_t i,n;
  Double_t x,y;
  TLatex *l;
  
  TH1F *g = (TH1F*)gPad->GetListOfPrimitives()->FindObject(name);
  for (i=0; i<5; i++) {
    l = new TLatex(g->GetBinCenter(i+1),g->GetBinContent(i+1)+0.01,Form("%.3f",y));
    l->SetTextSize(0.025);
    l->SetTextFont(42);
    l->SetTextAlign(21);
    l->Paint();
  }
}


