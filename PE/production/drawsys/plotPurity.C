#include "myStyle.h"
void plotPurity()
{
   SetMyStyle();
   TFile* fShape = new TFile("Nsigma_2_8_shapefit.root");
   TFile* fDef= new TFile("Nsigma_2_8_purity.root");

   TGraph*  gpurity_ptdef = (TGraph*)fDef->Get("gpurity_ptdef");
   // TGraph*  gpurity_ptdef = (TGraph*)fShape->Get("gptdef")->Clone("gpurity_ptshape");
   TGraph*  gpurity_ptsys = (TGraph*)fDef->Get("gpurity_ptsys");
   // TGraph* gpurity_ptshape = (TGraph*)fShape->Get("gptdef")->Clone("gpurity_ptshape");
   TGraph* gpurity_ptshape = (TGraph*)fShape->Get("gpurity_ptsys")->Clone("gpurity_ptshape");

   fShape->Close();
   fDef->Close();
  
   int nbins=0;
   //add the shape fitting in error consideration only when pt<=0.84
   double x[200],y[200],err1[200],err2[200],errL[200],errH[200];

   for (int ip=0;ip<gpurity_ptdef->GetN();ip++)
   {
     gpurity_ptdef->GetPoint( ip,x[nbins],y[nbins]);
     // cout <<x[nbins] << endl;
     if (x[nbins]>0.2) err1[nbins]=gpurity_ptsys->Eval(x[nbins]);
     else err1[nbins]=y[nbins];
     if (x[nbins]<=3 && x[nbins]>0.2)  {
       err2[nbins] = gpurity_ptshape->Eval(x[nbins]);
     }
     else err2[nbins]=y[nbins];

     errL[nbins]=0;
     errH[nbins]=0;
     // if (err1[nbins]<y[nbins]) errL[nbins] = pow(err1[nbins]-y[nbins], 2);
     // else errH[nbins] = pow(err1[nbins]-y[nbins],2);
     //
     // if (err2[nbins]<y[nbins]) errL[nbins] += pow(err2[nbins]-y[nbins], 2);
     // else errH[nbins] +=pow(err2[nbins]-y[nbins],2);

     // errH[nbins]=sqrt(errH[nbins]);
     // errL[nbins]=sqrt(errL[nbins]);

     errL[nbins] = err1[nbins]<err2[nbins]?err1[nbins]:err2[nbins];
     errL[nbins] = y[nbins]<errL[nbins]?y[nbins]:errL[nbins];
     errL[nbins] = y[nbins]-errL[nbins];
     errH[nbins] = err1[nbins]>err2[nbins]?err1[nbins]:err2[nbins];
     errH[nbins] = y[nbins]>errH[nbins]?y[nbins]:errH[nbins];
     errH[nbins] = errH[nbins]-y[nbins];

     nbins++;  
   }
   
   TGraphAsymmErrors* gpurityerr = new TGraphAsymmErrors(nbins,x,y,0,0,errL,errH);
   gpurityerr->SetName("gpurityerr");
   gpurityerr->SetFillStyle(1001);

   gpurityerr->SetFillColorAlpha(kBlue,0.2);
   gpurityerr->GetXaxis()->SetTitle("p_{T} [GeV/c]");;
   gpurityerr->GetYaxis()->SetTitle("purity");
   gpurityerr->GetYaxis()->SetRangeUser(0.2,1.15);
   gpurityerr->Draw("A3");
   // gpurityerr->Draw("psame");
   gpurity_ptdef->Draw("same");
   gpurity_ptdef->SetName("gpurity_ptdef");
    
   TLatex lat;
   lat.DrawLatexNDC( 0.65,0.45,"0-60%");  
   lat.DrawLatexNDC( 0.6,0.5,"Au+Au 54.4 GeV");  

   TBox* box1 = new TBox(0.4,0.2,0.65,1.1);
   box1->SetFillColorAlpha(kGreen,0.2);
   box1->SetFillStyle(1001);
   box1->Draw();
   TBox* box2 = new TBox(0.7,0.2,1.2,1.1);
   box2->SetFillColorAlpha(kGreen,0.2);
   box2->SetFillStyle(1001);
   box2->Draw();
   TLegend* leg = new TLegend(0.65,0.3,0.83,0.43);
   leg->AddEntry(box1,"excluded","f");
   leg->Draw();
     
   gpurityerr->Draw("3same");
   gpurity_ptdef->Draw("same");

  TFile* final = new TFile("final54purity.root","recreate");
  gpurity_ptdef->Write();
  gpurityerr->Write();


}
