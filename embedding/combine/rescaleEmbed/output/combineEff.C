#include "TCollection.h"
#include "TList.h"
#include <math>
#include <vector>
void combineEff()
{
  TList* col = new TList();
  TFile* f1 = new TFile("PhoEff_pi0.root");
  TEfficiency* eff1 = (TEfficiency*)f1->Get("EffRecoEff_2_8");
  eff1->SetDirectory(0);
  f1->Close();
  TFile* f2 = new TFile("PhoEff_eta.root");
  TEfficiency* eff2 = (TEfficiency*)f2->Get("EffRecoEff_2_8");
  eff2->SetDirectory(0);
  f2->Close();
  TFile* f3 = new TFile("PhoEff_gamma.root");
  TEfficiency* eff3 = (TEfficiency*)f3->Get("EffRecoEff_2_8");
  eff3->SetDirectory(0);
  f3->Close();
  

  col->AddFirst(eff1);
  col->Add(eff2);
  col->AddLast(eff3);

  TEfficiency* effcomb = new TEfficiency();
  double weight[3]={1,1,1};
  TGraphAsymmErrors* gcomb= effcomb->Combine(col);
  gcomb->Draw();
  //

}
