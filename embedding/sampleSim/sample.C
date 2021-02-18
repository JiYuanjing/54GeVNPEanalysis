#include "TRandom3.h"
void sample()
{
  TF1* funeff = new TF1("funeff","0.5/(TMath::Exp(1.5-x*1.5)+1)",0,4);
  TF1* funspectra = new TF1("funspectra","TMath::Exp(15-x*2)",0,4);
  // TF1* funspectra = new TF1("funspectra","4-(2-x)*(2-x)+2",0,4);
  TH1F* h1 = new TH1F("h1","h1",40,0,4);
  h1->Sumw2();
  TH1F* h2 = new TH1F("h2","h2",40,0,4);
  h2->Sumw2();
  double Nevt = pow(2,32);

  TRandom * rndm = new TRandom3(0);;
  rndm->SetSeed(0);

  double mcpt=-999,rcpt=-999;
  double smean=0.1, ssigma=0.02;
  double smean1=0, ssigma1=0.02;
  double smear=-999,eff=-999;
  for (long i=0;i<Nevt;i++)
  {
    mcpt = funspectra->GetRandom(); 
    eff = funeff->Eval(mcpt); 
    // eff = 0.3; 
    //smearing
    smear=0.;
    if (rndm->Uniform(1)<0.5/(TMath::Exp(1.5-mcpt*1.5)+1)+0.5) 
    {
      if (mcpt>rndm->Gaus(1,0.3)) smear = rndm->Gaus(smean,ssigma);
      else smear=-2;
      eff = 0.3;
    }
    else  smear = rndm->Gaus(smean1,ssigma1);
    // else  smear = -2;
    rcpt = mcpt +smear*mcpt;
    rcpt = mcpt +smear*mcpt;
    h1->Fill(rcpt);
    if (rndm->Uniform(1)<eff) h2->Fill(rcpt);
  }
  h2->Divide(h1);
  h2->Draw();
  // h1->Draw();
  // funeff->Draw();
  // funspectra->Draw();

}

//so from this simulation, we know that smearing is not the reason for the bump,
//even they have 2 components they should have a smooth distribution
//still because these kind of track have some turning in the final efficiency. 
//I think might due to they are parents.
