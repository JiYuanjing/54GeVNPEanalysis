#include "TPythia6Decayer.h"
#include "TPythia6.h"
#include "TParticle.h"
#include "TF1.h"
#include "TH3.h"
#include "TH1.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TProfile.h"
#include "TProfile3D.h"
#include "TProfile2D.h"
#include "TH2.h"
#include "TFile.h"

//constants

namespace Cuts
{
  double M_PION_0 = 0.134977;
  double M_Eta = 0.548;
  double M_D0 = 1.86484;
  int Pi0Id=111;
  int EtaId=221;
  int GammaId = 22;
  int ElectronId = 11;
  int D0Id = 421;
  int nCent = 9;
  std::pair<double,double> momentumrange(0, 10);
  double EtaRange = 1.5; //-1.5,1.5
  int const SpectraParPi0_centbin[9]={0,0,1,1,2,2,3,4,4};
};

TH1F* hEP;
TH1F* hCent;
TPythia6Decayer* pydecay;
TProfile2D* pPi0v2;
TH3D* hPi0v2;
TH2D* hPi0Spec;
TF1* fPi0v2[9];
TF1* fpispectra[9];
TProfile3D* pGammav2;
TH3D* hGammav2[9];
TH3D* hGammaSpec; 

void decayAndFill(int const kf, TLorentzVector* mother, TClonesArray& daughters,float EP ,int cent);
// void fill(TLorentzVector* mother,TLorentzVector* gamma1,TLorentzVector* gamma2, double ptweight, double phiweight,float EP, int cent);
void fill(TLorentzVector* mother,TLorentzVector* gamma1, double ptweight, double phiweight,float EP, int cent);
void fillGamma(TLorentzVector* gamma,double ptweight, double phiweight,float EP,int cent);
// void setDecayChannels(int const defirst,int const desecond,int const mdme);
void setDecayChannels();
void getKinematics(TLorentzVector& b, double const mass);
void getPi0Weight(TLorentzVector* mother,int cent,float EP,double& ptweight,double& phiweight);
float getDeltaPhi(float mphi,float EP);
void bookHists(int mode);
void getCentAndEP(int & cent,float & EP);
void initHists();
void Write(TFile* file);
