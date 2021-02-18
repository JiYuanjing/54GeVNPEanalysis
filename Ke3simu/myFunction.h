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
#include "TGraphAsymmErrors.h"
//constants

namespace Cuts
{
  double M_PION_0 = 0.134977;
  double M_Eta = 0.548;
  double M_D0 = 1.86484;
  double M_Kaon = 0.49367;
  int Pi0Id=111;
  int EtaId=221;
  int GammaId = 22;
  int ElectronId = 11;
  int D0Id = 421;
  int KaonId = 321;
  int nCent = 9;
  std::pair<double,double> momentumrange(0.3, 5);
  double EtaRange = 1.5; //-1.5,1.5
  int const SpectraParPi0_centbin[9]={0,0,1,1,2,2,3,4,4};
};

TH1F* hEP;
TH1F* hCent;
TPythia6Decayer* pydecay;
TProfile3D* pKaonV2;
// TH3D* hPi0v2;
TH3D* hKaonSpec;
TH1F* hDecayL;
TF1* fPi0v2[9];
TGraphAsymmErrors* gpispectra[9];
TF1* fpispectra1[9];
TF1* fpispectra2[9];
TProfile3D* pElectronV2NoCut;
// TH3D* hGammav2[9];
TH3D* hElectronSpecNoCut; 
TProfile3D* pElectronV2DcaCut;
// TH3D* hGammav2[9];
TH3D* hElectronSpecDcaCut; 

double funK0s(double *x, double *par);
void decayAndFill(int const kf, TLorentzVector* mother, TClonesArray& daughters,float EP ,int cent, int mode);
// void fill(TLorentzVector* mother,TLorentzVector* gamma1,TLorentzVector* gamma2, double ptweight, double phiweight,float EP, int cent);
// void fill(TLorentzVector* mother,TLorentzVector* gamma1, double ptweight, double phiweight,float EP, int cent);
void fillGamma(TLorentzVector* gamma,TVector3 dca, TVector3 v0,double ptweight, double phiweight,float EP,int cent);
// void setDecayChannels(int const defirst,int const desecond,int const mdme);
TVector3 smearPos(TLorentzVector const* rMom, TVector3 const& pos);
void setDecayChannels(int mode);
void getKinematics(TLorentzVector& b, double const mass);
void getPi0Weight(TLorentzVector* mother,int cent,float EP,double& ptweight,double& phiweight);
float getDeltaPhi(float mphi,float EP);
void bookHists(int mode);
void getCentAndEP(int & cent,float & EP);
void initHists();
void Write(TFile* file);
