#ifndef StMcAnaCuts_H
#define StMcAnaCuts_H

/* **************************************************
 *  Cuts namespace.
 *
 *  Authors:  Mustafa Mustafa (mmustafa@lbl.gov)
 *
 * **************************************************
 */
#include <vector>

#include "Rtypes.h"
#include "StEvent/StEnumerations.h"

namespace McAnaCuts
{
  std::vector<unsigned int> getAllTriggers()
  {
    std::vector<unsigned int> t;
    t.push_back(580001); 
    t.push_back(580021);

    return t;
  }

  std::vector<unsigned int> const interesting_triggers = getAllTriggers();

  float const mcTrackStartVtxR = 2.0; // maximum
  int const geantId = 2;
  float const TagEta = 0.8;
  float const eta = 0.8;
  float const partnerEta = 1;
  float const gDCA = 1.5;
  int const nHitsFit = 20;
  int const nHitsdEdx= 15;
  float const nFit2nMax = 0.52;
  float const  EEdcaDaughter = 1;
  // std::pair<float,float> phicut={1.75,1.25};
  std::pair<float,float> phicut={1.95,1.25};

  int const SpectraParPi0_centbin[9]={0,0,1,1,2,2,3,4,4};
  // double SpectraParPi0[4][5]={ {1052,0.408,0.091,1.021,10.98},
  //                              {698.6,0.3983,0.064,1.058,11.06},
  //                              {381.4,0.3689,0.03409,1.141,11.42},
  //                              {147.2,0.3094,0.1883,0.927,10.19}
  //                            };
   double const SpectraParPi0[5][5]={ {10.8807, -8.38707, 0.436178, -684.524, 0.519633},
                               {264.787, 0.349324, 0.0475088, 1.01322, 10.903},
                               {406.429, 0.361361, 0.0558035, 1.10099, 11.2414},
                               {682.788, 0.396488, 0.0690007, 1.06078, 11.0674},
                               {992.598, 0.406828, 0.0831998, 1.04464, 11.0797}
                             };
   double const gammaFromPi0[5]={7.6e7, 0.59,0.1,0.71,7.4}; 
   int const nparVz_mult = 7;
   float const parVz_mult[nparVz_mult]={435.9,-0.02413,-0.003707,0.0002204,1.487e-5,-2.95e-07,-1.866e-8};
   int const nCent = 9 ;
   float const Refmult_cent[nCent] = {7,16,31,54,89,138,205,299,361}; //refmult > par[i],  70-80%, 60-70%, ... ,0-5% 
  
 
  int const parentId = 1;
  int const gammaId = -1;  // specfic for gamma embedding, forbid X->gamma->e
  int const dauId1 = 2;
  int const dauId2 = 3; 
  StDedxMethod dedxMethod = kLikelihoodFitId;

  float const minPt=0.1; 
  float const partnerPt = 0.2;
  float const tagPt = 0.2;
  
  float vr = 2;
  float VzVpdvz = 3;
  int const maxNumberOfTriggers = 2;
}
#endif
