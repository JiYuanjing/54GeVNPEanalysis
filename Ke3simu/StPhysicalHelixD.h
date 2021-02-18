#ifndef ST_PHYSICAL_HELIX_D_HH
#define ST_PHYSICAL_HELIX_D_HH

#include "TVector3.h"
#include "StHelixD.h"

class StPhysicalHelixD : public StHelixD {
public:
    // Requires: momentum, origin, signed Magnetic Field
    //           and Charge of particle (+/- 1)
    StPhysicalHelixD(const TVector3&,
		     const TVector3&,
		     double, double);
    
    // curvature, dip angle, phase, origin, h
    StPhysicalHelixD(double, double, double,
		     const TVector3&, int h=-1);
    
    StPhysicalHelixD();
    ~StPhysicalHelixD();

    // Requires:  signed Magnetic Field
    TVector3 momentum(double) const;     // returns the momentum at origin
    TVector3 momentumAt(double, double) const; // returns momemtum at S
    int            charge(double)   const;     // returns charge of particle
    // 2d DCA to x,y point signed relative to curvature
    double curvatureSignedDistance(double x, double y) ;
    // 2d DCA to x,y point signed relative to rotation 
    double geometricSignedDistance(double x, double y) ;
    // 3d DCA to 3d point signed relative to curvature
    double curvatureSignedDistance(const TVector3&) ;
    // 3d DCA to 3d point signed relative to rotation
    double geometricSignedDistance(const TVector3&) ;
};
#endif
