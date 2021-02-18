#include <math.h>
#include "StPhysicalHelixD.h"
#include "PhysicalConstants.h" 
#include "SystemOfUnits.h"

StPhysicalHelixD::StPhysicalHelixD() { /* nop */}

StPhysicalHelixD::~StPhysicalHelixD() { /* nop */ }

StPhysicalHelixD::StPhysicalHelixD(const TVector3& p,
				   const TVector3& o,
				   double B, double q)
{
    mH = (q*B <= 0) ? 1 : -1;
    if(p.y() == 0 && p.x() == 0)
	setPhase((M_PI/4)*(1-2.*mH));
    else
	setPhase(atan2(p.y(),p.x())-mH*M_PI/2);
    setDipAngle(atan2(p.z(),p.Perp()));
    mOrigin = o;
    
#ifndef ST_NO_NAMESPACES
    {
	using namespace units;
#endif
	setCurvature(fabs((c_light*nanosecond/meter*q*B/tesla)/
			  (p.Mag()/GeV*mCosDipAngle)/meter));   
#ifndef ST_NO_NAMESPACES
    }
#endif
}

StPhysicalHelixD::StPhysicalHelixD(double c, double d, double phase,
				   const TVector3& o, int h)
    : StHelixD(c, d, phase, o, h) { /* nop */}


TVector3 StPhysicalHelixD::momentum(double B) const
{
    if(mSingularity)
	return(TVector3(0,0,0));
    else {
#ifndef ST_NO_NAMESPACES
	{
	    using namespace units;
#endif
	    double pt =	GeV*fabs(c_light*nanosecond/meter*B/tesla)/(fabs(mCurvature)*meter);
	    
	    return (TVector3(pt*cos(mPhase+mH*M_PI/2),   // pos part pos field
				   pt*sin(mPhase+mH*M_PI/2),
				   pt*tan(mDipAngle)));
#ifndef ST_NO_NAMESPACES
	}
#endif
    }
}

TVector3 StPhysicalHelixD::momentumAt(double S, double B) const
{
    // Obtain phase-shifted momentum from phase-shift of origin
    double xc = this->xcenter();
    double yc = this->ycenter();
    double rx = (y(S)-yc)/(mOrigin.y()-yc);
    double ry = (x(S)-xc)/(mOrigin.x()-xc);
    TVector3 a1 = this->momentum(B);
    TVector3 a(a1.x()*rx, a1.y()*ry, a1.z());
    return a;
}

int StPhysicalHelixD::charge(double B) const
{
    return (B > 0 ? -mH : mH);
}

double StPhysicalHelixD::geometricSignedDistance(double x, double y) 
{
    // Geometric signed distance
    double thePath = this->pathLength(x,y);
    TVector3 DCA2dPosition = this->at(thePath);
    DCA2dPosition.SetZ(0);
    TVector3 position(x,y,0);
    TVector3 DCAVec = (DCA2dPosition-position);
    TVector3 momVec;
    // Deal with straight tracks
    if (this->mSingularity) {
	momVec = this->at(1)- this->at(0);
	momVec.SetZ(0);
    }
    else {
	momVec = this->momentumAt(thePath,1./tesla); // Don't care about Bmag.  Helicity is what matters.
	momVec.SetZ(0);
	
    }
    
    double cross = DCAVec.x()*momVec.y() - DCAVec.y()*momVec.x();
    double theSign = (cross>=0) ? 1. : -1.;
    return theSign*DCAVec.Perp();
}

double StPhysicalHelixD::curvatureSignedDistance(double x, double y)
{
    // Protect against helicity 0 or zero field
    if (this->mSingularity || fabs(this->mH)<=0) {
	return (this->geometricSignedDistance(x,y));
    }
    else {
	return (this->geometricSignedDistance(x,y))/(this->mH);
    }
    
}

double StPhysicalHelixD::geometricSignedDistance(const TVector3& pos)
{
    double sdca2d = this->geometricSignedDistance(pos.x(),pos.y());
    double theSign = (sdca2d>=0) ? 1. : -1.;
    return (this->distance(pos))*theSign;
}

double StPhysicalHelixD::curvatureSignedDistance(const TVector3& pos)
{
    double sdca2d = this->curvatureSignedDistance(pos.x(),pos.y());
    double theSign = (sdca2d>=0) ? 1. : -1.;
    return (this->distance(pos))*theSign;
}
