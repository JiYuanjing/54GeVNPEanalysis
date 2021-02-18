#ifndef ST_HELIX_D_HH
#define ST_HELIX_D_HH

#include <math.h>
#include "TVector3.h"
#include "TObject.h"

struct pairD 
{
    double first;
    double second;
    pairD(double, double);
    pairD();
};

inline pairD::pairD(double a, double b) : first(a), second(b) {}
inline pairD::pairD() : first(0), second(0) {}

class StHelixD 
#ifdef __ROOT__
 : public TObject 
#endif
{
public:
    // curvature, dip angle, phase, origin, h
    StHelixD(double c, double dip, double phase,
	     const TVector3& o, int h=-1);
    
    StHelixD();   
    virtual ~StHelixD();
    // StHelixD(const StHelixD&);		// use default
    // StHelixD& operator=(const StHelixD&);	// use default

    double       dipAngle()   const;           
    double       curvature()  const;	// 1/R in xy-plane
    double       phase()      const;	// aziumth in xy-plane measured from ring center
    double       xcenter()    const;	// x-center of circle in xy-plane
    double       ycenter()    const;	// y-center of circle in xy-plane
    int          h()          const;	// -sign(q*B);
    
    const TVector3& origin() const;	// starting point

    void setParameters(double c, double dip, double phase, const TVector3& o, int h);
    
    double       x(double s)  const;
    double       y(double s)  const;
    double       z(double s)  const;

    TVector3 at(double s) const;

    // returns period length of helix
    double       period()       const;
    
    // path length at given r (cylindrical r)
    pairD        pathLength(double r)   const;
    
    // path length at given r (cylindrical r, cylinder axis at x,y)
    pairD        pathLength(double r, double x, double y);

    // path length at distance of closest approach to a given point
    double       pathLength(const TVector3& p) const;
    
    // path length at intersection with plane
    double       pathLength(const TVector3& r,
			    const TVector3& n) const;

    // path length at distance of closest approach in the xy-plane to a given point
    double       pathLength(double x, double y) const;
    
    // path lengths at dca between two helices 
    pairD        pathLengths(const StHelixD&) const;
    
    // minimal distance between point and helix
    double       distance(const TVector3&) const;    
    
    // checks for valid parametrization
    int          valid(double world = 1.e+5) const;
    
    // move the origin along the helix to s which becomes then s=0
    virtual void moveOrigin(double s);
    
protected:
    void setCurvature(double);	// performs also various checks   
    void setPhase(double);	        
    void setDipAngle(double);
    
    // value of S where distance in x-y plane is minimal
    double fudgePathLength(const TVector3&) const;
    
protected:
    int                    mSingularity;	// true for straight line case (B=0)
    TVector3               mOrigin;
    double                 mDipAngle;
    double                 mCurvature;
    double                 mPhase;
    int                    mH;			// -sign(q*B);

    double                 mCosDipAngle;
    double                 mSinDipAngle;
    double                 mCosPhase;
    double                 mSinPhase;

};

//
//     Non-member functions
//
int operator== (const StHelixD&, const StHelixD&);
int operator!= (const StHelixD&, const StHelixD&);
ostream& operator<<(ostream&, const StHelixD&);

//
//     Inline functions
//
inline int StHelixD::h() const {return mH;}

inline double StHelixD::dipAngle() const {return mDipAngle;}

inline double StHelixD::curvature() const {return mCurvature;}

inline double StHelixD::phase() const {return mPhase;}

inline double StHelixD::x(double s) const
{
    if (mSingularity)
	return mOrigin.x() - s*mCosDipAngle*mSinPhase;
    else
	return mOrigin.x() + (cos(mPhase + s*mH*mCurvature*mCosDipAngle)-mCosPhase)/mCurvature;
}
 
inline double StHelixD::y(double s) const
{
    if (mSingularity)
	return mOrigin.y() + s*mCosDipAngle*mCosPhase;
    else
	return mOrigin.y() + (sin(mPhase + s*mH*mCurvature*mCosDipAngle)-mSinPhase)/mCurvature;
}

inline double StHelixD::z(double s) const
{
    return mOrigin.z() + s*mSinDipAngle;
}

inline const TVector3& StHelixD::origin() const {return mOrigin;}

inline TVector3 StHelixD::at(double s) const
{
    return TVector3(x(s), y(s), z(s));
}

inline double StHelixD::pathLength(double x, double y) const
{
    return fudgePathLength(TVector3(x, y, 0));
}

#endif
