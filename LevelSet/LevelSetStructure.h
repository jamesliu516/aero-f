#ifndef _LEVEL_SET_STRUCTURE_H_
#define _LEVEL_SET_STRUCTURE_H_

#include "Vector3D.h"

class Domain;
template<class Scalar, int dim>
class DistSVec;
template <class Scalar> class Vec;
template <class Scalar> class DistVec;

/** Structure used to return levelset information */
struct LevelSetResult {
   Vec3D gradPhi;
   Vec3D normVel;

   LevelSetResult(double gpx, double gpy, double gpz,
                  double nvx, double nvy, double nvz) :
		  gradPhi(gpx, gpy, gpz), normVel(nvx, nvy, nvz) {
		     gradPhi *= 1.0/gradPhi.norm();
		   }

};

/** Abstract class for finding levelset information */
class LevelSetStructure {
  public:
    virtual LevelSetResult
       getLevelSetDataAtEdgeCenter(double t, int ni, int nj) = 0;
    virtual bool isActive(double t, int n) = 0; //!< Whether this node is active or ghost.
    virtual bool edgeIntersectsStructure(double t, int ni, int nj) const = 0;

    virtual void computePhi(Vec<double> &phi);

    Vec3D totalForce;
};

class DistLevelSetStructure {
  protected:
    int numLocSub;
    DistVec<double> *pseudoPhi;
  public:
    virtual ~DistLevelSetStructure() {}

    virtual void initialize(Domain *, DistSVec<double,3> &X) = 0;
    virtual LevelSetStructure & operator()(int subNum) const = 0;
    virtual void clearTotalForce();
    virtual Vec3D getTotalForce();

    virtual DistVec<double> &getPhi();

    double totalForce[3];
};

#endif
