#ifndef _LEVEL_SET_STRUCTURE_H_
#define _LEVEL_SET_STRUCTURE_H_

#include "Vector3D.h"

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

    Vec3D totalForce;
};

class DistLevelSetStructure {
  public:
    virtual LevelSetStructure & operator()(int subNum) const = 0;
};

#endif
