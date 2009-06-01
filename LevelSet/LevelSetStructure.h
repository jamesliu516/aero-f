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
   double alpha;
   double xi[3];
   int trNodes[3];
   Vec3D gradPhi;
   Vec3D normVel;

   LevelSetResult() {
     alpha = xi[0] = xi[1] = -1.0;
     trNodes[0] = trNodes[1] = trNodes[2] = -1;
     gradPhi = normVel = 0.0;
   }
   LevelSetResult(double gpx, double gpy, double gpz,
                  double nvx, double nvy, double nvz) :
		  gradPhi(gpx, gpy, gpz), normVel(nvx, nvy, nvz) {
		     gradPhi *= 1.0/gradPhi.norm();
                     alpha = -1.0;
                     xi[0] = xi[1] = xi[2] = -1.0;
                     trNodes[0] = trNodes[1] = trNodes[2] = -1;
		   }

   class iterator {
	   double *xip;
	   int *nodep;
   public:
	   iterator(double *x, int *n) : xip(x), nodep(n) { }
	   iterator & operator++() { xip++; nodep++; }
	   bool operator !=(const iterator &it) const {
		   return xip != it.xip || nodep != it.nodep;
	   }
	   double Ni() const { return *xip; }
	   int nodeNum() const { return *nodep; }
   };

   iterator begin() { return iterator(xi, trNodes); }
   iterator end() { return iterator(xi+3, trNodes+3); }
};

/** Abstract class for finding levelset information */
class LevelSetStructure {
  public:
    /** returns the normal and normal velocity at intersection between edge ni, nj and structure
     *
     * If ni, nj is not an edge of the fluid mesh, result is undefined.
     * */
    virtual LevelSetResult
       getLevelSetDataAtEdgeCenter(double t, int ni, int nj) = 0;
    virtual bool isActive(double t, int n) const = 0; //!< Whether this node is active or ghost.
    virtual bool edgeIntersectsStructure(double t, int ni, int nj) const = 0; //!< whether an edge between i and j intersects the structure

    /** creates an array of values which are positive inside the fluid and negative outside. */
    virtual void computePhi(Vec<double> &phi);

    virtual double isPointOnSurface(Vec3D, int, int, int) = 0;

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
    virtual Vec3D getTotalForce(const double pref);

    virtual DistVec<double> &getPhi();

    virtual Vec<Vec3D> &getStructPosition() = 0;
    virtual int getNumStructNodes() = 0;
    double totalForce[3];
};

#endif
