#ifndef _LEVEL_SET_STRUCTURE_H_
#define _LEVEL_SET_STRUCTURE_H_
 
#include "Vector3D.h"

class Domain;
class IoData;
template<class Scalar, int dim>
class DistSVec;
template <class Scalar> class Vec;
template <class Scalar, int dim> class SVec;
template <class Scalar> class DistVec;

/** Structure used to return levelset information */
struct LevelSetResult {
   double alpha;
   double xi[3];
   int trNodes[3];
   Vec3D gradPhi;
   Vec3D normVel; //NOTE: this is the velocity, NOT normal velocity.

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
    virtual bool isActive(double t, int n, int phase = 0) const = 0; //!< Whether this node is active or ghost.
    virtual bool wasActive(double t, int n, int phase = 0) const = 0; //!< Whether this node was active or ghost in the last iteration.
    virtual bool edgeIntersectsStructure(double t, int ni, int nj) const = 0; //!< whether an edge between i and j intersects the structure

    /** creates an array of values which are positive inside the fluid and negative outside. */
    virtual void computePhi(Vec<double> &phi);

    virtual double isPointOnSurface(Vec3D, int, int, int) = 0;

    virtual int numOfFluids() = 0; 

    virtual void findNodesNearInterface(SVec<double,3>&, SVec<double,3>&, SVec<double,3>&) = 0;
};

class DistLevelSetStructure {
  protected:
    int numLocSub;
    DistVec<double> *pseudoPhi;
    int numFluid;

  public:
    virtual ~DistLevelSetStructure() {}

    int numOfFluids() {return numFluid;}
    void setNumOfFluids(int nf) {numFluid = nf;}
    virtual void initialize(Domain *, DistSVec<double,3> &X, IoData &iod) = 0;
    virtual LevelSetStructure & operator()(int subNum) const = 0;

    virtual DistVec<double> &getPhi();
    virtual DistVec<int> &getStatus() {fprintf(stderr,"Not implemented yet!\n");} //TODO: to be fixed 

    virtual void updateStructure(Vec3D *Xs, Vec3D *Vs, int nNodes) = 0;
    virtual void recompute(double dtf, double dtfLeft, double dts) = 0;
    virtual Vec<Vec3D> &getStructPosition() = 0;
    virtual Vec<Vec3D> &getStructPosition_0() = 0;
    virtual Vec<Vec3D> &getStructPosition_n() = 0;
    virtual int getNumStructNodes() = 0;
};

#endif
