#ifndef _LEVEL_SET_STRUCTURE_H_
#define _LEVEL_SET_STRUCTURE_H_
 
#include "Vector3D.h"
#include <DistVector.h>
#include <Vector.h>

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
	   iterator & operator++() { xip++; nodep++; return (*this); }
	   bool operator !=(const iterator &it) const {
		   return xip != it.xip || nodep != it.nodep;
	   }
	   double Ni() const { return *xip; }
	   int nodeNum() const { return *nodep; }
   };

   iterator begin() { return iterator(xi, trNodes); }
   iterator end() { return iterator(xi+3, trNodes+3); }
};

/** Storing the closest point on the interface */
struct ClosestPoint {
  int mode; //-2: unknown, -1: far from interface, 0: face, 1: edge, 2: vertex.
  double dist; //this is the unsigned distance. always >= 0.
  int tracker[2]; // for mode=0: tracker[0] = tria Id; for mode=1: the two vertices; for mode=2: the vertex
  double xi1,xi2; // local coordinates.
  ClosestPoint() {mode=-2;}
  bool known() {return mode!=-2;}
  bool nearInterface() {return mode!=-1;}
};

/** Abstract class for finding levelset information */
class LevelSetStructure {
  protected:
    Vec<int> &status;
    Vec<double> &distance;
    Vec<bool> &is_swept;
    Vec<bool> &is_active;
    Vec<bool> &is_occluded;
    Vec<bool> &edge_intersects;
  public:
    LevelSetStructure(Vec<int>& status,Vec<double>& distance,Vec<bool>& is_swept,Vec<bool>& is_active,Vec<bool>& is_occluded,Vec<bool>& edge_intersects)
        : status(status),distance(distance),is_swept(is_swept),is_active(is_active),is_occluded(is_occluded),edge_intersects(edge_intersects)
    {}
    virtual ~LevelSetStructure()
    {}

    /** returns the normal and normal velocity at intersection between edge ni, nj and structure
     *
     * If ni, nj is not an edge of the fluid mesh, result is undefined.
     * */
    virtual LevelSetResult
       getLevelSetDataAtEdgeCenter(double t, int l, bool i_less_j) = 0;
    virtual bool withCracking() const = 0;
    virtual bool isNearInterface(double t, int n) const = 0;

    void forceOccluded(double t, int n) const                { is_swept[n] = true; is_occluded[n] = true; }
    int fluidModel(double t, int n) const                 { return status[n]; }
    double distToInterface(double t, int n) const         { return distance[n]; } 
    bool isSwept(double t, int n) const                   { return is_swept[n]; }
    bool isActive(double t, int n) const                  { return is_active[n]; }
    bool isOccluded(double t, int n) const                { return is_occluded[n]; }
    bool edgeIntersectsStructure(double t, int eij) const { return edge_intersects[eij]; }
    void computeSwept(Vec<int> &swept){
        for(int i = 0; i < swept.size(); ++i)
            swept[i] = is_swept[i] ? 1 : 0;
    }

    Vec<int> & getStatus() { return status; }

    virtual class CrackingSurface* 
      getCrackingSurface() { return NULL; }

    virtual double isPointOnSurface(Vec3D, int, int, int) = 0;

    virtual int numOfFluids() = 0; 

    virtual void findNodesNearInterface(SVec<double,3>&, SVec<double,3>&, SVec<double,3>&) = 0;
};

class DistLevelSetStructure {
  protected:
    DistVec<int> *status;
    DistVec<double> *distance;
    DistVec<bool> *is_swept;
    DistVec<bool> *is_active;
    DistVec<bool> *is_occluded;
    DistVec<bool> *edge_intersects;
  protected:
    int numLocSub;
    int numFluid;

  public:
    DistLevelSetStructure()
        : status(0), distance(0), is_swept(0), is_active(0), is_occluded(0), edge_intersects(0)
    {}
    virtual ~DistLevelSetStructure()
    {delete status;delete distance;delete is_swept;delete is_active;delete is_occluded;delete edge_intersects;}

    int numOfFluids() {return numFluid;}
    void setNumOfFluids(int nf) {numFluid = nf;}
    virtual void initialize(Domain *, DistSVec<double,3> &X, DistSVec<double,3> &Xn, IoData &iod, DistVec<int> *point_based_id = 0, DistVec<int>* oldStatus = 0) = 0;
    virtual LevelSetStructure & operator()(int subNum) const = 0;

    void getSwept(DistVec<int>& swept){
        for (int iSub=0; iSub<numLocSub; iSub++)
            (*this)(iSub).computeSwept((swept)(iSub));
    }

    DistVec<int> & getStatus()            const { return *status; }
    DistVec<double> & getDistance()       const { return *distance; }
    DistVec<bool> & getIsSwept()          const { return *is_swept; }
    DistVec<bool> & getIsActive()         const { return *is_active; }
    DistVec<bool> & getIsOccluded()       const { return *is_occluded; }
    DistVec<bool> & getIntersectedEdges() const { return *edge_intersects; }

    virtual DistVec<ClosestPoint> &getClosestPoints() = 0;
    virtual DistVec<ClosestPoint> *getClosestPointsPointer() = 0;
    virtual void setStatus(DistVec<int> nodeTag) = 0;                                

    virtual void updateStructure(double *Xs, double *Vs, int nNodes, int(*abc)[3]=0) = 0;
    virtual int recompute(double dtf, double dtfLeft, double dts, bool findStatus, bool retry = false) = 0;
    virtual Vec<Vec3D> &getStructPosition() = 0;
    virtual Vec<Vec3D> &getStructPosition_0() = 0;
    virtual Vec<Vec3D> &getStructPosition_n() = 0;
    virtual Vec<Vec3D> &getStructPosition_np1() = 0;
    virtual int getNumStructNodes() = 0;
    virtual int getNumStructElems() = 0;
    virtual int (*getStructElems())[3] = 0;
    virtual int getSurfaceID(int) = 0;
};

#endif
