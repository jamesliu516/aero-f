#ifndef _MULTIGRID_LEVEL_SET_STRUCTURE_H_
#define _MULTIGRID_LEVEL_SET_STRUCTURE_H_
 
#include "LevelSetStructure.h"

#include "Domain.h"

#include "MultiGridLevel.h"

/** Abstract class for finding levelset information */
class MultiGridLevelSetStructure : public LevelSetStructure {
  protected:
    Vec<int> &status;
    Vec<double> &distance;
    Vec<bool> &is_swept;
    Vec<bool> &is_active;
    Vec<bool> &is_occluded;
    Vec<bool> &edge_intersects;

    class DistMultiGridLevelSetStructure& distLSS;

    LevelSetStructure* parent;

    SubDomain &subD;
    EdgeSet &edges;

    int mySub;

    MultiGridLevel<double>* myLevel;

  public:
    MultiGridLevelSetStructure(class DistMultiGridLevelSetStructure& lss,
			       SubDomain& sub,
			       Vec<int>& status,Vec<double>& distance,Vec<bool>& is_swept,
			       Vec<bool>& is_active,Vec<bool>& is_occluded,
			       Vec<bool>& edge_intersects,
			       LevelSetStructure* parent,
			       int mySub,MultiGridLevel<double>* myLevel);

    ~MultiGridLevelSetStructure()
    {}

    /** returns the normal and normal velocity at intersection between edge ni, nj and structure
     *
     * If ni, nj is not an edge of the fluid mesh, result is undefined.
     * */
    LevelSetResult
       getLevelSetDataAtEdgeCenter(double t, int l, bool i_less_j);
    bool withCracking() const;
    bool isNearInterface(double t, int n) const;

    double isPointOnSurface(Vec3D, int, int, int);

    int numOfFluids();

    void findNodesNearInterface(SVec<double,3>&, SVec<double,3>&, SVec<double,3>&);

    void recompute();

    void computeEdgeCrossing();
};

class DistMultiGridLevelSetStructure : public DistLevelSetStructure {

  friend class MultiGridLevelSetStructure;

  // LSS for the parent level
  DistLevelSetStructure* parent;

  Communicator *com;

  MultiGridLevel<double>* myLevel;

  int numLocSub;

  MultiGridLevelSetStructure** subLSS;

  Domain* domain;

  public:
    DistMultiGridLevelSetStructure(IoData &iod, Communicator *comm,
				   DistLevelSetStructure* parent,
				   MultiGridLevel<double>*);

    virtual ~DistMultiGridLevelSetStructure()
      {}

    void initialize(Domain *, DistSVec<double,3> &X, DistSVec<double,3> &Xn, IoData &iod, DistVec<int> *point_based_id = 0, DistVec<int>* oldStatus = 0);
    LevelSetStructure & operator()(int subNum) const;

    DistVec<ClosestPoint> &getClosestPoints();
    DistVec<ClosestPoint> *getClosestPointsPointer();
    void setStatus(DistVec<int> nodeTag) = 0;                                

    void updateStructure(double *Xs, double *Vs, int nNodes, int(*abc)[3]=0) {

      parent->updateStructure(Xs, Vs, nNodes, abc);
    }


    int recompute(double dtf, double dtfLeft, double dts, bool findStatus, bool retry = false);


    Vec<Vec3D> &getStructPosition()  { return parent->getStructPosition(); }
    Vec<Vec3D> &getStructPosition_0()  { return parent->getStructPosition_0(); }
    Vec<Vec3D> &getStructPosition_n() { return parent->getStructPosition_n(); }
    Vec<Vec3D> &getStructPosition_np1() { return parent->getStructPosition_np1(); }
    int getNumStructNodes() { return parent->getNumStructNodes(); }
    int getNumStructElems() { return parent->getNumStructElems(); }
    int (*getStructElems())[3]  { return parent->getStructElems(); }

    int getSurfaceID(int k) const { return parent->getSurfaceID(k); }
};

#endif
