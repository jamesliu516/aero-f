#ifndef _PHYSBAMINTERSECT_H_
#define _PHYSBAMINTERSECT_H_

#include <string>
#include <list>

#include "LevelSet/LevelSetStructure.h"
#include "LevelSet/IntersectionFactory.h"
#include "PHYSBAM_INTERFACE.h"
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Vector.h>

using std::pair;
using std::map;
using std::list;
using PhysBAM::PhysBAMInterface;
using PhysBAM::ARRAY;
using PhysBAM::PAIR;
using PhysBAM::VECTOR;
using PhysBAM::IntersectionResult;

class Vec3D;
class Communicator;
class PhysBAMIntersector;
class SubDomain;
class EdgeSet;
class Timer;
class IoData;
template<class Scalar, int dim> class SVec;

class DistPhysBAMIntersector : public DistLevelSetStructure {

  friend class PhysBAMIntersector;

  typedef pair<int, int> iipair;
  typedef pair<int, bool> ibpair;
  typedef pair<iipair, ibpair> EdgePair;

  protected:
    int numStNodes, numStElems;
    double xMin, xMax, yMin, yMax, zMin, zMax; //a bounding box over the struct body
    DistSVec<double,3> *boxMax, *boxMin; //fluid node bounding boxes

    // struct node coords
    Vec3D *Xs;
    Vec3D *Xs0;
    Vec3D *Xs_n;
    Vec3D *Xs_np1;
    Vec<Vec3D> *solidX;   //pointer to Xs
    Vec<Vec3D> *solidXn;  //pointer to Xs_n
    Vec<Vec3D> *solidX0;  //pointer to Xs0

    int (*stElem)[3]; //structure elements (topology)

    Vec3D *Xsdot; //velocity

    DistVec<int> *status;   //node status
    DistVec<int> *status0;  //previous node status

    double *triSize;
    Vec3D *triNorms;  
    Vec3D *nodalNormal; //memory allocated only if interpolatedNormal == true

    DistSVec<double,3> *X; //pointer to fluid node coords

    // parameters from input
    double tolerance;
    bool interpolatedNormal;

    // sophisticated stuff...
    Communicator *com;
    PhysBAMIntersector **intersector;
    PhysBAMInterface<double> *physInterface;
    Domain *domain;

    DistVec<bool> *poly; //true if a node lies in n>2 subdomains
    void findPoly();
    void buildSolidNormals();
    void getBoundingBox();
    void findInAndOut();
    void finishStatusByPoints(IoData &iod);

  public: //TODO: a lot of them can be moved to "protected".
    DistPhysBAMIntersector(double tol);
    void init(std::string structureFileName, std::string structureFileName);
    ~DistPhysBAMIntersector();

    double getTolerance() const { return tolerance; }
    EdgePair makeEdgePair(int,int,int);
    bool checkTriangulatedSurface();
    void initializePhysBAM();

    void initialize(Domain *, DistSVec<double,3> &X, IoData &iod, bool interpNormal);
    void updateStructure(Vec3D *xs, Vec3D *Vs, int nNodes);
    void updatePhysBAMInterface(Vec3D *particles, int size);
    void recompute(double dtf, double dtfLeft, double dts);

    LevelSetStructure & operator()(int subNum) const;

    PhysBAMInterface<double> &getInterface() { return *physInterface; }
    const Vec3D &getSurfaceNorm(int i) const {return triNorms[i]; }
    const Vec3D &getNodalNorm(int i) const {if (!nodalNormal) {fprintf(stderr,"ERROR: nodal normal not initialized!\n");exit(-1);} return nodalNormal[i];}

    Vec<Vec3D> &getStructPosition() { return *solidX; }
    Vec<Vec3D> &getStructPosition_0() { return *solidX0; }
    Vec<Vec3D> &getStructPosition_n() { return *solidXn; }
    DistVec<int> &getStatus() { return *status; }
    int getNumStructNodes () { return numStNodes; }
    int getNumStructElems () { return numStElems; }
};

class PhysBAMIntersector : public LevelSetStructure {
  public:
    static const int UNDECIDED = -1, INSIDE = 0, OUTSIDE = 1; //INSIDE: inside real fluid, OUTSIDE: ~~

  protected:

    int globIndex;
    int *locToGlobNodeMap;

    std::map<int,IntersectionResult<double> > CrossingEdgeRes;
    std::map<int,IntersectionResult<double> > ReverseCrossingEdgeRes;

    DistPhysBAMIntersector &distIntersector;
    Vec<int> &status; //<! Whether a node is inside the fluid domain or not
    Vec<int> &status0; //<! status at the previous time-step.
    EdgeSet &edges;
    int nFirstLayer;

    /** compute projection of any point on the plane of a triangle.
     * send back the signed distance and barycentric coordinates of the projection point. */
    void projection(Vec3D, int, double&, double&, double&);

  public:
    PhysBAMIntersector(SubDomain &, SVec<double, 3> &X, Vec<int> &status, Vec<int> &status0, DistPhysBAMIntersector &);
    int numOfFluids() {return distIntersector.numOfFluids();}
    void reset(); //<! set status0=status and reset status and nFirstLayer.
    void getClosestTriangles(SVec<double,3> &X, SVec<double,3> &boxMin, SVec<double,3> &boxMax, Vec<int> &tId, Vec<double> &dist);
    /** Function to compute a signed distance and normal estimates for nodes that are next to the structure
     *
     * results are for the subdomain only */
    void computeFirstLayerNodeStatus(Vec<int> tId, Vec<double> dist);
    /** compute the status (inside / outside) for all the nodes. */
    void findIntersections(SVec<double,3>&X);
    /** find intersections for each edge that has nodes with different statuses */
    double isPointOnSurface(Vec3D pt, int N1, int N2, int N3);
    /** check the distance of apoint to a surface defined by a triangle. (used for debug only) */ 
    void floodFill(SubDomain& sub, int& nUndecided);
    void noCheckFloodFill(SubDomain& sub, int& nUndecided);
    int findNewSeedsAfterMerging(Vec<int>& status_temp, Vec<bool>& poly, int& nUndecided);
    int findSeedsByPoints(SubDomain& sub, SVec<double,3>& X, list< pair<Vec3D,int> > P, int& nUndecided);
    void finishStatusByHistory(SubDomain& sub);
    void printFirstLayer(SubDomain& sub, SVec<double,3>& X, int TYPE = 1);

    LevelSetResult
    getLevelSetDataAtEdgeCenter(double t, int ni, int nj);
    bool isActive(double t, int n, int phase = 0) const;
    bool wasActive(double t, int n, int phase = 0) const;
    bool edgeIntersectsStructure(double t, int ni, int nj) const;
    void findNodesNearInterface(SVec<double, 3>&, SVec<double, 3>&, SVec<double, 3>&) {}

};

#endif
