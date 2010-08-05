#ifndef _PHYSBAMINTERSECT_H_
#define _PHYSBAMINTERSECT_H_

#include <map>
#include <string>
#include <list>
#include <set>

#include "../LevelSet/LevelSetStructure.h"
#include "../Vector.h"
#include "PHYSBAM_INTERFACE.h"
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

using std::pair;
using std::map;
using std::list;
using std::set;
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
    DistSVec<double,3> *boxMax, *boxMin; //fluid node bounding boxes
    DistVec<double> *distance;
    DistVec<int> *tId;

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

    Vec3D *triNorms;
    Vec3D *nodalNormal; //memory allocated only if interpolatedNormal == true

    DistSVec<double,3> *X; //pointer to fluid node coords

    // parameters from input
    bool interpolatedNormal;

    // sophisticated stuff...
    Communicator *com;
    PhysBAMIntersector **intersector;
    PhysBAMInterface<double> *globPhysInterface;
    Domain *domain;

    DistVec<bool> *poly; //true if a node lies in n>2 subdomains
    void findPoly();
    void buildSolidNormals();
    void expandScope();
    void findInAndOut();
    void finishStatusByPoints(IoData &iod);
    void updateStructCoords(double, double);

  public: //TODO: a lot of them can be moved to "protected".
    DistPhysBAMIntersector(IoData &iod, Communicator *comm, int nNodes = 0, double (*xyz)[3] = 0, int nElems = 0, int (*abc)[3] = 0);
    ~DistPhysBAMIntersector();

    void init(char *meshfile, char *restartfile);
    void init(int nNodes, double (*xyz)[3], int nElems, int (*abc)[3], char *restartSolidSurface);

    EdgePair makeEdgePair(int,int,int);
    bool checkTriangulatedSurface();
    void initializePhysBAM();

    void initialize(Domain *, DistSVec<double,3> &X, IoData &iod);
    void updateStructure(Vec3D *xs, Vec3D *Vs, int nNodes);
    void updatePhysBAMInterface();
    void recompute(double dtf, double dtfLeft, double dts);

    LevelSetStructure & operator()(int subNum) const;

    PhysBAMInterface<double> &getInterface() { return *globPhysInterface; }
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

  friend class DistPhysBAMIntersector;

  public:
    static const int UNDECIDED = -1, INSIDE = 0, OUTSIDE = 1; //INSIDE: inside real fluid, OUTSIDE: ~~

  protected:
    int globIndex;
    int *locToGlobNodeMap;

    SubDomain *subD;
    EdgeSet &edges;
    DistPhysBAMIntersector &distIntersector;

    set<int> *package;
    map<int,int> sub2pack;
    Connectivity &nodeToSubD;

    set<int> scope;
    int *iscope; 
    map<int,int> n2p; //node Id (index from 0) -> particle Id (index from 0 NOT 1!!!)
    list<int> particle; //particles in the scope

    map<int,IntersectionResult<double> > CrossingEdgeRes;
    map<int,IntersectionResult<double> > ReverseCrossingEdgeRes;
    Vec<int> &status; //<! Whether a node is inside the fluid domain or not
    Vec<int> &status0; //<! status at the previous time-step.

    PhysBAMInterface<double> *physInterface;
    PhysBAM::GEOMETRY_PARTICLES<PhysBAM::VECTOR<double,3> > *physbam_solids_particle;
    ARRAY<VECTOR<int,3> > *physbam_stElem;
    PhysBAM::TRIANGLE_MESH *physbam_triangle_mesh;
    PhysBAM::TRIANGULATED_SURFACE<double> *physbam_triangulated_surface;

    void reset(); //<! set status0=status and reset status and nFirstLayer.
    void rebuildPhysBAMInterface(Vec3D *Xs, int nsNodes, int (*sElem)[3], int nsElem);
    void getClosestTriangles(SVec<double,3> &X, SVec<double,3> &boxMin, SVec<double,3> &boxMax, Vec<int> &tId, Vec<double> &dist, bool useScope);
    void computeFirstLayerNodeStatus(Vec<int> tId, Vec<double> dist);
    void finishStatusByHistory(SubDomain& sub);
    void findIntersections(SVec<double,3>&X, bool useScope);

    int buildScopeTopology(int (*sElem)[3], int nsElems);
    void projection(Vec3D, int, double&, double&, double&);
    void floodFill(SubDomain& sub, int& nUndecided);
    void noCheckFloodFill(SubDomain& sub, int& nUndecided);
    int findNewSeedsAfterMerging(Vec<int>& status_temp, Vec<bool>& poly, int& nUndecided);
    int findSeedsByPoints(SubDomain& sub, SVec<double,3>& X, list< pair<Vec3D,int> > P, int& nUndecided);
    void addToPackage(int node, int trID);

    // for debug 
    int nFirstLayer;
    double isPointOnSurface(Vec3D pt, int N1, int N2, int N3);
    void printFirstLayer(SubDomain& sub, SVec<double,3>& X, int TYPE = 1);

  public:
    PhysBAMIntersector(SubDomain &, SVec<double, 3> &X, Vec<int> &status, Vec<int> &status0, DistPhysBAMIntersector &);
    ~PhysBAMIntersector();

    int numOfFluids()                                            {return distIntersector.numOfFluids();}
    bool isActive(double t, int n, int phase = 0) const          {return status[n]==phase;}
    bool wasActive(double t, int n, int phase = 0) const         {return status0[n]==phase;}
    bool edgeIntersectsStructure(double t, int ni, int nj) const {return status[ni]!=status[nj];}

    LevelSetResult getLevelSetDataAtEdgeCenter(double t, int ni, int nj);
    void findNodesNearInterface(SVec<double, 3>&, SVec<double, 3>&, SVec<double, 3>&) {/* pure virtual in LevelSet */}

};

#endif
