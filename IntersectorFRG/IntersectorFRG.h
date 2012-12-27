#ifndef _INTERSECTORFRG_H_
#define _INTERSECTORFRG_H_

#include <list>
#include <map>
#include <set>
#include <string>

#include <Edge.h>
#include <LevelSet/LevelSetStructure.h>
#include <Vector.h>

#include "PHYSBAM_INTERFACE.h"
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>

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
class IntersectorFRG;
class SubDomain;
class EdgeSet;
class Timer;
class IoData;
template<class Scalar, int dim> class SVec;

class DistIntersectorFRG : public DistLevelSetStructure {

  friend class IntersectorFRG;

  typedef pair<int, int> iipair;
  typedef pair<int, bool> ibpair;
  typedef pair<iipair, ibpair> EdgePair;
  using DistLevelSetStructure::status;
  using DistLevelSetStructure::distance;
  using DistLevelSetStructure::is_swept;
  using DistLevelSetStructure::is_active;
  using DistLevelSetStructure::edge_intersects;

  protected:
    int numStNodes, numStElems;
    DistSVec<double,3> *boxMax, *boxMin; //fluid node bounding boxes
    DistVec<int> *tId;

    bool twoPhase; //including fluid-shell-fluid and fluid-solid
    
    // struct node coords
    Vec3D *Xs;
    Vec3D *Xs0;
    Vec3D *Xs_n;
    Vec3D *Xs_np1;
    Vec<Vec3D> *solidX;   //pointer to Xs
    Vec<Vec3D> *solidXn;  //pointer to Xs_n
    Vec<Vec3D> *solidX0;  //pointer to Xs0

    int (*stElem)[3];     //structural elements (element to node connectivity)
    set<int> *node2node;  //structural node to node connectivity
    set<int> *node2elem;  //structural node to element (triangle) connectivity
    
    Vec3D *Xsdot; //velocity

    DistVec<int> *status0;  //previous node status

    Vec3D *triNorms;
    Vec3D *nodalNormal; //memory allocated only if interpolatedNormal == true

    DistSVec<double,3> *X; //pointer to fluid node coords

    // parameters from input
    bool interpolatedNormal;

    // sophisticated stuff...
    Communicator *com;
    IntersectorFRG **intersector;
    PhysBAMInterface<double> *globPhysInterface;
    Domain *domain;

    DistVec<bool> *poly; //true if a node lies in n>2 subdomains (not used!)
    void findPoly(); //not used!
    void buildConnectivity();
    void buildSolidNormals();
    void expandScope();
    void findInAndOut();
    void finishStatusByPoints(IoData &iod, DistVec<int> *point_based_id = 0);
    void finalizeStatus(); //contains a local communication
    void updateStructCoords(double, double);

  public: //TODO: a lot of them can be moved to "protected".
    DistIntersectorFRG(IoData &iod, Communicator *comm, int nNodes = 0, double *xyz = 0, int nElems = 0, int (*abc)[3] = 0);
    virtual ~DistIntersectorFRG();

    void init(char *meshfile, char *restartfile, double XScale);
    void init(int nNodes, double *xyz, int nElems, int (*abc)[3], char *restartSolidSurface);

    EdgePair makeEdgePair(int,int,int);
    bool checkTriangulatedSurface();
    void initializePhysBAM();

    void initialize(Domain *, DistSVec<double,3> &X, IoData &iod, DistVec<int> *point_based_id = 0);
    void updateStructure(double* xs, double *Vs, int nNodes, int(*abc)[3]=0);
    void updatePhysBAMInterface();
    int recompute(double dtf, double dtfLeft, double dts, bool findStatus, bool retry = false);

    LevelSetStructure & operator()(int subNum) const;

    PhysBAMInterface<double> &getInterface() { return *globPhysInterface; }
    const Vec3D &getSurfaceNorm(int i) const {return triNorms[i]; }
    const Vec3D &getNodalNorm(int i) const {if (!nodalNormal) {fprintf(stderr,"ERROR: nodal normal not initialized!\n");exit(-1);} return nodalNormal[i];}

    Vec<Vec3D> &getStructPosition() { return *solidX; }
    Vec<Vec3D> &getStructPosition_0() { return *solidX0; }
    Vec<Vec3D> &getStructPosition_n() { return *solidXn; }
    DistVec<ClosestPoint> * getClosestPointsPointer() {return NULL;}
    DistVec<ClosestPoint> & getClosestPoints() {
      fprintf(stderr,"ERROR: closest point not stored in IntersectorFRG.\n");exit(-1);
      DistVec<ClosestPoint> *toto = new DistVec<ClosestPoint>(status->info()); return *toto;}
    void setStatus(DistVec<int> nodeTag) { *status = nodeTag; } //for reset after failSafe

    int getNumStructNodes () { return numStNodes; }
    int getNumStructElems () { return numStElems; }
};

class IntersectorFRG : public LevelSetStructure {

  friend class DistIntersectorFRG;

  using LevelSetStructure::status;
  using LevelSetStructure::distance;
  using LevelSetStructure::is_swept;
  using LevelSetStructure::is_active;
  using LevelSetStructure::edge_intersects;

  public:
    static const int OUTSIDE = -2, UNDECIDED = -1, INSIDE = 0; //INSIDE: inside real fluid, OUTSIDE: ~~
    static int OUTSIDECOLOR;

  protected:
    int globIndex;
    int *locToGlobNodeMap;

    SubDomain *subD;
    EdgeSet &edges;
    DistIntersectorFRG &distIntersector;

    set<int> *package;
    map<int,int> sub2pack;
    Connectivity &nodeToSubD;

    set<int> scope;
    int *iscope; 
    map<int,int> n2p; //node Id (index from 0) -> particle Id (index from 0 NOT 1!!!)
    list<int> particle; //particles in the scope

    map<int,IntersectionResult<double> > CrossingEdgeRes;
    map<int,IntersectionResult<double> > ReverseCrossingEdgeRes;
    Vec<int> &status0; //<! status at the previous time-step.

    PhysBAMInterface<double> *physInterface;
    bool testIsActive(double t, int n) const                         {return (status[n]>=0 && status[n]!=OUTSIDECOLOR);}

    void reset(const bool retry); //<! set status0=status and reset status and nFirstLayer.
    void rebuildPhysBAMInterface(Vec3D *Xs, int nsNodes, int (*sElem)[3], int nsElem);
    void getClosestTriangles(SVec<double,3> &X, SVec<double,3> &boxMin, SVec<double,3> &boxMax, Vec<int> &tId, Vec<double> &dist, bool useScope);
    void computeFirstLayerNodeStatus(Vec<int> tId, Vec<double> dist);
    bool finishStatusByHistory(SubDomain& sub);
    int findIntersections(SVec<double,3>&X, bool useScope);

    int buildScopeTopology(int (*sElem)[3], int nsElems);
    void projection(Vec3D, int, double&, double&, double&);
    void floodFill(SubDomain& sub, int& nUndecided);
    void noCheckFloodFill(SubDomain& sub, int& nUndecided);
    int findNewSeedsAfterMerging(Vec<int>& status_temp, Vec<bool>& poly, int& nUndecided); //not used!
    int findNewSeedsAfterMerging(SVec<int,2>& status_and_weight, int& nUndecided);
    int findSeedsByPoints(SubDomain& sub, SVec<double,3>& X, list< pair<Vec3D,int> > P, int& nUndecided);
    void addToPackage(int node, int trID);

    // for debug 
    int nFirstLayer;
    double isPointOnSurface(Vec3D pt, int N1, int N2, int N3);
    void printFirstLayer(SubDomain& sub, SVec<double,3>& X, int TYPE = 1);

  public:
    IntersectorFRG(SubDomain &, SVec<double, 3> &X, Vec<int> &status0, DistIntersectorFRG &);
    ~IntersectorFRG();

    int numOfFluids()                                            {return distIntersector.numOfFluids();}
    bool isNearInterface(double t, int n) const                  {return false;}
    bool withCracking() const                                    {return false;}

    LevelSetResult getLevelSetDataAtEdgeCenter(double t, int l, bool i_less_j);
    void findNodesNearInterface(SVec<double, 3>&, SVec<double, 3>&, SVec<double, 3>&) {/* pure virtual in LevelSet */}

};

#endif
