#ifndef _INTERSECTORPHYSBAM_H_
#define _INTERSECTORPHYSBAM_H_

#include <string>
#include <list>
#include <Vector.h>

#include "../LevelSet/LevelSetStructure.h"

#include "PHYSBAM_INTERFACE.h"
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <set>
#include <map>

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
class IntersectorPhysBAM;
class FloodFill;
class SubDomain;
class EdgeSet;
class Timer;
class IoData;
class CrackingSurface;
template<class Scalar, int dim> class SVec;

class DistIntersectorPhysBAM : public DistLevelSetStructure {

  friend class IntersectorPhysBAM;

  typedef pair<int, int> iipair;
  typedef pair<int, bool> ibpair;
  typedef pair<iipair, ibpair> EdgePair;

  protected:
    int numStNodes, numStElems;
    int totStNodes, totStElems;
    double xMin, xMax, yMin, yMax, zMin, zMax; //a bounding box over the struct body
    DistSVec<double,3> *boxMax, *boxMin; //fluid node bounding boxes

    FloodFill* floodFill;

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

    CrackingSurface *cracking; //only a pointer.

  public:
    DistVec<int> *status;   //node status
    DistVec<int> *status0;  //previous node status
    DistVec<bool> *occluded_node; // identifies an occluded node
    DistVec<bool> *swept_node;    // identifies an invalidated node

    double *triSize;
    Vec3D *triNorms;
    Vec3D *nodalNormal; //memory allocated only if interpolatedNormal == true

    DistSVec<double,3> *X; //pointer to fluid node coords

    // parameters from input
    bool interpolatedNormal;

    // sophisticated stuff...
    Communicator *com;
    IntersectorPhysBAM **intersector;
    PhysBAMInterface<double> *physInterface;
    Domain *domain;

    void buildSolidNormals();
    void getBoundingBox();

    void findActiveNodesUsingFloodFill(const DistVec<bool>& tId,const list<pair<Vec3D,int> >& points);
    void findActiveNodes(const DistVec<bool>& tId);

  public: //TODO: a lot of them can be moved to "protected".
    DistIntersectorPhysBAM(IoData &iod, Communicator *comm, int nNodes = 0, double *xyz = 0, int nElems = 0, int (*abc)[3] = 0, CrackingSurface *cs = 0);
    ~DistIntersectorPhysBAM();

    void init(char *meshfile, char *restartfile, double XScale);
    void init(int nNodes, double *xyz, int nElems, int (*abc)[3], char *restartSolidSurface);

    EdgePair makeEdgePair(int,int,int);
    bool checkTriangulatedSurface();
    void initializePhysBAM();

    void initialize(Domain *, DistSVec<double,3> &X, IoData &iod, DistVec<int>* point_based_id = 0);
    void updateStructure(double *xs, double *Vs, int nNodes, int (*abc)[3]=0);
    void updateCracking(int (*abc)[3]);
    void expandScope();
    void updatePhysBAMInterface(Vec3D *particles, int size,const DistSVec<double,3>& fluid_nodes,const bool fill_scope=false);
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

class IntersectorPhysBAM : public LevelSetStructure {
  friend class DistIntersectorPhysBAM;

  public:
    static const int OUTSIDE = -2, UNDECIDED = -1, INSIDE = 0; //INSIDE: inside real fluid, OUTSIDE: not a fluid
    static int OUTSIDECOLOR;

    Vec<bool> edgeIntersections;

    int locIndex,globIndex;
    int *locToGlobNodeMap;

    std::map<int,IntersectionResult<double> > CrossingEdgeRes;
    std::map<int,IntersectionResult<double> > ReverseCrossingEdgeRes;

    SubDomain &subD;
    EdgeSet &edges;
    DistIntersectorPhysBAM &distIntersector;

    std::set<int> *package;
    map<int,int> sub2pack;
    Connectivity &nodeToSubD;

    Vec<int> &status; //<! Whether a node is inside the fluid domain or not
    Vec<int> &status0; //<! status at the previous time-step.
    Vec<bool> &occluded_node; //<! Whether a node is occluded by the solid surface.
    Vec<bool> &swept_node; //<! Whether a node is swept by the solid surface, over the present time-step.
    int nFirstLayer;
    ARRAY<int> reverse_mapping,forward_mapping;
    ARRAY<VECTOR<double,3> > xyz;

    int hasCloseTriangle(SVec<double,3>& X,SVec<double,3> &boxMin, SVec<double,3> &boxMax, Vec<bool> &tId);
    int findIntersections(SVec<double,3>& X,Vec<bool>& tId,Communicator&);
    int computeSweptNodes(SVec<double,3>& X, Vec<bool>& tId,Communicator&);

  public:
    IntersectorPhysBAM(SubDomain &, SVec<double, 3> &X, Vec<int> &status, Vec<int> &status0, Vec<bool>& occluded_node,Vec<bool>& swept_node, DistIntersectorPhysBAM &);
    ~IntersectorPhysBAM();
    int numOfFluids() {return distIntersector.numOfFluids();}

    void reset(); //<! set status0=status and reset status and nFirstLayer.
    /** find intersections for each edge that has nodes with different statuses */
    double isPointOnSurface(Vec3D pt, int N1, int N2, int N3);
    /** check the distance of apoint to a surface defined by a triangle. (used for debug only) */ 
    void printFirstLayer(SubDomain& sub, SVec<double,3>& X, int TYPE = 1);

    LevelSetResult
    getLevelSetDataAtEdgeCenter(double t, int ni, int nj);
    bool edgeIntersectsStructure(double t, int ni, int nj) const;
    bool edgeIntersectsStructure(double t, int edge_num) const;
    void findNodesNearInterface(SVec<double, 3>&, SVec<double, 3>&, SVec<double, 3>&) {}

    bool isActive(double t, int n) const {return (status[n] >= 0 && status[n]!=OUTSIDECOLOR);}
    bool isOccluded(double t, int n) const {return occluded_node[n];}
    bool isSwept(double t, int n) const {return swept_node[n];}
    int fluidModel(double t, int n) const {return status[n];}

  private:
    void addToPackage(const int i,const int candidate);
};

#endif
