#ifndef _PHYSBAMINTERSECT_H_
#define _PHYSBAMINTERSECT_H_

#include <string>

#include "LevelSet/LevelSetStructure.h"
#include "LevelSet/IntersectionFactory.h"
#include "PHYSBAM_INTERFACE.h"
#include "Particles/SOLIDS_PARTICLE.h"
#include "Grids/TRIANGLE_MESH.h"
#include "Geometry/TRIANGULATED_SURFACE.h"
#include <Vector.h>

using std::pair;
using std::map;
using PhysBAM::PhysBAMInterface;
using PhysBAM::LIST_ARRAY;
using PhysBAM::PAIR;
using PhysBAM::VECTOR;
using PhysBAM::IntersectionResult;

class Vec3D;
class Communicator;
class PhysBAMIntersector;
class SubDomain;
class EdgeSet;
template<class Scalar, int dim> class SVec;

class DistPhysBAMIntersector : public DistLevelSetStructure {
  typedef pair<int, int> iipair;
  typedef pair<int, bool> ibpair;
  typedef pair<iipair, ibpair> EdgePair;

  protected:
  public: // For debugging
    int length_solids_particle_list, length_triangle_list;
    int (*triangle_list)[3];
    double xMin, xMax, yMin, yMax, zMin, zMax;

    Vec3D *solids_particle_list;
    Vec3D *solids_particle_list0;
    Vec3D *solids_particle_list_n;
    Vec3D *solids_particle_list_nPlus1;
    Vec3D *solidVel;

    Vec<Vec3D> *solidX;
    Vec<Vec3D> *solidXn;
    Vec<Vec3D> *solidX0;

    Vec3D *triNorms;
    double *triSize;
    Communicator *com;
    PhysBAMIntersector **intersector;

    PhysBAMInterface<double> *physInterface;
    Domain *domain;
    DistSVec<double,3> *X;
    double tolerance;
    double insidePointTol;

    // A point known to be inside of the closed structural surface.
    Vec3D insidePoint;

    void buildSolidNormals();
    void getBoundingBox();
  public:
    DistPhysBAMIntersector(double tol);
    void init(std::string structureFileName, std::string structureFileName);

    double getTolerance() const { return tolerance; }
    EdgePair makeEdgePair(int,int,int);
    bool checkTriangulatedSurface();
    void initializePhysBAM();

    void initialize(Domain *, DistSVec<double,3> &X);
    void updateStructure(Vec3D *Xs, Vec3D *Vs, int nNodes);
    void updatePhysBAMInterface(Vec3D *particles, int size);
    void recompute(double dtf, double dtfLeft, double dts);

    LevelSetStructure & operator()(int subNum) const;

    PhysBAMInterface<double> &getInterface() { return *physInterface; }
    const Vec3D &getSurfaceNorm(int i) const {return triNorms[i]; }
    const Vec3D getInsidePoint() const { return insidePoint; }

    Vec<Vec3D> &getStructPosition() { return *solidX; }
    Vec<Vec3D> &getStructPosition_0() { return *solidX0; }
    Vec<Vec3D> &getStructPosition_n() { return *solidXn; }
    int getNumStructNodes () { return length_solids_particle_list; }
};

class PhysBAMIntersector : public LevelSetStructure {
  public:
    static const int UNDECIDED = -1, INSIDE = 0, OUTSIDE = 1; //INSIDE: inside real fluid, OUTSIDE: ~~

  protected:
  public: // For debug

    int globIndex;
    int *locToGlobNodeMap;

    std::map<int,IntersectionResult<double> > CrossingEdgeRes;
    std::map<int,IntersectionResult<double> > ReverseCrossingEdgeRes;

    DistPhysBAMIntersector &distIntersector;
    Vec<int> status; //<! Whether a node is inside the fluid domain or not
    Vec<int> status0; //<! status at the previous time-step.
    EdgeSet &edges;
    int nFirstLayer;

    /** compute projection of any point on the plane of a triangle.
     * send back the signed distance and barycentric coordinates of the projection point. */
    void projection(Vec3D, int, double&, double&, double&);

    double isPointOnSurface(Vec3D pt, int N1, int N2, int N3);
    /** check the distance of apoint to a surface defined by a triangle. (used for debug only) */ 

  public:
    PhysBAMIntersector(SubDomain &, SVec<double, 3> &X, DistPhysBAMIntersector &);
    void reset(); //<! set status0=status and reset status and nFirstLayer.
    void getClosestTriangles(SVec<double,3> &X, SVec<double,3> &boxMin, SVec<double,3> &boxMax, Vec<int> &tId, Vec<double> &dist);
    /** Function to compute a signed distance and normal estimates for nodes that are next to the structure
     *
     * results are for the subdomain only */
    void computeFirstLayerNodeStatus(Vec<int> tId, Vec<double> dist);
    /** compute the status (inside / outside) of nodes close to interface. */
    void fixUntouchedSubDomain(SVec<double,3>&X);
    /** if a subdomain has no intersection, find if it's (entirely) inside or outside. */
    void finishNodeStatus(SubDomain& sub, SVec<double,3>&X);
    /** compute the status (inside / outside) for all the nodes. */
    void findIntersections(SVec<double,3>&X);
    /** find intersections for each edge that has nodes with different statuses */

    LevelSetResult
    getLevelSetDataAtEdgeCenter(double t, int ni, int nj);
    bool isActive(double t, int n) const;
    bool wasActive(double t, int n) const;
    bool edgeIntersectsStructure(double t, int ni, int nj) const;

};

#endif
