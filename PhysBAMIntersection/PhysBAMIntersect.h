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
  protected:
  public: // For debugging
    int length_solids_particle_list, length_triangle_list;
    int (*triangle_list)[3];
    double xMin, xMax, yMin, yMax, zMin, zMax;
    Vec3D *solids_particle_list;
    Vec<Vec3D> *solidX;
    Vec3D *triNorms;
    double *triSize;
    Communicator *com;
    PhysBAMIntersector **intersector;

    PhysBAMInterface<double> *physInterface;
    Domain *domain;
    DistSVec<double,3> *X;
    double tolerance;

    // A point known to be inside of the closed structural surface.
    Vec3D insidePoint;

    void buildSolidNormals();
    void getBoundingBox();
  public:
    DistPhysBAMIntersector(double tol);
    void init(std::string structureFileName);

    double getTolerance() const { return tolerance; }
    bool checkTriangulatedSurface();
    void initializePhysBAM();

    void initialize(Domain *, DistSVec<double,3> &X);
    LevelSetStructure & operator()(int subNum) const;

    PhysBAMInterface<double> &getInterface() { return *physInterface; }
    const Vec3D &getSurfaceNorm(int i) const {return triNorms[i]; }
    const Vec3D getInsidePoint() const { return insidePoint; }

    Vec<Vec3D> &getStructPosition() { return *solidX; }
    int getNumStructNodes () { return length_solids_particle_list; }
};

class PhysBAMIntersector : public LevelSetStructure {
  public:
    static const int UNDECIDED = -1, INSIDE = 0, OUTSIDE = 1;

  protected:
  public: // For debug

    int globIndex;
    int *locToGlobNodeMap;
    int *nodeMap;

    std::map<int,IntersectionResult<double> > CrossingEdgeRes;
    std::map<int,IntersectionResult<double> > ReverseCrossingEdgeRes;

    DistPhysBAMIntersector &distIntersector;
    Vec<int> status; //<! Whether a node is inside the fluid domain or not
    EdgeSet &edges;
    int nFirstLayer;

//    void updatePhi(int p, int q, IntersectionResult<double> &res, SVec<double, 3> &X, Vec<double> &phi,
//                   SVec<double, 3> &normApprox, Vec<double> &weightSum);

    /** compute projection of any point on the plane of a triangle.
     * send back the signed distance and barycentric coordinates of the projection point. */
    void projection(Vec3D, int, double&, double&, double&);

  public:
    PhysBAMIntersector(SubDomain &, SVec<double, 3> &X, DistPhysBAMIntersector &);
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

//    void computeLocalPseudoPhi(SVec<double,3> &X, SVec<double,3> &n, Vec<double> &lWeight);
    /** complete the computation of the signed distance after communication with neighboring subdomains
     * has been done
     *
     * status is created and the locNorm contains our estimate of the normal for nodes near the structure*/
//    void finishPseudoPhi(SubDomain &sub, SVec<double,3> &X, SVec<double,3> &n, Vec<double> &lWeight);

    LevelSetResult
    getLevelSetDataAtEdgeCenter(double t, int ni, int nj);
    bool isActive(double t, int n);
    bool edgeIntersectsStructure(double t, int ni, int nj) const;

};

#endif
