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
    Vec<Vec3D> *solidX;
    Vec3D *triNorms;
    int (*triToTri)[3];
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
    EdgePair makeEdgePair(int,int,int);
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
    void projection(Vec3D, int, double&, double&, double&);
    /** compute projection of any point on the plane of a triangle.
     * send back the signed distance and barycentric coordinates of the projection point. */

  public: // For debug
    int globIndex;
    int *locToGlobNodeMap;
    int *nodeMap;
    std::map<int,IntersectionResult<double> > secondIntersection;
    
    DistPhysBAMIntersector &distIntersector;
    Vec<int> status; //<! Whether a node is inside the fluid domain or not
    Vec<double> &phi; //<! Pseudo phi value
    Vec<Vec3D> locNorm;
    EdgeSet &edges;
    LIST_ARRAY<PAIR<VECTOR<int,2>,IntersectionResult<double> > > edgeRes;
    int nIntersect;
    void updatePhi(int p, int q, IntersectionResult<double> &res, SVec<double, 3> &X, Vec<double> &phi,
                   SVec<double, 3> &normApprox, Vec<double> &weightSum);
    
    double checkPointOnSurface(Vec3D pt, int N1, int N2, int N3) {
      Vec<Vec3D> &solidX = distIntersector.getStructPosition();
      Vec3D X1 = solidX[N1];
      Vec3D X2 = solidX[N2];
      Vec3D X3 = solidX[N3];
 
      Vec3D normal = (X2-X1)^(X3-X1);
      normal /=  normal.norm();
     
      return fabs((pt-X1)*normal)/((pt-X1).norm());
    }

  public:
    PhysBAMIntersector(SubDomain &, SVec<double, 3> &X, Vec<double> &phi, DistPhysBAMIntersector &);
    /** Function to compute a signed distance and normal estimates for nodes that are next to the structure
     *
     * results are for the subdomain only */
    int closestTriangle(Vec3D, int, double&);
    /** find the closest triangle to a point, starting with a candidate. 
     * returns the triangle ID and the signed distance between the point and THE PLANE of the triangle.*/ 
    void computeLocalPseudoPhi(SVec<double,3> &X, SVec<double,3> &n, Vec<double> &lWeight);
    /** complete the computation of the signed distance after communication with neighboring subdomains
     * has been done
     *
     * status is created and the locNorm contains our estimate of the normal for nodes near the structure*/
    void finishPseudoPhi(SubDomain &sub, SVec<double,3> &X, SVec<double,3> &n, Vec<double> &lWeight);

    LevelSetResult
    getLevelSetDataAtEdgeCenter(double t, int ni, int nj);
    bool isActive(double t, int n);
    bool edgeIntersectsStructure(double t, int ni, int nj) const;


};

#endif
