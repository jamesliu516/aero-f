#ifndef _LEVEL_SET_STRUCTURE_H_
#define _LEVEL_SET_STRUCTURE_H_

#include "Vector3D.h"

class Domain;
template<class Scalar, int dim>
class DistSVec;
template <class Scalar> class Vec;
template <class Scalar> class DistVec;

/** Structure used to return levelset information */
struct LevelSetResult {
   Vec3D gradPhi;
   Vec3D normVel;

   LevelSetResult(double gpx, double gpy, double gpz,
                  double nvx, double nvy, double nvz) :
		  gradPhi(gpx, gpy, gpz), normVel(nvx, nvy, nvz) {
		     gradPhi *= 1.0/gradPhi.norm();
		   }

};

// ------- test for load transfer ------------------------------
class PistonSurface1D {
  public:
    Vec3D nodes[4];
    int elems[2][3];
    Vec3D nodalForce[4];

    PistonSurface1D() { 
      nodes[0] = Vec3D(0.0, -0.25, -0.25);
      nodes[1] = Vec3D(0.0, 0.25, -0.25);
      nodes[2] = Vec3D(0.0, -0.25, 0.25);
      nodes[3] = Vec3D(0.0, 0.25, 0.25);
      elems[0][0] = 0; elems[0][1] = 1; elems[0][2] = 2;
      elems[1][0] = 1; elems[1][1] = 2; elems[1][2] = 3;
    }
    ~PistonSurface1D(){}

    void addForce(Vec3D X, Vec3D force) {
      for (int iElem=0; iElem<2; iElem++) {
        Vec3D xprod = (nodes[elems[iElem][1]]-nodes[elems[iElem][0]])^
                      (nodes[elems[iElem][2]]-nodes[elems[iElem][0]]);
        double area = 0.5*xprod.norm();
        X[0] = nodes[0][0];
        double coord[3];
        for (int iNode = 0; iNode<3; iNode++) {
          xprod = (nodes[elems[iElem][(iNode+1)%3]]-X)^
                  (nodes[elems[iElem][(iNode+2)%3]]-X);
          coord[iNode] = 0.5*xprod.norm()/area;
        }
        if (fabs(coord[0]+coord[1]+coord[2]-1.0)<1e-6) {
          nodalForce[elems[iElem][0]] += coord[0]*force;
          nodalForce[elems[iElem][1]] += coord[1]*force;
          nodalForce[elems[iElem][2]] += coord[2]*force;
          break;
        }
        if (iElem==1) fprintf(stderr,"ERROR in load transfer. Abort...\n");
      }
    }
    void clearForce() {for (int i=0; i<4; i++) nodalForce[i] = 0.0;} 
    Vec3D getTotalForce() {return (nodalForce[0]+nodalForce[1]+nodalForce[2]+nodalForce[3]);} 
}; 
// --------------------------------------------------------------

/** Abstract class for finding levelset information */
class LevelSetStructure {
  public:
    /** returns the normal and normal velocity at intersection between edge ni, nj and structure
     *
     * If ni, nj is not an edge of the fluid mesh, result is undefined.
     * */
    virtual LevelSetResult
       getLevelSetDataAtEdgeCenter(double t, int ni, int nj) = 0;
    virtual bool isActive(double t, int n) = 0; //!< Whether this node is active or ghost.
    virtual bool edgeIntersectsStructure(double t, int ni, int nj) const = 0; //!< whether an edge between i and j intersects the structure

    /** creates an array of values which are positive inside the fluid and negative outside. */
    virtual void computePhi(Vec<double> &phi);

    Vec3D totalForce;
    PistonSurface1D solidSurface;
    void sendNodalForceToStruct(Vec3D X, Vec3D force) {solidSurface.addForce(X,force);}
};

class DistLevelSetStructure {
  protected:
    int numLocSub;
    DistVec<double> *pseudoPhi;
  public:
    virtual ~DistLevelSetStructure() {}

    virtual void initialize(Domain *, DistSVec<double,3> &X) = 0;
    virtual LevelSetStructure & operator()(int subNum) const = 0;
    virtual void clearTotalForce();
    virtual Vec3D getTotalForce(const double pref);

    virtual DistVec<double> &getPhi();

    double totalForce[3];
};

#endif
