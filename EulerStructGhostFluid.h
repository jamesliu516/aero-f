// This is an subdomain interface to communicate with PhysBAM for ghost solid simulation. 

#ifndef _EULER_STRUCT_GHOST_FLUID_H_
#define _EULER_STRUCT_GHOST_FLUID_H_


#include <Elem.h>
#include <SubDomain.h>
#include <Vector.h>
#include <Vector3D.h>
#include <DistVector.h>
#include "../PhysBAM/Projects/Charbel/AERO_INTERFACE_1.h"

class Elem;
class ElemSet;
class Edge;
class EdgeSet;
class Node;
class NodeSet;
class TriangulatedSurface;
class SubDomain;
template<class Scalar, int dim> class SVec;

class EulerStructGhostFluid {
  TriangulatedSurface* triaSurf;
  
  Vec<double> &philevel;
  SVec<double,3> &gradPhilevel;
  Vec3D* n_phi; 

  Vec<bool> &realNodeTag0;
  Vec<bool> &realNodeTag; 
 
  SubDomain* subDomain;

  PhysBAM::AERO_INTERFACE_1<double> *PhysBAM_Interface;

  int* nodeTag;
  int numChosenNodes;
  int* tempNodeList;
  int numChosenElems;
  int* tempElemList;  //stores the global index of chosen elems.
  int* chainMail;  //stores the index of tets covering the interface.
  PhysBAM::COMPRESSIBLE_FLUID_PARTICLE<PhysBAM::VECTOR<double,3> > *compressible_fluid_particle;
  PhysBAM::TETRAHEDRALIZED_VOLUME<double> *tetrahedralized_volume;
  PhysBAM::SOLIDS_PARTICLE<PhysBAM::VECTOR<double,3> > *solids_particle;
  PhysBAM::LIST_ARRAY<PhysBAM::VECTOR<int,4> > *elem_list;
  PhysBAM::TETRAHEDRON_MESH *tetrahedron_mesh;
  
//  int* CFPtoGlobNodeMap;

public:
  EulerStructGhostFluid(SubDomain*, Vec<double> &, SVec<double,3> &, Vec<bool> &, Vec<bool> &); 
  ~EulerStructGhostFluid();
  void getTriangulatedSurfaceFromFace(SVec<double,3> &);  //get coords of triangles.
  void getTriangulatedSurfaceFromFace();  //get node numbers of triangles.
  void getPhiFromModule(SVec<double,3> &, bool);
  void getPhiFromModule(SVec<double,3> &, double, double, double, double, double, double, bool);
  Vec3D getGradPhiAtPosition(Vec3D);
  void getGradPhiFromModule(SVec<double,3> &);
  void getGradPhiFromModule(SVec<double,3> &, double, double, double, double, double, double);
  int (*getTriangleNodeNumFromTriaSurf() const)[3];
  int getNumTriangleFromTriaSurf();
  double specifyBandwidth();
  void constructInterface(PhysBAM::GRID_3D<double>);
  void computeLevelSet(PhysBAM::TRIANGULATED_SURFACE<double> &);
  void initializePhysBAMMPI(int,double,int);
  void computeGhostNodes();
  void updateRealNodeTag(SVec<double,3> &,double,double,double,double,double,double);
  void recomputeLevelSet(Vec3D);

  bool insideOutside(double *position, const double xmin, const double xmax,
                     const double ymin, const double ymax, const double zmin, const double zmax);//check if a point is inside or outside a box. used in getTetNearInterface.
  bool insideOutside(Vec3D position, const double xmin, const double xmax,
                     const double ymin, const double ymax, const double zmin, const double zmax);
 
  template<int dim>
  void getTetNearInterface(SVec<double,3> &, const double, const double, const double,
                           const double, const double, const double, SVec<double,dim> &, 
                           const double, bool);
  template<int dim>
  void updateGhostNodes(SVec<double,dim> &, const double, Vec3D);
  template<int dim>
  void updateCFP(SVec<double,dim>&, int, int);
  template<int dim>
  void computeTotalForce(Vec3D&, int(*)[3], int, Vec3D*, double, SVec<double,3> &, SVec<double,dim> &);
  template<int dim>
  void populateNewFluid(SVec<double,3> &, SVec<double,dim> &);
};

#ifdef TEMPLATE_FIX
#include <EulerStructGhostFluid.C>
#endif

#endif
