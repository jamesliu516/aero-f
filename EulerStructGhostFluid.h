#ifndef _EULER_STRUCT_GHOST_FLUID_H_
#define _EULER_STRUCT_GHOST_FLUID_H_

#include <Vector.h>
#include <Vector3D.h>
#include <DistVector.h>
#include "AERO_INTERFACE_1.h"

#include "LevelSetStructure.h"

template<class Scalar, int dim> class SVec;
class EulerStructGhostFluid : public LevelSetStructure {

  Vec<double> &philevel;
  SVec<double,3> &gradPhilevel;
  SVec<double, 3> &x;
  PhysBAM::AERO_INTERFACE_1<double> *PhysBAM_Interface;

  PhysBAM::TETRAHEDRALIZED_VOLUME<double> *tetrahedralized_volume;
  /* Coordinates of nodes of the structure */
  PhysBAM::SOLIDS_PARTICLE<PhysBAM::VECTOR<double,3> > *solids_particle;
  PhysBAM::LIST_ARRAY<PhysBAM::VECTOR<int,3> > *triangle_list;
  PhysBAM::LIST_ARRAY<PhysBAM::VECTOR<int,4> > *elem_list;
  PhysBAM::TETRAHEDRON_MESH *tetrahedron_mesh;
  PhysBAM::COMPRESSIBLE_FLUID_PARTICLE<PhysBAM::VECTOR<double,3> > *compressible_fluid_particle;

public:
  int* nodeTag; //tag nodes in bounding box.
  int numChosenNodes;
  int* tempNodeList;
  int numChosenElems;
  int (*tempElemList)[4];  //stores the global index of chosen elems.


  EulerStructGhostFluid(Vec<double> &, SVec<double,3> &, int, int);
  ~EulerStructGhostFluid();

  void initializePhysBAM(PhysBAM::GRID_3D<double>);
  void initializePhysBAMMPI(int*, int,double);
  void computeLevelSet(PhysBAM::TRIANGULATED_SURFACE<double> &);
  void getPhiFromModule(SVec<double,3> &, double, double, double, double, double, double);
  void getGradPhiFromModule(SVec<double,3> &, double, double, double, double, double, double);
  Vec3D getGradPhiAtPosition(Vec3D);
  void getTetNearInterface(SVec<double,3> &);

  Vec<double> *getPhilevelPointer() {return &philevel;}
  SVec<double,3>* getGradPhilevelPointer() {return &gradPhilevel;}
  int getClosestStructureFace(Vec3D position, Vec3D& projection, double& distance) {
/*  //commented because not compatible with current version of PhysBAM  
    if (!PhysBAM_Interface) {fprintf(stderr,"ERROR:PhysBAM_Interface not initialized!\n"); return -1;}
    PhysBAM::VECTOR<double,3> physbam_position(position[0],position[1],position[2]);
    PhysBAM::VECTOR<double,3> physbam_projection;
    int index;
    index = PhysBAM_Interface->Closest_Boundary_Point(physbam_position, physbam_projection, distance);
    for (int i=0; i<3; i++) projection[i] = physbam_projection[i+1];
    return index;*/
    return 0;
  }

  LevelSetResult
       getLevelSetDataAtEdgeCenter(double t, int ni, int nj);
  double phiAtNode(double t, int n);
};

#ifdef TEMPLATE_FIX
#include <EulerStructGhostFluid.C>
#endif

#endif





