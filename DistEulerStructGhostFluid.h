// Introduction
// 1. Class "DistEulerStructGhostFluid" is an interface to communicate with PhysBAM for ghost solid simulation. It's stored in TsDesc. The interface on PhysBAM's side is AERO_INTERFACE_1. (PhysBAM is sometimes called a "module".)
// 2. The major functions of DistEulerStructGhostFluid includes
//    (1) load in triangulated surface from a file (prepareForCommunication). 
//    (2) create a Cartesian grid (constructInterface).
//    (3) compute Levelset and its gradient on Cartesian grid. (computeLevelset).
//    (4) interpolate levelset and its gradient from Cartesian grid to fluid grid. (getPhiFromModule and getGradPhiFromModule.) 
//    (5) construct a tet grid in PhysBAM's format by tets inside the bounding box (this grid is sometimes called G', or Gprime) (getTetNearInterface)
//    (6) update ghost nodes within a bandwidth by mirroring. (computeGhostNodes and updateGhostNodes)


#ifndef _DIST_EULER_STRUCT_GHOST_FLUID_H_
#define _DIST_EULER_STRUCT_GHOST_FLUID_H_

#include <string.h>
#include <EulerStructGhostFluid.h>
#include <DistVector.h>
#include <Domain.h>
#include <VarFcn.h>
#include "../PhysBAM/Projects/Charbel/AERO_INTERFACE_1.h"
#include "../PhysBAM/Public_Library/Geometry/TETRAHEDRALIZED_VOLUME.h"
#include <IoData.h>

class SubDomain;
class VarFcn;
class Timer;

class DistEulerStructGhostFluid {

//  PhysBAM::AERO_INTERFACE_1<double> *PhysBAM_Interface;  
  EulerStructGhostFluid **subESGF;
  SubDomain** subDomain;

  bool givenBB; //given boundingbox or not?
  bool forcedMotion;
  int numLocSub;
  char* solidsurface;
  double Xmin, Xmax, Ymin, Ymax, Zmin, Zmax;
  double bandwidth;

  const int numGlobNodes;
    
  int (*triangle_list)[3]; // list of triangles on the interface.
  int length_triangle_list; //Note: should always be the length of triangle_list!
  Vec3D* solids_particle_list; 
 //double (*solids_particle_list)[3]; //particles(nodes) on the interface.
  int length_solids_particle_list; //Note: should aloways be the length of solids_particle_list!
  int* tag;
  
  DistVec<double> *philevel;
  DistSVec<double,3> *gradPhilevel; 
  DistSVec<double,1> *PhiWrite;
  
  Timer* timer;
  Communicator *com;

  // for Forced Motion with a constant translating velocity.
  Vec3D vel, dsp;
  DistVec<bool> *realNodeTag0;  //for the previous time-step
  DistVec<bool> *realNodeTag;   //for the current time-step


protected:
  void getBoundingBox(double &, double &, double &, double &, double &, double &);
  void specifyBoundingBox(DistSVec<double,3>*);
  double specifydx(Domain*, DistSVec<double,3>*);
  double specifyBandwidth();

  void getTriangulatedSurfaceFromFace(DistSVec<double, 3> &); // store coords under each subESGF
  void getTriangulatedSurfaceFromFace();  //store the node numbers under each subESGF
  void constructInterface(int, int, int, double, double, double, double, double, double);//PhysBAM_Interf
  void prepareForCommunication( DistSVec<double, 3> &); //get triangle_list, solids_particle_list
  void prepareForCommunication(); //get triangle list, solids_particle_list directly from file.

  void computeLevelSet();//construct the triangulated_surface and call function Compute_Level_Set.
  void initializePhysBAMMPI(int,double);
  void getPhiFromModule(DistSVec<double,3> &, bool);
  void getPhiFromModule(DistSVec<double,3> &, double, double, double, double, double, double, bool);
//  double getPhiFromModule(double, double, double);//get Phi at arbitrary position
  void getGradPhiFromModule(DistSVec<double,3> &);
  void getGradPhiFromModule(DistSVec<double,3> &, double, double, double, double, double, double);

  void computeGhostNodes(); //call PhysBAM to populate ghost nodes on compreesible_fluid_particle.

//  bool insideOutside(double *position, const double xmin, const double xmax,
//                     const double ymin, const double ymax, const double zmin, const double zmax);//check if a point is inside or outside a box. used in getTetNearInterface.
  void getMinAndMax(Vec3D*, int& size, double*output);
  bool checkTriangulatedSurface();
  bool updateStructureDynamics(double);
  void recomputeLevelSet(Vec3D);
  void updateRealNodeTag(DistSVec<double,3> &, double, double, double, double, double, double);

  template<int dim>
  void getTetNearInterface(DistSVec<double,3> &, const double, const double, const double, const double, const double, const double, DistSVec<double,dim> &); //construct tet_volume.
  template<int dim>
  void updateCFP(DistSVec<double,dim> &); //update rho, u, e on compressible_fluid_particle. for 1 CPU!!
  template<int dim>
  void updateGhostNodes(DistSVec<double,dim> &); //update U at ghost nodes. for 1 CPU only!
  template<int dim>
  void populateNewFluid(DistSVec<double,3> &, DistSVec<double,dim> &); 

public:

  DistEulerStructGhostFluid(Domain*, IoData&);
  ~DistEulerStructGhostFluid();

  DistVec<double>* getPhilevelPointer() {return philevel;}
  DistSVec<double,3>* getGradPhilevelPointer() {return gradPhilevel;}

  template<int dim>
  void setupCommunication(Domain*, DistSVec<double,3>*, DistSVec<double,dim> &);
  template<int dim>
  void updateGhostFluid(DistSVec<double,3>*, DistSVec<double,dim> &, Vec3D&, double);
//  template<int dim>
//  void updateGhostFluid(DistSVec<double,3>*, DistSVec<double,dim>&);//for moving solid.
  template<int dim>
  void computeTotalForce(Vec3D&, DistSVec<double,3> &, DistSVec<double,dim>&);
  
};

#ifdef TEMPLATE_FIX
#include <DistEulerStructGhostFluid.C>
#endif

#endif
