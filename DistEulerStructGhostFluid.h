#ifndef _DIST_EULER_STRUCT_GHOST_FLUID_H_
#define _DIST_EULER_STRUCT_GHOST_FLUID_H_

#include <string.h>
#include <EulerStructGhostFluid.h>
#include <DistVector.h>
#include <VarFcn.h>
#include "../PhysBAM/Projects/Charbel/AERO_INTERFACE_1.h"
#include "../PhysBAM/Public_Library/Geometry/TETRAHEDRALIZED_VOLUME.h"
#include <IoData.h>

class SubDomain;
class Domain;
class VarFcn;
class Timer;

class DistEulerStructGhostFluid {

  EulerStructGhostFluid **subESGF;
  SubDomain** subDomain;

  bool givenBB; //given boundingbox or not?
  int numLocSub;
  char* solidsurface; //solid surface file.
  double Xmin, Xmax, Ymin, Ymax, Zmin, Zmax; //coord. of bounding box.
  double bandwidth; 

  const int numGlobNodes;

  int (*triangle_list)[3]; // list of triangles on the interface.
  int length_triangle_list; //length of triangle_list!
  Vec3D* solids_particle_list; //particles(nodes) on the interface.
  int length_solids_particle_list; //Note: should aloways be the length of solids_particle_list!
  int* tag;

  DistVec<double> *philevel;
  DistSVec<double,3> *gradPhilevel;

  Timer* timer;
  Communicator *com;

  FILE* forceFile;  
protected:
  void specifyBoundingBox(DistSVec<double,3>*);
  double specifydx(Domain*, DistSVec<double,3>*);
  void prepareForCommunication(); //get triangle list, solids_particle_list directly from file.
  void initializePhysBAM(int, int, int, double, double, double, double, double, double);//PhysBAM_Interf
  double specifyBandwidth();
  void initializePhysBAMMPI();
  void getTetNearInterface(DistSVec<double,3> &, const double, const double, const double, const double, const double, const double); //construct tet_volume.

  void computeLevelSet();//construct the triangulated_surface and call function Compute_Level_Set.
  void getPhiFromModule(DistSVec<double,3> &, double, double, double, double, double, double);
  void getGradPhiFromModule(DistSVec<double,3> &, double, double, double, double, double, double);

  void getMinAndMax(Vec3D*, int& size, double*output);
  bool checkTriangulatedSurface();


public:
  double totalForce[3];
  double pref;
  
  DistEulerStructGhostFluid(Domain*, IoData&);
  ~DistEulerStructGhostFluid();

  EulerStructGhostFluid &operator() (int i) const {return *subESGF[i];}
  EulerStructGhostFluid* getSubESGFPointer(int i) {return subESGF[i];}

  DistVec<double>* getPhilevelPointer() {return philevel;}
  DistSVec<double,3>* getGradPhilevelPointer() {return gradPhilevel;}

  int getClosestStructureFace(Vec3D position, Vec3D& projection, double& distance)
  {return (subESGF[0]->getClosestStructureFace(position, projection, distance));} //currently only works for 1 proc.

  void clearTotalForce();
  Vec3D getTotalForce();

/*  void clearTotalForce() 
  {
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; iSub++) 
      subESGF[iSub]->totalForce[0] =  subESGF[iSub]->totalForce[1] = subESGF[iSub]->totalForce[2] = 0.0;
    totalForce[0] = totalForce[1] = totalForce[2] = 0.0;
  }

  Vec3D getTotalForce()  
  {
    double subDomainTotalForce[numLocSub][3];
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; iSub++) {
      subDomainTotalForce[iSub][0] = subESGF[iSub]->totalForce[0];
      subDomainTotalForce[iSub][1] = subESGF[iSub]->totalForce[1];
      subDomainTotalForce[iSub][2] = subESGF[iSub]->totalForce[2];
    }
    for (int iSub=0; iSub<numLocSub; iSub++) {
      totalForce[0] += subDomainTotalForce[iSub][0];
      totalForce[1] += subDomainTotalForce[iSub][1];
      totalForce[2] += subDomainTotalForce[iSub][2];
    }
    com->globalSum(3, totalForce);

    com->fprintf(stderr,"Total Force on Structure Surface = [%f, %f, %f]\n", totalForce[0]*pref, totalForce[1]*pref, totalForce[2]*pref);
    com->fprintf(forceFile, "%lf %lf %lf\n", totalForce[0]*pref, totalForce[1]*pref, totalForce[2]*pref);
  }
*/
  template<int dim>
  void setupCommunication(Domain*, DistSVec<double,3>*, DistSVec<double,dim> &);

  

};

#ifdef TEMPLATE_FIX
#include <DistEulerStructGhostFluid.C>
#endif

#endif

