#include <cstdio>
#include "Ghost/EulerStructGhostFluid.h"
//---------------------------------------------------------------------

EulerStructGhostFluid::EulerStructGhostFluid(Vec<double> & phi_input, SVec<double,3> &gradphi_input, int numNodes, int numElems): x(gradphi_input), philevel(phi_input), gradPhilevel(gradphi_input)
{
  PhysBAM_Interface = 0;

  nodeTag = new int[numNodes];
  numChosenNodes = 0;
  tempNodeList = new int[numNodes];
  numChosenElems = 0;
  tempElemList = new int[numElems][4];

  tetrahedralized_volume = 0;
  solids_particle = 0;
  elem_list = 0;
  tetrahedron_mesh = 0;
  compressible_fluid_particle = 0;

  totalForce[0] = totalForce[1] = totalForce[2] = 0.0;
}

//---------------------------------------------------------------------------------

EulerStructGhostFluid::~EulerStructGhostFluid()
{
  if (PhysBAM_Interface) delete PhysBAM_Interface;
  if (nodeTag) delete[] nodeTag;
  if (tempNodeList) delete[] tempNodeList;
  if (tempElemList) delete[] tempElemList;
}

//---------------------------------------------------------------------------------

void EulerStructGhostFluid::initializePhysBAM(PhysBAM::GRID_3D<double> grid3d)
{
  if (PhysBAM_Interface) delete PhysBAM_Interface;
  PhysBAM_Interface = new PhysBAM::AERO_INTERFACE_1<double> (grid3d);
}

//---------------------------------------------------------------------------------

void EulerStructGhostFluid::initializePhysBAMMPI(int* locToGlob, int numNodes,double depth)
{
  if (!tetrahedralized_volume){fprintf(stderr,"ERROR: haven't built mesh G'. Aborting...\n"); exit(-1);}
  PhysBAM::ARRAY<int> CFPtoGlobNodeMap(numChosenNodes);
  for (int iParticle=1; iParticle<=numChosenNodes; iParticle++)
    CFPtoGlobNodeMap(iParticle)=locToGlob[tempNodeList[iParticle-1]]+1;

//  PhysBAM::VECTOR<double,3> dimCartMeshMin, dimCartMeshMax;
//  PhysBAM_Interface->Initialize_MPI(*tetrahedralized_volume,CFPtoGlobNodeMap,numNodes,depth, dimCartMeshMin, dimCartMeshMax);
  PhysBAM_Interface->Initialize_MPI(*tetrahedralized_volume,CFPtoGlobNodeMap,numNodes,depth);
}

//--------------------------------------------------------------------------------------------

// PhysBAM::TRIANGULATED_SURFACE<double> contains the solid surface triangles
void EulerStructGhostFluid::computeLevelSet(PhysBAM::TRIANGULATED_SURFACE<double> &physbam_triangulated_surface)
{
  if (!PhysBAM_Interface) {fprintf(stderr,"ERROR: PhysBAM not initialized. Aborting.\n"); exit(-1);}
//  if (!compressible_fluid_particle) {fprintf(stderr,"Particle list not constructed. Aborting.\n"); exit(-1);}
//  PhysBAM_Interface->Compute_Level_Set(physbam_triangulated_surface, *compressible_fluid_particle);
  PhysBAM_Interface->Compute_Level_Set(physbam_triangulated_surface);
}

//--------------------------------------------------------------------------------------------

void EulerStructGhostFluid::getPhiFromModule(SVec<double,3> &X, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
  for (int i=0; i<X.size(); i++){
    philevel[i] = 10.0;
    if (X[i][0]<xmin || X[i][0]>xmax || X[i][1]<ymin || X[i][1]>ymax || X[i][2]<zmin || X[i][2]>zmax)
      continue;
    PhysBAM::VECTOR<double,3> position(X[i][0], X[i][1], X[i][2]);
    philevel[i] = PhysBAM_Interface->Phi(position);
  }
}

//-----------------------------------------------------------------------------------------

void EulerStructGhostFluid::getGradPhiFromModule(SVec<double,3> &X, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
  for (int i=0; i<X.size(); i++){
    gradPhilevel[i][0] = gradPhilevel[i][1] = gradPhilevel[i][2] = 0.0;

    if (X[i][0]<xmin || X[i][0]>xmax || X[i][1]<ymin || X[i][1]>ymax || X[i][2]<zmin || X[i][2]>zmax)
      continue;

    PhysBAM::VECTOR<double,3> position(X[i][0], X[i][1], X[i][2]);
    PhysBAM::VECTOR<double,3> gradPhi = PhysBAM_Interface->Gradient(position);
    gradPhilevel[i][0] = gradPhi[1]; gradPhilevel[i][1] = gradPhi[2]; gradPhilevel[i][2] = gradPhi[3];
  }
}

//----------------------------------------------------------------------------------------

Vec3D EulerStructGhostFluid::getGradPhiAtPosition(Vec3D pos)
{
  PhysBAM::VECTOR<double,3> position(pos[0], pos[1], pos[2]);
  PhysBAM::VECTOR<double,3> gradPhi = PhysBAM_Interface->Gradient(position);
  Vec3D grad(gradPhi[1], gradPhi[2], gradPhi[3]);
  return grad;
}

//-----------------------------------------------------------------------------------------

void EulerStructGhostFluid::getTetNearInterface(SVec<double,3> &subX)
{
  if (!numChosenNodes || !numChosenElems || !tempNodeList || !tempElemList) return; //outside bounding box.
  //1.construct SOLIDS_PARTICLE.
  if (solids_particle) delete solids_particle;
  solids_particle = new PhysBAM::SOLIDS_PARTICLE<PhysBAM::VECTOR<double,3> >;
  solids_particle->Add_Particles(numChosenNodes);
  for (int iParticle=0; iParticle<numChosenNodes; iParticle++) {
    double* position = subX[tempNodeList[iParticle]];
    solids_particle->X(iParticle+1) = PhysBAM::VECTOR<double,3>(position[0], position[1], position[2]);
  }

  //2.construct LIST_ARRAY<VECTOR<int,4> >. (tetrahedron_list)
  if (elem_list) delete elem_list;
  elem_list = new PhysBAM::LIST_ARRAY<PhysBAM::VECTOR<int,4> >;
  for (int iElem=0; iElem<numChosenElems; iElem++)
    elem_list->Append(PhysBAM::VECTOR<int,4>(tempElemList[iElem][0]+1,tempElemList[iElem][1]+1,
                                             tempElemList[iElem][2]+1,tempElemList[iElem][3]+1));

  //3.construct TETRAHEDRON_MESH and then TETRAHEDRALIZED_VOLUME.
  if (tetrahedron_mesh) delete tetrahedron_mesh;
  tetrahedron_mesh = new PhysBAM::TETRAHEDRON_MESH(numChosenNodes, *elem_list);
  if (tetrahedralized_volume) delete tetrahedralized_volume;
  tetrahedralized_volume = new PhysBAM::TETRAHEDRALIZED_VOLUME<double>(*tetrahedron_mesh, *solids_particle);

  //5.construct Compressible_Fluid_Particle_List.
  if (compressible_fluid_particle) delete compressible_fluid_particle;
  compressible_fluid_particle = new PhysBAM::COMPRESSIBLE_FLUID_PARTICLE<PhysBAM::VECTOR<double,3> > ;
  compressible_fluid_particle->Add_Particles(numChosenNodes);
  for (int iParticle=0; iParticle<numChosenNodes; iParticle++) {
    if (nodeTag[tempNodeList[iParticle]] != iParticle) {
      fprintf(stderr,"error in constructing COMPRESSIBLE_FLUID_PARTICLE!\n");
      exit(1);
    }
    double* position = subX[tempNodeList[iParticle]];
    compressible_fluid_particle->X(iParticle+1) = PhysBAM::VECTOR<double,3>(position[0], position[1], position
[2]);
  }

}

//------------------------------------------------------------------------------------
LevelSetResult
EulerStructGhostFluid::getLevelSetDataAtEdgeCenter(double t, int ni, int nj) {
    Vec3D center(0.5*(x[ni][0]+x[nj][0]), 0.5*(x[ni][1]+x[nj][1]),
                     0.5*(x[ni][2]+x[nj][2]));
    Vec3D gp = getGradPhiAtPosition(center);
    return LevelSetResult(gp[0], gp[1], gp[2], 0, 0, 0);
}


bool EulerStructGhostFluid::isActive(double t, int n) {
   return philevel[n] >= 0;
}




bool EulerStructGhostFluid::edgeIntersectsStructure(double t, int ni, int nj) const {
   return (philevel[ni] >= 0 && philevel[nj] <0) || (philevel[ni] < 0 && philevel[nj] >= 0);//philevel[ni]*philevel[nj] < 0;
}





