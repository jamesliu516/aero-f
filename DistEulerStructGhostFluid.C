#include <DistEulerStructGhostFluid.h>
#include <stdio.h>

//-----------------------------------------------------------------------------------
/*template<int dim>
void DistEulerStructGhostFluid::mirroring(DistSVec<double, dim> &U)
{
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->mirroring(U(iSub));
}
*/
//------------------------------------------------------------------------------------
template<int dim>
void DistEulerStructGhostFluid::getTetNearInterface(DistSVec<double,3> &X, 
                                                    const double xmin, const double xmax, const double ymin,
                                                    const double ymax, const double zmin, const double zmax,
                                                    DistSVec<double,dim> &U)
{
  fprintf(stderr, "DistEulerStructGhostFluid::getTetNearInterface called.\n");
  specifyBandwidth();  bandwidth*= 2.0; //TODO: should be removed later. doesn't impact solutions.
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->getTetNearInterface(X(iSub), xmin, xmax, ymin, ymax, zmin, zmax, U(iSub), bandwidth, true);
}

//-----------------------------------------------------------------------------------

template<int dim>
void DistEulerStructGhostFluid::updateCFP(DistSVec<double,dim> &U) //update rho, u, e on compressible_fluid_particle. 
{
  fprintf(stderr,"DistEulerStructGhostFluid::updateCFP called.\n");
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->updateCFP(U(iSub), com->cpuNum(), com->size());
}

//-----------------------------------------------------------------------------------

template<int dim>
void DistEulerStructGhostFluid::updateGhostNodes(DistSVec<double,dim> &U) //update U at ghost cells. 
{
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->updateGhostNodes(U(iSub), bandwidth, vel);

  fprintf(stderr,"Ghost fluid state vectors updated.\n");

}

//------------------------------------------------------------------------------------

template<int dim>
void DistEulerStructGhostFluid::updateGhostFluid(DistSVec<double,3> *X, DistSVec<double,dim> &U, Vec3D& totalForce, double dt)
{
    bool structure_moved = updateStructureDynamics(dt);
    if (structure_moved) {
      recomputeLevelSet(dt*vel);  //pass in the displacement within one time-step.
      getPhiFromModule(*X,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax,true); //interpolate phi onto tet grid.
      getTetNearInterface(*X, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, U); //find tet within a band.
      populateNewFluid(*X,U); //go over "tetNearInterface". 
    }
   
       
    updateCFP(U);
    computeGhostNodes();
    updateGhostNodes(U); 

//   fprintf(stderr,"\n");
//    computeTotalForce(totalForce, *X, U);
//    fprintf(stderr,"***************** FORCE ******************\n");
//    fprintf(stderr,"(fx,fy,fz) = %f, %f, %f.\n", totalForce[0], totalForce[1], totalForce[2]);
//    fprintf(stderr,"******************************************\n"); 
}

//---------------------------------------------------------------------------------------

template<int dim>
void DistEulerStructGhostFluid::setupCommunication(Domain *domain, DistSVec<double,3>*X, DistSVec<double,dim> &U)
{
    if (vel.norm()) fprintf(stderr,"FORCED MOTION WITH CONSTANT VELOCITY: velocity = (%f, %f, %f).\n", vel[0], vel[1], vel[2]);  
    if (solidsurface)  prepareForCommunication();
    else {getTriangulatedSurfaceFromFace(); prepareForCommunication(*X); }
    if (!givenBB) specifyBoundingBox(X);
    double dx = specifydx(domain,X);
//    int m = static_cast<int>((Xmax - Xmin) / dx);
//    int n = static_cast<int>((Ymax - Ymin) / dx);
//    int mn = static_cast<int>((Zmax - Zmin) / dx);
    int m,n,mn;
    m = 128; n = mn = 43;   //TODO: this is just for debug. need to get a global dx.
    fprintf(stderr,"For Cart. Mesh, dx = %lf, dy = %lf, dz = %lf.\n", (Xmax-Xmin)/m, (Ymax-Ymin)/n, (Zmax-Zmin)/mn);
    constructInterface(m,n,mn,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax);

    double time = timer->getTime();
    computeLevelSet();
    time = timer->getTime() - time;
    fprintf(stderr,"Compute Levelset Time: %f.\n", time);

    getPhiFromModule(*X,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax,true);
    updateRealNodeTag(*X,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax); // called twice to fill both realNodeTag0 and realNodeTag.

    getTetNearInterface(*X, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, U);

    fprintf(stderr,"BANDWIDTH: %lf. numGlobNodes = %d.\n", bandwidth, numGlobNodes);
    
    initializePhysBAMMPI(numGlobNodes,bandwidth); //send the total number of nodes on the entire grid.
//    initializePhysBAMMPI(1000,bandwidth);
    computeGhostNodes();

    updateGhostNodes(U);

}

//----------------------------------------------------------------------------------------

template<int dim>
void DistEulerStructGhostFluid::computeTotalForce(Vec3D &forceToWrite, DistSVec<double,3> &X, DistSVec<double,dim> &U)
{
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->computeTotalForce(forceToWrite, triangle_list, length_triangle_list, solids_particle_list, length_solids_particle_list, X(iSub), U(iSub));  //TODO:now only works for 1-proc.

}

//----------------------------------------------------------------------------------------

template<int dim>
void DistEulerStructGhostFluid::populateNewFluid(DistSVec<double,3> &X, DistSVec<double,dim> &U)
{
  fprintf(stderr,"DistEulerStructGhostFluid::populateNewFluid(...) called.\n");
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->populateNewFluid(X(iSub), U(iSub));  
}



















