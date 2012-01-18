#include <LevelSet.h>
#include <fstream>
#include <cmath>

#include "IoData.h"
#include "Domain.h"
#include "TimeData.h"
#include <LevelSetStructure.h>

#include <OneDimensionalSolver.h>
using namespace std;

//-------------------------------------------------------------------------
template<int dimLS>
LevelSet<dimLS>::LevelSet(IoData &iod, Domain *dom):
  Phin(dom->getNodeDistInfo()), Phinm1(dom->getNodeDistInfo()), Phinm2(dom->getNodeDistInfo()),
  Psi(dom->getNodeDistInfo()), dt(dom->getNodeDistInfo()), PsiRes(dom->getNodeDistInfo()),
  w(dom->getNodeDistInfo()), domain(dom), Phi0(dom->getNodeDistInfo()),
  Tag(dom->getNodeDistInfo())
{

  com = domain->getCommunicator();
  numLocSub = domain->getNumLocSub();

  data = new TimeData(iod);

  bandlevel = iod.mf.bandlevel;
  localtime = iod.mf.localtime;
  subIt = iod.mf.subIt;
  cfl_psi = iod.mf.cfl;
  typeTracking = iod.mf.typeTracking;
  copy = iod.mf.copy;
  conv_eps = iod.mf.eps;
  diff = bool(iod.mf.outputdiff);
  if(typeTracking == MultiFluidData::GRADIENT){
    com->fprintf(stdout, "***Warning: if reinitialization in band --> problem!\n ***         You need to reinitialize in the whole domain with this method\n");
  }

  lsgrad  = new DistNodalGrad<dimLS,double>(iod,dom,2);

}
//-------------------------------------------------------------------------
template<int dimLS>
LevelSet<dimLS>::~LevelSet()
{
  if(data) delete data;
  if(lsgrad) delete lsgrad;
}

//---------------------------------------------------------------------------------------------------------
template<int dimLS>
template<int dim>
void LevelSet<dimLS>::setup(const char *name, DistSVec<double,3> &X, DistSVec<double,dim> &U,
                            DistSVec<double,dimLS> &Phi, IoData &iod, FluidSelector* fs, VarFcn* vf,
                            DistVec<ClosestPoint> *closest, DistVec<int> *fsId,int lsMethod)
{

  map<int, FluidModelData *>::iterator it = iod.eqs.fluidModelMap.dataMap.find(1);
  if(it == iod.eqs.fluidModelMap.dataMap.end()){
    fprintf(stderr, "*** Error: no FluidModel[1] was specified\n");
    exit(1);
  }
  if((iod.eqs.fluidModel.fluid  == FluidModelData::PERFECT_GAS ||
      iod.eqs.fluidModel.fluid  == FluidModelData::STIFFENED_GAS) &&
     it->second->fluid == FluidModelData::LIQUID)
    invertGasLiquid = -1.0;
  else invertGasLiquid = 1.0;

  // CONVENTION: phi >= 0 --- fluid1 --- 'normal' volID
  //             phi <  0 --- fluid2 --- 'special' volID
  // Initialization of Phi is done through the use of volumeID and
  // through the knowledge of a geometric shape (with its position).

  Phi0 = -1.0;
  setupPhiVolumesInitialConditions(iod, Phi0);
  setupPhiOneDimensionalSolution(iod,X,U,Phi0,fs,vf);
  setupPhiMultiFluidInitialConditions(iod,X, Phi0);
  if(closest && fsId)
    setupPhiFluidStructureInitialConditions(iod,X,Phi0,*closest,*fsId);
  if (lsMethod == 0)
    primitiveToConservative(Phi0, Phi, U);
  else
    Phi = Phi0;

  Phin   = Phi;
  Phinm1 = Phin;
  Phinm2 = Phinm1;

  if (name[0] != 0) {
    DistSVec<double,dimLS> ReadPhi(domain->getNodeDistInfo());
    domain->readVectorFromFile(name, 0, 0, ReadPhi);
    Phi  = ReadPhi;
    Phin = ReadPhi;

    if (data->use_nm1){
      DistSVec<double,dimLS> ReadPhi1(domain->getNodeDistInfo());
      data->exist_nm1 = domain->readVectorFromFile(name, 1, 0, ReadPhi1);
      Phinm1 = ReadPhi1;
    }

    if (data->use_nm2){
      DistSVec<double,dimLS> ReadPhi2(domain->getNodeDistInfo());
      data->exist_nm2 = domain->readVectorFromFile(name, 2, 0, ReadPhi2);
      Phinm2 = ReadPhi2;
    }
  }

  if (data->use_nm1 && !data->exist_nm1){
    Phinm1 = Phin;
  }
  if (data->use_nm2 && !data->exist_nm2){
    Phinm2 = Phinm1;
  }

  // determine which level-sets must be updated and reinitialized
  //exit(-1);
  if (lsMethod == 0)
    conservativeToPrimitive(Phi,Phi0,U);
  else
    Phi0 = Phi;

  double minDist[dimLS];
  double maxDist[dimLS];
  Phi0.min(minDist);
  Phi0.max(maxDist);
  for(int idim=0; idim<dimLS; idim++){
    if(fabs(minDist[idim])==fabs(maxDist[idim]) && fabs(minDist[idim]) == 1)
      trueLevelSet[idim] = false;
    else trueLevelSet[idim] = true;
    com->fprintf(stdout, "minDist[%d] = %e and maxDist[%d] = %e ==> %d\n", idim, minDist[idim], idim, maxDist[idim], trueLevelSet[idim]);
  }

  if(closest && fsId) {//cracking...
    if(dimLS!=1){fprintf(stderr,"ERROR: Multi-Phase FSI w/ Cracking supports only one level-set! dimLS = %d.\n", dimLS);exit(-1);}
    trueLevelSet[0] = true;
  }

  // for reinitialization testing
  Phi0 = Phi;

}

//-------------------------------------------------------------------------

template<int dimLS>
void LevelSet<dimLS>::setupPhiVolumesInitialConditions(IoData &iod, DistSVec<double,dimLS> &Phi){

  // loop on all Volumes to setup Phi_0 
  // Phi_0[iPhase] is 1 where VarFcn[iPhase+1] is present and -1 elsewhere
  // except for the first phase (varFcn[0]) that corresponds to all previous Phi being negative
  if(!iod.volumes.volumeMap.dataMap.empty()){
    map<int, VolumeData *>::iterator volIt;
    for (volIt=iod.volumes.volumeMap.dataMap.begin(); volIt!=iod.volumes.volumeMap.dataMap.end();volIt++)
      if(volIt->second->type==VolumeData::FLUID){
        com->fprintf(stdout, "Processing initialization of LevelSet for volume %d paired with EOS %d\n", volIt->first, volIt->second->fluidModelID);
        domain->setupPhiVolumesInitialConditions(volIt->first, volIt->second->fluidModelID,Phi);
      }
  }

}

//---------------------------------------------------------------------------------------

template<int dimLS>
template<int dim>
void LevelSet<dimLS>::setupPhiOneDimensionalSolution(IoData &iod, DistSVec<double,3> &X, DistSVec<double,dim> &U, DistSVec<double,dimLS> &Phi, FluidSelector* fs, VarFcn* vf){

  OneDimensional::read1DSolution(iod, U, 
				 &Phi, fs,
				 vf, X,*domain,
				 OneDimensional::ModePhi);
}

static inline double min(double a,double b, double c, double d, double e, double f) {

  if (a <= b && a <= c && a <= d && a <= e && a <= f)
    return a;
  else if (b <= c && b <= d && b <= e && b <= f)
    return b;
  else if (c <= d && c <= e && c <= f)
    return c;
  else if (d <= e && d <= f)
    return d;
  else if (e <= f)
    return e;
  else
    return f;  
}

static inline double maxp(double a,double b, double c,double as,double bs,double cs) {

  if (a >= b && a >= c)
    return a*as;
  else if (b >= c)
    return b*bs;
  else
    return c*cs;
}
//---------------------------------------------------------------------------------------

template<int dimLS>
void LevelSet<dimLS>::setupPhiMultiFluidInitialConditions(IoData &iod, DistSVec<double,3> &X, DistSVec<double,dimLS> &Phi){

  // Note that Phi was already initialized either to -1.0 everywhere
  // (if there were no volumes specified in the input file)
  // or to +1.0 and to -1.0.

  // Assumption: one geometric object surface does not cross another geometric
  //             object surface. In other words, there cannot exist a point
  //             in space where three or more fluids are in contact.
  // Assumption: a geometric surface always separates the main fluid
  //             given by varFcn[0] = iod.eqs.fluidModel =
  //             iod.eqs.fluidModelMap.dataMap.begin() from
  //             another fluid, at least one of the fluid involved is the main fluid

  DistSVec<double,dimLS> Distance(domain->getNodeDistInfo());
  Distance = -1.0;

  if(!iod.mf.multiInitialConditions.planeMap.dataMap.empty()){
    map<int, PlaneData *>::iterator planeIt;
    for(planeIt  = iod.mf.multiInitialConditions.planeMap.dataMap.begin();
        planeIt != iod.mf.multiInitialConditions.planeMap.dataMap.end();
        planeIt++){
      com->fprintf(stdout, "Processing initialization of LevelSet for plane %d paired with EOS %d\n", planeIt->first, planeIt->second->fluidModelID);
      double nx = planeIt->second->nx;
      double ny = planeIt->second->ny;
      double nz = planeIt->second->nz;
      double norm = sqrt(nx*nx+ny*ny+nz*nz);
      nx /= norm; ny /= norm; nz /= norm;

#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {

        SVec<double,dimLS> &phi(Phi(iSub));
        SVec<double,dimLS> &distance(Distance(iSub));
        SVec<double,3>     &x  (X(iSub));

        for (int i=0; i<x.size(); i++){
          double scalar = nx*(x[i][0]-planeIt->second->cen_x)+
                          ny*(x[i][1]-planeIt->second->cen_y)+
                          nz*(x[i][2]-planeIt->second->cen_z); // positive in the direction of n
          // set level-set only if this fluid has its own level-set (nothing to do for fluidModelID = 0)
          if(planeIt->second->fluidModelID > 0){
            int fluidId = planeIt->second->fluidModelID-1;
            if(scalar>0.0) phi[i][fluidId] = 1.0;
            if(distance[i][fluidId]<0.0) distance[i][fluidId] = fabs(scalar);
            else                         distance[i][fluidId] = fmin(fabs(scalar),distance[i][fluidId]);
          }
        }
      }
    }
  }

  if(!iod.mf.multiInitialConditions.sphereMap.dataMap.empty()){
    map<int, SphereData *>::iterator sphereIt;
    for(sphereIt  = iod.mf.multiInitialConditions.sphereMap.dataMap.begin();
        sphereIt != iod.mf.multiInitialConditions.sphereMap.dataMap.end();
        sphereIt++){
      com->fprintf(stdout, "Processing initialization of LevelSet for sphere %d paired with EOS %d\n", sphereIt->first, sphereIt->second->fluidModelID);
      double x0 = sphereIt->second->cen_x;
      double y0 = sphereIt->second->cen_y;
      double z0 = sphereIt->second->cen_z;
      double r0 = sphereIt->second->radius;
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {

        SVec<double,dimLS> &phi(Phi(iSub));
        SVec<double,dimLS> &distance(Distance(iSub));
        SVec<double,3>     &x  (X(iSub));

        for (int i=0; i<x.size(); i++){
          double scalar = sqrt((x[i][0] - x0)*(x[i][0] - x0) + 
                          (x[i][1] - y0)*(x[i][1] - y0) + 
                          (x[i][2] - z0)*(x[i][2] - z0)) - r0;//positive outside the sphere
          if(sphereIt->second->fluidModelID > 0){
            int fluidId = sphereIt->second->fluidModelID-1;
            if(scalar<0.0) phi[i][fluidId] = 1.0;
            if(distance[i][fluidId]<0.0) distance[i][fluidId] = fabs(scalar);
            else                         distance[i][fluidId] = fmin(fabs(scalar),distance[i][fluidId]);
          }
        }
      }
    }
  }

  if(!iod.mf.multiInitialConditions.prismMap.dataMap.empty()){
    map<int, PrismData *>::iterator prismIt;
    for(prismIt  = iod.mf.multiInitialConditions.prismMap.dataMap.begin();
        prismIt != iod.mf.multiInitialConditions.prismMap.dataMap.end();
        prismIt++){
      //com->fprintf(stdout, "Processing initialization of LevelSet for sphere %d paired with EOS %d\n", sphereIt->first, sphereIt->second->fluidModelID);
      double x0 = prismIt->second->cen_x;
      double y0 = prismIt->second->cen_y;
      double z0 = prismIt->second->cen_z;
      double w_x = prismIt->second->w_x;
      double w_y = prismIt->second->w_y;
      double w_z = prismIt->second->w_z;
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {

        SVec<double,dimLS> &phi(Phi(iSub));
        SVec<double,dimLS> &distance(Distance(iSub));
        SVec<double,3>     &x  (X(iSub));

        for (int i=0; i<x.size(); i++){
          double s[3] = {(x[i][0]-x0)/w_x*0.5,(x[i][1]-y0)/w_y*0.5,(x[i][2]-z0)/w_z*0.5};
	  double scalar = maxp(fabs(s[0]),fabs(s[1]),fabs(s[2]),1.0,1.0,1.0);
          if(prismIt->second->fluidModelID > 0){
            int fluidId = prismIt->second->fluidModelID-1;
            if(prismIt->second->inside(x[i][0],x[i][1],x[i][2])) 
	      phi[i][fluidId] = 1.0;
	    if(distance[i][fluidId]<0.0)
	      distance[i][fluidId] = scalar;
	    else
	      distance[i][fluidId] = fmin( scalar, distance[i][fluidId]);
          }
        }
      }
    }
  }

  // phi[fluidId] is initialized with -1 and +1 giving the location of fluid fluidId
  // distance[fluidId] is the positive distance to the closest interface between fluid fluidId and fluid 0.
  double minDist[dimLS];
  double maxDist[dimLS];
  Distance.min(minDist);
  Distance.max(maxDist);
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {

    SVec<double,dimLS> &phi(Phi(iSub));
    SVec<double,dimLS> &distance(Distance(iSub));
    for (int idim=0; idim<dimLS; idim++){
      // multiply only when distance is not -1 everywhere
      if(minDist[idim] != maxDist[idim])
        for (int i=0; i<phi.size(); i++)
          phi[i][idim] *= distance[i][idim];
      //else
      //  fprintf(stdout, "LevelSet[%d] does not need to be updated since it seems enclosed in a volume!\n", idim);
    }
  }

}

//---------------------------------------------------------------------------------------

template<int dimLS>
void LevelSet<dimLS>::setupPhiFluidStructureInitialConditions(IoData &iod, DistSVec<double,3> &X, DistSVec<double,dimLS> &Phi, 
                                                              DistVec<ClosestPoint> &closest, DistVec<int> &status)
{
  trueLevelSet[0] = true; //this is a 'true' level-set.
  //initialize the level-set near FS interface  
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {
      SVec<double,dimLS> &phi(Phi(iSub));
      Vec<ClosestPoint> &clo(closest(iSub));
      Vec<int> &stat(status(iSub));
      for(int i=0; i<phi.size(); i++) {
        switch (stat[i]) {
          case 0 : //outside the structure
            phi[i][0] = clo[i].nearInterface() ? -1.0*clo[i].dist : -1.0;
            break;
          case 1 : //inside the structure
            phi[i][0] = clo[i].nearInterface() ?  1.0*clo[i].dist :  1.0;
            break;
          case 2 : //occluded
            phi[i][0] = 0.0;
            break;
          default :
            fprintf(stderr,"ERROR: Status cannot be %d.\n", stat[i]);
            exit(-1);
        }
      }
    }

  //call "reinitialize"
  reinitializeLevelSet(X, Phi, false);
}

//---------------------------------------------------------------------------------------------------------
template<int dimLS>
void LevelSet<dimLS>::checkTrueLevelSetUpdate(DistSVec<double,dimLS> &dPhi)
{
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {
      SVec<double,dimLS> &dphi(dPhi(iSub));
      for(int i=0; i<dphi.size(); i++)
        for(int idim=0; idim<dimLS; idim++)
          if(!trueLevelSet[idim])
            dphi[i][idim] = 0.0;
    }
}

//---------------------------------------------------------------------------------------------------------
template<int dimLS>
void LevelSet<dimLS>::update(DistSVec<double,dimLS> &Phi)
{

  if (data->use_nm2 && data->exist_nm1) {
    Phinm2 = Phinm1;
    data->exist_nm2 = true;
  }
  if (data->use_nm1) {
    Phinm1 = Phin;
    data->exist_nm1 = true;
  }
  Phin = Phi;
}

//---------------------------------------------------------------------------------------------------------
template<int dimLS>
void LevelSet<dimLS>::writeToDisk(char *name)
{

  domain->writeVectorToFile(name, 0, 0.0, Phin);

  if (data->use_nm1){
    domain->writeVectorToFile(name, 1, 0.0, Phinm1);
  }

  if (data->use_nm2){
    domain->writeVectorToFile(name, 2, 0.0, Phinm2);
  }

}
//---------------------------------------------------------------------------------------------------------
template<int dimLS>
template<int dim>
void LevelSet<dimLS>::conservativeToPrimitive(DistSVec<double,dimLS> &Cons, DistSVec<double,dimLS> &Prim, 
	                               DistSVec<double,dim> &U)
{
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub){
    double (*u)[dim] = U.subData(iSub);
    double (*prim)[dimLS] = Prim.subData(iSub);
    double (*cons)[dimLS] = Cons.subData(iSub);
    for (int i=0; i<U.subSize(iSub); i++)
      for (int idim=0; idim<dimLS; idim++)
        prim[i][idim] = cons[i][idim]/u[i][0];
  }

}

//-------------------------------------------------------------------------
template<int dimLS>
template<int dim>
void LevelSet<dimLS>::primitiveToConservative(DistSVec<double,dimLS> &Prim, DistSVec<double,dimLS> &Cons,
				       DistSVec<double,dim> &U)
{
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub){
    double (*u)[dim] = U.subData(iSub);
    double (*prim)[dimLS] = Prim.subData(iSub);
    double (*cons)[dimLS] = Cons.subData(iSub);
    for (int i=0; i<U.subSize(iSub); i++)
      for (int idim=0; idim<dimLS; idim++)
        cons[i][idim] = prim[i][idim]*u[i][0];
  }

}

//-------------------------------------------------------------------------
template<int dimLS>
void LevelSet<dimLS>::reinitializeLevelSet(DistSVec<double,3> &X, DistSVec<double,dimLS> &Phi, bool copylv2)
{
  // XXX reinitializeLevelSetPDE(geoState,X,ctrlVol,U,Phi);
  reinitializeLevelSetFM(X,Phi,copylv2);

}
//-------------------------------------------------------------------------
template<int dimLS>
void LevelSet<dimLS>::reinitializeLevelSetFM(DistSVec<double,3> &X, DistSVec<double,dimLS> &Phi, bool copylv2)
{

  // reinitialize each LevelSet separately
  for(int idim=0; idim<dimLS; idim++){
    if(!trueLevelSet[idim]) continue;
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub){
      double (*phi)[dimLS] = Phi.subData(iSub);
      double (*psi)[1]     = Psi.subData(iSub);
      for (int i=0; i<Phi.subSize(iSub); i++)
        psi[i][0] = phi[i][idim];
    }
    //com->fprintf(stdout, "reinitializing phi[%d] to distance (phi[%d].norm = %e)\n", idim, idim, Psi.norm());

    // initialize Psi
    Psi = 1.0;

    // tag nodes that are close to interface up to level 'bandlevel'
    int level = 0;
    for (level=0; level<bandlevel; level++)
      domain->TagInterfaceNodes(idim,Tag,Phi,level);


    //compute distance for Psi<=0.5 and level<bandlevel
    // first, the nodes with tag = 1 (ie closest nodes to interface)
    domain->computeDistanceCloseNodes(idim,Tag,X,*lsgrad,Phi,Psi,copy);
    // second, the other nodes
    // we proceed layer by layer, going from one layer to the other
    // when a layer is converged.
    double res,resn,resnm1;
    double eps = conv_eps;
    int it=0;
    for (level=2; level<bandlevel; level++){
      com->fprintf(stdout, "*** level = %d\n", level);
      res = 1.0;
      resn = 1.0; resnm1 = 1.0;
      it = 0;
      while(res>eps || it<=1){
        resnm1 = resn;
        if(level==2 && !copylv2)
          domain->computeDistanceLevelNodes(idim,Tag,level,X,Psi,resn,Phi,MultiFluidData::FALSE);
        else
          domain->computeDistanceLevelNodes(idim,Tag,level,X,Psi,resn,Phi,copy);
        res = fabs((resn-resnm1)/(resn+resnm1));
        it++;
      }
    }

    domain->getSignedDistance(idim,Psi,Phi);
    //com->fprintf(stdout, "after reinitialization of phi[%d] to distance (phi[%d].norm = %e)\n", idim, idim, Psi.norm());

    // set Phi to the new distance function
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub){
      double (*phi)[dimLS] = Phi.subData(iSub);
      double (*psi)[1]     = Psi.subData(iSub);
      for (int i=0; i<Phi.subSize(iSub); i++){
        if(diff)
          phi[i][idim] -= psi[i][0];
        else
          phi[i][idim] = psi[i][0];
      }
    }

  }
}
//---------------------------------------------------------------------------------------------------------
