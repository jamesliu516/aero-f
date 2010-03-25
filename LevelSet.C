#include <LevelSet.h>

#include <math.h>

#include "IoData.h"
#include "Domain.h"
#include "TimeData.h"

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
    fprintf(stdout, "***Warning: if reinitialization in band --> problem!\n ***         You need to reinitialize in the whole domain with this method\n");
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
                     DistSVec<double,dimLS> &Phi, IoData &iod)
{

  map<int, FluidModelData *>::iterator it = iod.eqs.fluidModelMap.dataMap.find(1);
  if(it == iod.eqs.fluidModelMap.dataMap.end()){
    fprintf(stderr, "*** Error: no FluidModel[1] was specified\n");
    exit(1);
  }
  if(iod.eqs.fluidModel.fluid  == FluidModelData::GAS &&
     it->second->fluid == FluidModelData::LIQUID)
    invertGasLiquid = -1.0;
  else invertGasLiquid = 1.0;

  // CONVENTION: phi >= 0 --- fluid1 --- 'normal' volID
  //             phi <  0 --- fluid2 --- 'special' volID
  // Initialization of Phi is done through the use of volumeID and
  // through the knowledge of a geometric shape (with its position).

  Phi0 = -1.0;
  setupPhiVolumesInitialConditions(iod, Phi0);
  setupPhiMultiFluidInitialConditions(iod,X, Phi0);
  primitiveToConservative(Phi0, Phi, U);

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

//---------------------------------------------------------------------------------------------------------

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

  // phi[fluidId] is initialized with -1 and +1 giving the location of fluid fluidId
  // distance[fluidId] is the positive distance to the closest interface between fluid fluidId and fluid 0.
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {

    SVec<double,dimLS> &phi(Phi(iSub));
    SVec<double,dimLS> &distance(Distance(iSub));
    for (int i=0; i<phi.size(); i++)
      for (int idim=0; idim<dimLS; idim++)
        phi[i][idim] *= distance[i][idim];
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
template<int dim>
void LevelSet<dimLS>::reinitializeLevelSet(DistGeoState &geoState, 
				    DistSVec<double,3> &X, DistVec<double> &ctrlVol,
				    DistSVec<double,dim> &U, DistSVec<double,dimLS> &Phi)
{
  com->fprintf(stdout, "reinitializing\n");
  if(false)
    reinitializeLevelSetPDE(geoState,X,ctrlVol,U,Phi);
  else if(true)
    reinitializeLevelSetFM(geoState,X,ctrlVol,U,Phi);

}
//-------------------------------------------------------------------------
template<int dimLS>
template<int dim>
void LevelSet<dimLS>::reinitializeLevelSetPDE(DistGeoState &geoState, 
				    DistSVec<double,3> &X, DistVec<double> &ctrlVol,
				    DistSVec<double,dim> &U, DistSVec<double,dimLS> &Phi)
{

  /* solving reinitialization equation for the level set
  ** dphi/dt + sign(phi)*(abs(grad(phi))-1.0) = 0.0
  ** where phi is the 'primitive phi'
  **
  ** numerically, since we need to iterate over fictitious time
  ** to converge the solution, we solve
  ** dpsi/dtau + sign(phi)*(abs(grad(phi))-1.0) = 0.0
  ** and psi(tau = 0.0) = phi
  **
  ** it is not always possible to solve this equation using
  ** local time stepping, especially far from the interface.
  ** Since the distance function is needed only close to the
  ** interface, one solves this equation for nodes close to
  ** that interface (those nodes are tagged first). The
  ** remaining nodes are set to an arbitrary value. Here
  ** we chose to give them the highest value among its neighbors.
  */

  // initialize Psi
  Psi = Phi;

  // tag nodes that are close to interface up to level 'levelTot'
  int level = 0;
  for (level=0; level<bandlevel; level++)
    //domain->TagInterfaceNodes(Tag,Phi0,level);
    domain->TagInterfaceNodes(Tag,Phi,level);

  // steady state solution over tagged nodes
  //computeSteadyState(geoState, X, ctrlVol, U, Phi0); // for testing only!!
  computeSteadyState(geoState, X, ctrlVol, U, Phi);

  // psi is set to max(values of neighbours with tag>0) if tag is 0
  bool lastlevel = false;
  while(!lastlevel){
    domain->FinishReinitialization(Tag,Psi,level);
    level += 1;
    lastlevel = (Tag.min()==0 ? false : true);
  }

  //DistSVec<double,dimLS> testsign(domain->getNodeDistInfo());
  //testsign = Psi;
  //testsign *= Phi;
  //domain->checkNodePhaseChange(testsign);

  // set Phi to the new distance function
  Phi = Psi;

}

//-------------------------------------------------------------------------
template<int dimLS>
template<int dim>
void LevelSet<dimLS>::computeSteadyState(DistGeoState &geoState,
                                  DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                  DistSVec<double,dim> &U, DistSVec<double,dimLS> &Phi)
{

  bool lastIt = false;
  int iteration = 0;
  double res0;

  //need to compute only for tagged nodes =1 ie for nodes close to interface!
  if(typeTracking != MultiFluidData::LINEAR)
    lsgrad->compute(geoState.getConfig(),X,Psi);

  while(!lastIt){

    domain->computePsiResidual(X, *lsgrad, Phi, Psi, Tag,
		               w, dt, PsiRes, localtime, typeTracking);

    lastIt = checkConvergencePsi(iteration, res0);

    if(localtime){
      dt *= cfl_psi;
      PsiRes *= dt;
    }else{
      double dt_glob = cfl_psi*dt.min();
      PsiRes *= dt_glob;
    }

    Psi -= PsiRes;
    
    iteration++;

  }

}

//-------------------------------------------------------------------------
template<int dimLS>
template<int dim>
void LevelSet<dimLS>::reinitializeLevelSetFM(DistGeoState &geoState, 
				    DistSVec<double,3> &X, DistVec<double> &ctrlVol,
				    DistSVec<double,dim> &U, DistSVec<double,dimLS> &Phi)
{

  // initialize Psi
  Psi = 0.1;

  // tag nodes that are close to interface up to level 'bandlevel'
  int level = 0;
  for (level=0; level<bandlevel; level++)
    domain->TagInterfaceNodes(Tag,Phi,level);


  //compute distance for Psi<=0.5 and level<bandlevel
  // first, the nodes with tag = 1 (ie closest nodes to interface)
  domain->computeDistanceCloseNodes(Tag,X,*lsgrad,Phi,Psi,copy);
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
      domain->computeDistanceLevelNodes(Tag,level,X,Psi,resn,Phi,copy);
      res = fabs((resn-resnm1)/(resn+resnm1));
      it++;
    }
  }

  domain->getSignedDistance(Psi,Phi);

  // set Phi to the new distance function
  if(diff)
    Phi -= Psi;
  else
    Phi = Psi;
}
//---------------------------------------------------------------------------------------------------------
template<int dimLS>
bool LevelSet<dimLS>::checkConvergencePsi(int iteration, double &res0)
{

  double eps = 1.0e-6;
  double res = sqrt(PsiRes*PsiRes);

  if(iteration == 0)
    res0 = res;

  double target = eps*res0;
  if(res < target || iteration >= subIt || res < 1.e-12){
    com->fprintf(stdout, "*** ReinitLS: SubIt = %d (init = %e, res = %e, target = %e)\n", iteration, res0, res, target);
    return true;
  }

  return false;

}
