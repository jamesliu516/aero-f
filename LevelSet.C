#include <LevelSet.h>

#include <math.h>



//-------------------------------------------------------------------------
template<int dim>
void LevelSet::setup(char *name, DistSVec<double,3> &X, DistVec<double> &Phi,
		     DistSVec<double,dim> &U, IoData &iod)
{

  if(iod.eqs.fluidModel.fluid  == FluidModelData::GAS && 
     iod.eqs.fluidModel2.fluid == FluidModelData::LIQUID)
    invertGasLiquid = -1.0;
  else invertGasLiquid = 1.0;

  double dist, r, xb, yb, zb;
  xb   = iod.mf.icd.s1.cen_x;
  yb   = iod.mf.icd.s1.cen_y;
  zb   = iod.mf.icd.s1.cen_z;
  r    = iod.mf.icd.s1.r;

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub){
    double (*x)[3] = X.subData(iSub);
    double (*u)[dim] = U.subData(iSub);
    double (*phi) = Phi.subData(iSub);
    for (int i=0; i<X.subSize(iSub); i++){
      //for bubble
      if(iod.mf.problem == MultiFluidData::BUBBLE){
        phi[i] = invertGasLiquid*(sqrt( (x[i][0] -xb)*(x[i][0] -xb)  +
                                        (x[i][1] -yb)*(x[i][1] -yb)  +
                                        (x[i][2] -zb)*(x[i][2] -zb))  -r);
        //phi[i] = 2.0*sin(phi[i]/4.0)+1.3;
      }else if(iod.mf.problem == MultiFluidData::SHOCKTUBE){
      //for shock tube (comments: cf LevelSetCore.C)
          phi[i] = x[i][0] - r;
          phi[i] = (xb*x[i][0]+yb*x[i][1]+zb*x[i][2]+r)/sqrt(xb*xb+yb*yb+zb*zb);
      }
      phi[i] *= u[i][0];
    }
  }
  Phin   = Phi;
  Phinm1 = Phin;
  Phinm2 = Phinm1;

  if (name[0] != 0) {
    DistSVec<double,1> ReadPhi(domain->getNodeDistInfo());
    domain->readVectorFromFile(name, 0, 0, ReadPhi);
    DistVec<double> PhiRead(domain->getNodeDistInfo(), reinterpret_cast<double (*)>(ReadPhi.data()));
    Phi  = PhiRead;
    Phin = PhiRead;

    if (data->use_nm1){
      DistSVec<double,1> ReadPhi1(domain->getNodeDistInfo());
      data->exist_nm1 = domain->readVectorFromFile(name, 1, 0, ReadPhi1);
      DistVec<double> PhiRead1(domain->getNodeDistInfo(), reinterpret_cast<double (*)>(ReadPhi1.data()));
      Phinm1 = PhiRead1;
    }

    if (data->use_nm2){
      DistSVec<double,1> ReadPhi2(domain->getNodeDistInfo());
      data->exist_nm2 = domain->readVectorFromFile(name, 2, 0, ReadPhi2);
      DistVec<double> PhiRead2(domain->getNodeDistInfo(), reinterpret_cast<double (*)>(ReadPhi2.data()));
      Phinm2 = PhiRead2;
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
template<int dim>
void LevelSet::conservativeToPrimitive(DistVec<double> &Cons, DistVec<double> &Prim, 
	                               DistSVec<double,dim> &U)
{
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub){
    double (*u)[dim] = U.subData(iSub);
    double (*prim) = Prim.subData(iSub);
    double (*cons) = Cons.subData(iSub);
    for (int i=0; i<U.subSize(iSub); i++)
      prim[i] = cons[i]/u[i][0];
  }

}

//-------------------------------------------------------------------------
template<int dim>
void LevelSet::primitiveToConservative(DistVec<double> &Prim, DistVec<double> &Cons,
				       DistSVec<double,dim> &U)
{
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub){
    double (*u)[dim] = U.subData(iSub);
    double (*prim) = Prim.subData(iSub);
    double (*cons) = Cons.subData(iSub);
    for (int i=0; i<U.subSize(iSub); i++)
      cons[i] = prim[i]*u[i][0];
  }

}

//-------------------------------------------------------------------------
template<int dim>
void LevelSet::reinitializeLevelSet(DistGeoState &geoState, 
				    DistSVec<double,3> &X, DistVec<double> &ctrlVol,
				    DistSVec<double,dim> &U, DistVec<double> &Phi)
{
  com->fprintf(stdout, "reinitializing\n");
  if(false)
    reinitializeLevelSetPDE(geoState,X,ctrlVol,U,Phi);
  else if(true)
    reinitializeLevelSetFM(geoState,X,ctrlVol,U,Phi);

}
//-------------------------------------------------------------------------
template<int dim>
void LevelSet::reinitializeLevelSetPDE(DistGeoState &geoState, 
				    DistSVec<double,3> &X, DistVec<double> &ctrlVol,
				    DistSVec<double,dim> &U, DistVec<double> &Phi)
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

  //fprintf(stdout, "computeSteadyState done\n");
  // psi is set to max(values of neighbours with tag>0) if tag is 0
  bool lastlevel = false;
  while(!lastlevel){
    domain->FinishReinitialization(Tag,Psi,level);
    level += 1;
    lastlevel = (Tag.min()==0 ? false : true);
  }

  DistSVec<double,1> testsign(domain->getNodeDistInfo());
  testsign = Psi;
  testsign *= Phi;
  domain->checkNodePhaseChange(testsign);

  // set Phi to the new distance function
  DistVec<double> distance(domain->getNodeDistInfo(), reinterpret_cast<double (*)>(Psi.data()));
  Phi = distance;

}

//-------------------------------------------------------------------------
template<int dim>
void LevelSet::computeSteadyState(DistGeoState &geoState,
                                  DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                  DistSVec<double,dim> &U, DistVec<double> &Phi)
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
template<int dim>
void LevelSet::reinitializeLevelSetFM(DistGeoState &geoState, 
				    DistSVec<double,3> &X, DistVec<double> &ctrlVol,
				    DistSVec<double,dim> &U, DistVec<double> &Phi)
{

  // initialize Psi
  Psi = 0.1;

  // tag nodes that are close to interface up to level 'bandlevel'
	com->fprintf(stdout, "Tagging nodes\n");
  int level = 0;
  for (level=0; level<bandlevel; level++)
    domain->TagInterfaceNodes(Tag,Phi,level);


  //compute distance for Psi<=0.5 and level<bandlevel
  // first, the nodes with tag = 1 (ie closest nodes to interface)
	com->fprintf(stdout, "computing distance for close nodes to interface\n");
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
  DistVec<double> distance(domain->getNodeDistInfo(), reinterpret_cast<double (*)>(Psi.data()));
  if(diff)
    Phi -= distance;
  else
    Phi = distance;
}
