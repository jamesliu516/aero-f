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
        //phi[i] = xb*sin(zb*(x[i][0]-0.50005))+yb;
        //phi[i] = fabs(x[i][0]-xb) - r;
			}
      phi[i] *= u[i][0];
      //fprintf(stdout, "%e %e %e %e\n", x[i][0],dist, phi[i], u[i][0]);
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
      fprintf(stdout, "setup LS: data->unse_nm1\n");
      DistSVec<double,1> ReadPhi1(domain->getNodeDistInfo());
      data->exist_nm1 = domain->readVectorFromFile(name, 1, 0, ReadPhi1);
      DistVec<double> PhiRead1(domain->getNodeDistInfo(), reinterpret_cast<double (*)>(ReadPhi1.data()));
      Phinm1 = PhiRead1;
    }

    if (data->use_nm2){
      fprintf(stdout, "setup LS: data->unse_nm2\n");
      DistSVec<double,1> ReadPhi2(domain->getNodeDistInfo());
      data->exist_nm2 = domain->readVectorFromFile(name, 2, 0, ReadPhi2);
      DistVec<double> PhiRead2(domain->getNodeDistInfo(), reinterpret_cast<double (*)>(ReadPhi2.data()));
      Phinm2 = PhiRead2;
    }
  }

  if (data->use_nm1 && !data->exist_nm1){
    fprintf(stdout, "setup LS: data->use_nm1 && !data->exist_nm1\n");
    Phinm1 = Phin;
  }
  if (data->use_nm2 && !data->exist_nm2){
    fprintf(stdout, "setup LS: data->use_nm2 && !data->exist_nm2\n");
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

  fprintf(stdout, "reinitializing LevelSet\n");
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

  //lsgrad->compute(geoState.getConfig(), X, ctrlVol, Psi);
  //lsgrad->compute(X,Psi);
  //need to compute only for tagged nodes =1 ie for nodes close to interface!

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

