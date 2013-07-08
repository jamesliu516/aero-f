#include <VectorSet.h>
#include <LevelSet/LevelSetStructure.h>
#include <MultiGridSegTsDesc.h>

template <int dim,int neq1,int neq2>
MultiGridSegTsDesc<dim,neq1,neq2>::
MultiGridSegTsDesc(IoData & iod, GeoSource & gs,  Domain * dom) :
  ImplicitSegTsDesc<dim,neq1,neq2>(iod,gs,dom) {

  memset(numSmooths_post,0,sizeof(numSmooths_post));
  numSmooths_pre[0] = 1;
  numSmooths_pre[1] = 2;
  numSmooths_pre[2] = 3;
  numSmooths_pre[3] = 3;
  numSmooths_pre[4] = 3;
  numSmooths_pre[5] = 3;

 
  prolong_relax_factor = iod.mg.prolong_relax_factor;
  restrict_relax_factor = iod.mg.restrict_relax_factor;

  if (iod.mg.cycle_scheme == MultiGridData::VCYCLE)
    mc = 1;
  else if (iod.mg.cycle_scheme == MultiGridData::WCYCLE)
    mc = 2;     

  smoothWithGMRES = (iod.mg.mg_smoother == MultiGridData::MGGMRES);

  globalIt = 0;

  mgMvp1 = NULL;
  pKernel = NULL;
  mgSpaceOp = NULL;
  mgKspSolver = NULL;
  smoothingMatrices1 = NULL;
  smoothingMatrices2 = NULL;
}

template <int dim,int neq1,int neq2>
MultiGridSegTsDesc<dim,neq1,neq2>::
~MultiGridSegTsDesc() {

  if (mgMvp1)
    delete mgMvp1;
  if (pKernel)
    delete pKernel;
  if (mgSpaceOp)
    delete mgSpaceOp;
  if (mgKspSolver)
    delete mgKspSolver;
  if (smoothingMatrices1)
    delete smoothingMatrices1;
  if (smoothingMatrices2)
    delete smoothingMatrices2;
}


template <int dim,int neq1,int neq2>
void MultiGridSegTsDesc<dim,neq1,neq2>::
setupTimeStepping(DistSVec<double,dim> *U0, IoData &iod) {

  ImplicitSegTsDesc<dim,neq1,neq2>::setupTimeStepping(U0,iod);

  pKernel = new MultiGridKernel<double>(this->domain, *this->geoState,
                                        iod,iod.mg.num_multigrid_levels);

  pKernel->initialize(dim,neq1,neq2);

  mgSpaceOp = 
    new MultiGridSpaceOperator<double,dim>(iod, this->domain, this->spaceOp, pKernel,
                                           this->spaceOp1);

  mgMvp1 = new MultiGridMvpMatrix<double,neq1>(this->domain,pKernel);

  mgKspSolver = new MultiGridKspSolver<double,neq1,double>(this->domain, iod.ts.implicit.newton.ksp.ns,
                                                          pKernel);

  pKernel->setUseVolumeWeightedAverage(iod.mg.restrictMethod == MultiGridData::VOLUME_WEIGHTED);

  if (!smoothWithGMRES) {
    typename MultiGridSmoothingMatrix<double,neq1>::SmoothingMode s;
    if (iod.mg.mg_smoother == MultiGridData::MGRAS)
      s = MultiGridSmoothingMatrix<double,neq1>::RAS; 
    if (iod.mg.mg_smoother == MultiGridData::MGJACOBI)
      s = MultiGridSmoothingMatrix<double,neq1>::BlockJacobi;
    if (iod.mg.mg_smoother == MultiGridData::MGLINEJACOBI)
      s = MultiGridSmoothingMatrix<double,neq1>::LineJacobi;
    smoothingMatrices1 = new MultiGridSmoothingMatrices<double,neq1>(pKernel, s);
    
    typename MultiGridSmoothingMatrix<double,neq2>::SmoothingMode s2;
    if (iod.mg.mg_smoother == MultiGridData::MGRAS)
      s2 = MultiGridSmoothingMatrix<double,neq2>::RAS; 
    if (iod.mg.mg_smoother == MultiGridData::MGJACOBI)
      s2 = MultiGridSmoothingMatrix<double,neq2>::BlockJacobi;
    if (iod.mg.mg_smoother == MultiGridData::MGLINEJACOBI)
      s2 = MultiGridSmoothingMatrix<double,neq2>::LineJacobi;
    smoothingMatrices2 = new MultiGridSmoothingMatrices<double,neq2>(pKernel,s2);

  }

  V.init(pKernel);
  res.init(pKernel);
  R.init(pKernel);
  R1.init(pKernel);
  R2.init(pKernel);
  F.init(pKernel);
  Forig.init(pKernel);
  U.init(pKernel);
  Uold.init(pKernel); 
  dx.init(pKernel);
  dx1.init(pKernel);
  dx2.init(pKernel);

  update_tmp.init(pKernel);  
}

template <int dim,int neq1,int neq2>
void MultiGridSegTsDesc<dim,neq1,neq2>::
smooth0(DistSVec<double,dim>& x,int steps) {

  int i;
  double dummy = 0.0;
  this->updateStateVectors(x, 0);
  for (i = 0; i < steps; ++i) {

    this->computeTimeStep(globalIt,&dummy, x);
    ++globalIt;
    this->computeFunction(0, x, R(0));
    this->computeJacobian(0, x, R(0));
    if (smoothWithGMRES)
      this->setOperators(x);
    else {
      smoothingMatrices1->acquire(*this->GetJacobian1());
      smoothingMatrices2->acquire(*this->GetJacobian2());
    }
  
    R(0) *= -1.0;
    if (smoothWithGMRES)
      this->solveLinearSystem(0, R(0),dx(0));
    else {
      R(0).split(R1(0), R2(0));
      smoothingMatrices1->apply(0, dx1, R1);
      smoothingMatrices2->apply(0, dx2, R2);
      dx(0).merge(dx1(0),dx2(0));
    }
    
    x += dx(0);
    
    this->updateStateVectors(x, 0);
    this->monitorConvergence(0, x);
    R(0) = -1.0*this->getCurrentResidual();
    
  }
  if (i == 0) {
    this->monitorConvergence(0, x);
    R(0) = -1.0*this->getCurrentResidual();
  }
  double one = 1.0;
  if (globalIt%500 == 1) 
    this->domain->writeVectorToFile("myResidual", globalIt/500, globalIt, R(0), &one);
}

template <int dim,int neq1,int neq2>
void MultiGridSegTsDesc<dim,neq1,neq2>::
smooth(int lvl, MultiGridDistSVec<double,dim>& x,
       DistSVec<double,dim>& f,int steps) {

  int i;
  for (i = 0; i < steps; ++i) {

    this->varFcn->conservativeToPrimitive(x(lvl), V(lvl));

    mgSpaceOp->updateStateVectors(lvl,x);

    mgSpaceOp->computeTimeStep(lvl,this->data->cfl*pow(0.75,lvl),
                               V);
 
    mgSpaceOp->computeResidual(lvl, x, V, res);
    mgSpaceOp->computeJacobian(lvl, x, V, *mgMvp1);
    R(lvl) = f-1.0*res(lvl);
    R(lvl).split(R1(lvl), R2(lvl));
    if (smoothWithGMRES) {
      mgKspSolver->solve(lvl, *mgMvp1, R1, dx1);
    }
    else {
      smoothingMatrices1->acquire(lvl, *mgMvp1);
      smoothingMatrices1->apply(lvl, dx1, R1); 
    }
    dx2(lvl) = 0.0;
    dx(lvl).merge(dx1(lvl), dx2(lvl));
    
    x(lvl) += dx(lvl);
   
    this->varFcn->conservativeToPrimitive(x(lvl), V(lvl));

    pKernel->fixNegativeValues(lvl,V(lvl), x(lvl), dx(lvl), f,Forig(lvl), this->varFcn,
                               mgSpaceOp->getOperator(lvl));

  }
 
  if (steps > 0) {   
    mgSpaceOp->computeResidual(lvl, x, V, R, false);
    R(lvl) = f-R(lvl);
  }
/*
  if (lvl == 1) {
    pKernel->getLevel(lvl)->writePVTUSolutionFile("r.sol",x(lvl));
    MPI_Barrier(MPI_COMM_WORLD);
    exit(0);
  }
*/
}

template <int dim,int neq1,int neq2>
void MultiGridSegTsDesc<dim,neq1,neq2>::cycle(int lvl, DistSVec<double,dim>& f,
                                    MultiGridDistSVec<double,dim>& x) {

  if (lvl == 0) { 
    smooth0(x(lvl), numSmooths_pre[0]);
    mgSpaceOp->setupBcs(this->getSpaceOperator()->getDistBcData());
  }
  else
    smooth(lvl,x, f,  numSmooths_pre[lvl]);

  if (lvl < pKernel->numLevels()-1) {

    pKernel->Restrict(lvl+1, x(lvl), U(lvl+1));
    //this->domain->getCommunicator()->fprintf(stderr,"Restricting residual...\n");
    //fflush(stderr);
    //MPI_Barrier(MPI_COMM_WORLD);
    pKernel->Restrict(lvl+1, R(lvl), R(lvl+1));
    Uold(lvl+1) = U(lvl+1);
    
    this->varFcn->conservativeToPrimitive(U(lvl+1), V(lvl+1));
    mgSpaceOp->computeResidual(lvl+1, U, V, F, false);
    pKernel->applyFixes(lvl+1, R(lvl+1));
    //pKernel->fixNegativeValues(lvl+1,V(lvl+1), U(lvl+1), dx(lvl+1), F(lvl+1), this->varFcn);
    /*if (lvl == 0 && globalIt % 25 == 0) {

      pKernel->getLevel(lvl+1)->writePVTUSolutionFile("myR",R(lvl+1));
    }*/
    Forig(lvl+1) = F(lvl+1);
    F(lvl+1) += R(lvl+1)*restrict_relax_factor;
    for (int i = 0; i < mc; ++i)
      cycle(lvl+1, F(lvl+1), U);
    
    update_tmp(lvl) = 0.0;
    pKernel->Prolong(lvl+1, Uold(lvl+1), U(lvl+1), update_tmp(lvl), prolong_relax_factor);

    pKernel->applyFixes(lvl,update_tmp(lvl));
    x(lvl) += update_tmp(lvl);

  }
  
  if (lvl == 0) 
    smooth0(x(lvl), numSmooths_post[0]);
  else
    smooth(lvl,x, f,  numSmooths_post[lvl]);
}

template <int dim,int neq1,int neq2>
void MultiGridSegTsDesc<dim,neq1,neq2>::cycle(DistSVec<double,dim>& x) {

  F(0) = 0.0;

  U(0) = x;
  
  cycle(0, F(0), U);

  x = U(0); 
}

template class MultiGridSegTsDesc<6,5,1>;
template class MultiGridSegTsDesc<7,5,2>;

