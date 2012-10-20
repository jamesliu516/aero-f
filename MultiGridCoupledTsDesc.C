#include <VectorSet.h>
#include <LevelSet/LevelSetStructure.h>
#include <MultiGridCoupledTsDesc.h>

template <int dim>
MultiGridCoupledTsDesc<dim>::
MultiGridCoupledTsDesc(IoData & iod, GeoSource & gs,  Domain * dom) :
  ImplicitCoupledTsDesc<dim>(iod,gs,dom) {

  memset(numSmooths_post,0,sizeof(numSmooths_post));
  numSmooths_pre[0] = 1;
  numSmooths_pre[1] = 2;
  numSmooths_pre[2] = 3;
  numSmooths_pre[3] = 3;
  numSmooths_pre[4] = 3;
  numSmooths_pre[5] = 3;

  smoothWithGMRES = (iod.mg.mg_smoother == MultiGridData::MGGMRES);

  prolong_relax_factor = iod.mg.prolong_relax_factor;
  restrict_relax_factor = iod.mg.restrict_relax_factor;

  if (iod.mg.cycle_scheme == MultiGridData::VCYCLE)
    mc = 1;
  else if (iod.mg.cycle_scheme == MultiGridData::WCYCLE)
    mc = 2;     

  globalIt = 0;
}

template <int dim>
MultiGridCoupledTsDesc<dim>::
~MultiGridCoupledTsDesc() {

}


template <int dim>
void MultiGridCoupledTsDesc<dim>::
setupTimeStepping(DistSVec<double,dim> *U0, IoData &iod) {

  ImplicitCoupledTsDesc<dim>::setupTimeStepping(U0,iod);

  pKernel = new MultiGridKernel<double>(this->domain, *this->geoState,
                                        iod,iod.mg.num_multigrid_levels);

  pKernel->initialize(dim,dim,0);

  mgSpaceOp = 
    new MultiGridSpaceOperator<double,dim>(iod, this->domain, this->spaceOp, pKernel);

  mgMvp = new MultiGridMvpMatrix<double,dim>(this->domain,pKernel);

  mgKspSolver = new MultiGridKspSolver<double,dim,double>(this->domain, iod.ts.implicit.newton.ksp.ns,
                                                          pKernel);

  if (!smoothWithGMRES) {
    typename MultiGridSmoothingMatrix<double,dim>::SmoothingMode s;
    if (iod.mg.mg_smoother == MultiGridData::MGRAS)
      s = MultiGridSmoothingMatrix<double,dim>::RAS; 
    if (iod.mg.mg_smoother == MultiGridData::MGJACOBI)
      s = MultiGridSmoothingMatrix<double,dim>::BlockJacobi;
    if (iod.mg.mg_smoother == MultiGridData::MGLINEJACOBI)
      s = MultiGridSmoothingMatrix<double,dim>::LineJacobi;
    smoothingMatrices = new MultiGridSmoothingMatrices<double,dim>(pKernel, s);

  }

  V.init(pKernel);
  res.init(pKernel);
  R.init(pKernel);
  F.init(pKernel);
  Forig.init(pKernel);
  U.init(pKernel);
  Uold.init(pKernel); 
  dx.init(pKernel);
}

template <int dim>
void MultiGridCoupledTsDesc<dim>::
smooth0(DistSVec<double,dim>& x,int steps) {

  int i;
  double dummy = 0.0;
  updateStateVectors(x, 0);
  for (i = 0; i < steps; ++i) {

    computeTimeStep(globalIt,&dummy, x);
    ++globalIt;
    computeFunction(0, x, R(0));
    computeJacobian(0, x, R(0));
    if (smoothWithGMRES)
      setOperators(x);
    else
      smoothingMatrices->acquire( *this->GetJacobian());
  
    R(0) *= -1.0;
    if (smoothWithGMRES)
      solveLinearSystem(0, R(0),dx(0));
    else
      smoothingMatrices->apply(0, dx, R);
    
    x += dx(0);
    
    updateStateVectors(x, 0);
    monitorConvergence(0, x);
    R(0) = -1.0*this->getCurrentResidual();
    
  }
  if (i == 0) {
    monitorConvergence(0, x);
    R(0) = -1.0*this->getCurrentResidual();
  }
}

template <int dim>
void MultiGridCoupledTsDesc<dim>::
smooth(int lvl, MultiGridDistSVec<double,dim>& x,
       DistSVec<double,dim>& f,int steps) {

  int i;
  for (i = 0; i < steps; ++i) {

    this->varFcn->conservativeToPrimitive(x(lvl), V(lvl));

    mgSpaceOp->updateStateVectors(lvl,x);

    mgSpaceOp->computeTimeStep(lvl,this->data->cfl*pow(0.5,lvl),
                               V);
 
    mgSpaceOp->computeResidual(lvl, x, V, res);
    mgSpaceOp->computeJacobian(lvl, x, V, *mgMvp);
    R(lvl) = f-res(lvl);
    if (smoothWithGMRES)
      mgKspSolver->solve(lvl, *mgMvp, R, dx);  
    else {
      smoothingMatrices->acquire(lvl, *mgMvp);
      smoothingMatrices->apply(lvl, dx, R); 
    }
    
    x(lvl) += dx(lvl);
   
    this->varFcn->conservativeToPrimitive(x(lvl), V(lvl));
    pKernel->fixNegativeValues(lvl,V(lvl), x(lvl), dx(lvl), f,Forig(lvl), this->varFcn);
    mgSpaceOp->computeResidual(lvl, x, V, R, false);
    R(lvl) = f-R(lvl);
  }
}

template <int dim>
void MultiGridCoupledTsDesc<dim>::cycle(int lvl, DistSVec<double,dim>& f,
                                        MultiGridDistSVec<double,dim>& x) {

  if (lvl == 0) { 
    smooth0(x(lvl), numSmooths_pre[0]);
    mgSpaceOp->setupBcs(this->getSpaceOperator()->getDistBcData());
  }
  else
    smooth(lvl,x, f,  numSmooths_pre[lvl]);

  if (lvl < pKernel->numLevels()-1) {

    pKernel->Restrict(lvl+1, x(lvl), U(lvl+1));
    pKernel->Restrict(lvl+1, R(lvl), R(lvl+1));
    Uold(lvl+1) = U(lvl+1);
    
    this->varFcn->conservativeToPrimitive(U(lvl+1), V(lvl+1));
    mgSpaceOp->computeResidual(lvl+1, U, V, F, false);
    pKernel->applyFixes(lvl+1, R(lvl+1));
    F(lvl+1) += R(lvl+1)*restrict_relax_factor;
    for (int i = 0; i < mc; ++i)
      cycle(lvl+1, F(lvl+1), U);
    
    pKernel->Prolong(lvl+1, Uold(lvl+1), U(lvl+1), x(lvl), prolong_relax_factor);
  }
  
  if (lvl == 0) 
    smooth0(x(lvl), numSmooths_post[0]);
  else
    smooth(lvl,x, f,  numSmooths_post[lvl]);
}

template <int dim>
void MultiGridCoupledTsDesc<dim>::cycle(DistSVec<double,dim>& x) {

  F(0) = 0.0;

  U(0) = x;

  cycle(0, F(0), U);  

  x = U(0);
}

template class MultiGridCoupledTsDesc<5>;
template class MultiGridCoupledTsDesc<6>;
template class MultiGridCoupledTsDesc<7>;

