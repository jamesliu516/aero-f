#include <ImplicitTsDesc.h>

#include <DistTimeState.h>
#include <GeoSource.h>
#include <SpaceOperator.h>
#include <Domain.h>
#include <MatVecProd.h>
#include <KspSolver.h>
#include <NewtonSolver.h>

//------------------------------------------------------------------------------

template<int dim>
ImplicitTsDesc<dim>::ImplicitTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  TsDesc<dim>(ioData, geoSource, dom)
{

  tag = 0;

  failSafeNewton = ioData.ts.implicit.newton.failsafe;
  maxItsNewton = ioData.ts.implicit.newton.maxIts;
  epsNewton = ioData.ts.implicit.newton.eps;  

  //this->timeState = new DistTimeState<dim>(ioData, this->varFcn, this->domain, this->V);
  this->timeState = new DistTimeState<dim>(ioData, this->spaceOp, this->varFcn, this->domain, this->V);
  ns = new NewtonSolver<ImplicitTsDesc<dim> >(this);

}

//------------------------------------------------------------------------------

template<int dim>
ImplicitTsDesc<dim>::~ImplicitTsDesc()
{

  if (tag) delete tag;
  if (ns) delete ns;

}

//------------------------------------------------------------------------------

template<int dim>
int ImplicitTsDesc<dim>::solveNonLinearSystem(DistSVec<double,dim> &U)
{

  int its = ns->solve(U);

  return its;

}

//------------------------------------------------------------------------------
// this function evaluates (Aw),t + F(w,x,v)
template<int dim>
void ImplicitTsDesc<dim>::computeFunction(int it, DistSVec<double,dim> &Q, 
					  DistSVec<double,dim> &F)
{
  // XML
  //spaceOp->applyBCsToSolutionVector(Q);

  this->spaceOp->computeResidual(*this->X, *this->A, Q, F, this->timeState);

  this->timeState->add_dAW_dt(it, *this->geoState, *this->A, Q, F);

  this->spaceOp->applyBCsToResidual(Q, F);

}

//------------------------------------------------------------------------------
template<int dim>
void ImplicitTsDesc<dim>::recomputeFunction(DistSVec<double,dim> &Q,
                                            DistSVec<double,dim> &rhs)
{
  this->spaceOp->recomputeRHS(*this->X, Q, rhs);
}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
KspPrec<neq> *ImplicitTsDesc<dim>::createPreconditioner(PcData &pcdata, Domain *dom)
{

  KspPrec<neq> *_pc = 0;

  if (pcdata.type == PcData::IDENTITY)
    _pc = new IdentityPrec<neq>();
  else if (pcdata.type == PcData::JACOBI)
    _pc = new JacobiPrec<Scalar,neq>(DiagMat<Scalar,neq>::DENSE, dom);
  else if (pcdata.type == PcData::AS || 
	   pcdata.type == PcData::RAS || 
	   pcdata.type == PcData::ASH || 
	   pcdata.type == PcData::AAS)
    _pc = new IluPrec<Scalar,neq>(pcdata, dom);

  return _pc;

}

//------------------------------------------------------------------------------

template<int dim>
template<int neq>
KspSolver<DistSVec<double,neq>, MatVecProd<dim,neq>, KspPrec<neq>, Communicator> *
ImplicitTsDesc<dim>::createKrylovSolver(const DistInfo &info, KspData &kspdata, 
					MatVecProd<dim,neq> *_mvp, KspPrec<neq> *_pc, 
					Communicator *_com)
{

  KspSolver<DistSVec<double,neq>, MatVecProd<dim,neq>, 
    KspPrec<neq>, Communicator> *_ksp = 0;

  if (kspdata.type == KspData::RICHARDSON)
    _ksp = new RichardsonSolver<DistSVec<double,neq>, MatVecProd<dim,neq>, 
      KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  else if (kspdata.type == KspData::CG)
    _ksp = new CgSolver<DistSVec<double,neq>, MatVecProd<dim,neq>, 
      KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  else if (kspdata.type == KspData::GMRES)
    _ksp = new GmresSolver<DistSVec<double,neq>, MatVecProd<dim,neq>, 
      KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  else if (kspdata.type == KspData::GCR)
     _ksp = new GcrSolver<DistSVec<double,neq>, MatVecProd<dim,neq>,
       KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
                                                        
  return _ksp;

}

//------------------------------------------------------------------------------

template<int dim>
template<int neq>
KspSolver<DistVec<double>, MatVecProd<dim,1>, KspPrec<neq,double>, Communicator, double> *
ImplicitTsDesc<dim>::createKrylovSolverLS(const DistInfo &info, KspData &kspdata,
                                        MatVecProd<dim,1> *_mvpLS, KspPrec<neq,double> *_pc,
                                        Communicator *_com)
{
                                                                                                      
  KspSolver<DistVec<double>, MatVecProd<dim,1>,
    KspPrec<neq>, Communicator> *_ksp = 0;

  if (kspdata.type == KspData::RICHARDSON)
    _ksp = new RichardsonSolver<DistVec<double>, MatVecProd<dim,1>,
                           KspPrec<neq>, Communicator>(info, kspdata, _mvpLS, _pc, _com);
  else if (kspdata.type == KspData::CG)
    _ksp = new CgSolver<DistVec<double>, MatVecProd<dim,1>,
                           KspPrec<neq>, Communicator>(info, kspdata, _mvpLS, _pc, _com);
  else if (kspdata.type == KspData::GMRES)
    _ksp = new GmresSolver<DistVec<double>, MatVecProd<dim,1>,
                           KspPrec<neq>, Communicator>(info, kspdata, _mvpLS, _pc, _com);
  else if (kspdata.type == KspData::GCR)
    _ksp = new GcrSolver<DistVec<double>, MatVecProd<dim,1>,
                           KspPrec<neq>, Communicator>(info, kspdata, _mvpLS, _pc, _com);                                                                                                                                                                 
  return _ksp;
                                                                                                      
}

//------------------------------------------------------------------------------
template<int dim>
int ImplicitTsDesc<dim>::checkFailSafe(DistSVec<double,dim>& U)
{

  if (!failSafeNewton) return 0;

  if (!tag)
    tag = new DistSVec<bool,2>(this->getVecInfo());

  this->domain->checkFailSafe(this->varFcn, U, *tag);
  this->spaceOp->fix(*tag);

  return failSafeNewton;

}

//------------------------------------------------------------------------------
template<int dim>
void ImplicitTsDesc<dim>::resetFixesTag()
{

  this->spaceOp->resetTag();

}

//------------------------------------------------------------------------------
template<int dim>
double ImplicitTsDesc<dim>::reinitLS(DistVec<double> &Phi, DistSVec<double,dim> &U, int iti)
{
  return this->spaceOp->reinitLS(*this->X, Phi, U, iti);
}
