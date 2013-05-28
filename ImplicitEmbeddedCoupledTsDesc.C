#include "IoData.h"
#include "GeoSource.h"
#include "Domain.h"
#include "LevelSet.h"
#include "DistTimeState.h"

#include <MatVecProd.h>
#include <KspSolver.h>
#include <SpaceOperator.h>
#include <NewtonSolver.h>


#ifdef TYPE_PREC
#define PrecScalar TYPE_PREC
#else
#define PrecScalar double
#endif

//------------------------------------------------------------------------------

template<int dim>
ImplicitEmbeddedCoupledTsDesc<dim>::
ImplicitEmbeddedCoupledTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  ImplicitEmbeddedTsDesc<dim>(ioData, geoSource, dom)
{
  ImplicitData &implicitData = ioData.ts.implicit;
  
  // MatVecProd, Prec and Krylov solver for Euler equations
  if (implicitData.mvp == ImplicitData::FD)
    mvp = new MatVecProdFD<dim,dim>(implicitData,this->timeState, this->geoState,
                                      this->spaceOp,this->domain,ioData);
  else if (implicitData.mvp == ImplicitData::H1)
    mvp = new MatVecProdH1<dim,double,dim>(this->timeState, this->spaceOp, this->domain,ioData);
  else{
    this->com->fprintf(stdout, "*** Error: MatVecProdH2 is not available\n");
    exit(1);
  }
  
  pc = ImplicitEmbeddedTsDesc<dim>::template 
    createPreconditioner<dim>(implicitData.newton.ksp.ns.pc, this->domain);
  
  ksp = this->createKrylovSolver(this->getVecInfo(), implicitData.newton.ksp.ns, mvp, pc, this->com);
  
  typename MatVecProd<dim,dim>::_fsi fsi = {

    this->distLSS,
    &this->nodeTag,
    this->riemann,
    this->linRecAtInterface,
    this->viscSecOrder,
    this->Nsbar,
    &this->Wtemp,
    this->riemannNormal,
    this->ghostPoints,
  };

  mvp->AttachStructure(fsi);
  
  if (this->modifiedGhidaglia)
    mvp->attachHH(this->embeddedU);
  
}

//------------------------------------------------------------------------------

template<int dim>
ImplicitEmbeddedCoupledTsDesc<dim>::~ImplicitEmbeddedCoupledTsDesc()
{
  if (mvp)   delete mvp;
  if (pc)    delete pc;
  if (ksp)   delete ksp;

}

//------------------------------------------------------------------------------
template<int dim>
void ImplicitEmbeddedCoupledTsDesc<dim>::computeJacobian(int it, DistSVec<double,dim> &Q,
							DistSVec<double,dim> &F)
{

  MatVecProdH1<dim,double,dim> *mvph1 = dynamic_cast<MatVecProdH1<dim,double,dim> *>(mvp);
  if (mvph1)  {
    mvph1->clearGhost(); 
  }

  mvp->evaluate(it,*(this->X) ,*(this->A), Q, F);

  if (this->modifiedGhidaglia)
    mvp->evaluateHH(*this->hhResidual, *this->bcData->getBoundaryStateHH());

  mvph1 = dynamic_cast<MatVecProdH1<dim,double,dim> *>(mvp);
//  if (mvph1 && this->ghostPoints) 
//    this->domain->populateGhostJacobian(*this->ghostPoints,Q, this->varFcn, *this->distLSS, this->nodeTag,*mvph1);

}
//------------------------------------------------------------------------------
template<int dim>
void ImplicitEmbeddedCoupledTsDesc<dim>::setOperators(DistSVec<double,dim> &Q)
{
  
  DistMat<double,dim> *_pc = dynamic_cast<DistMat<double,dim> *>(pc);
  
  if (_pc) {
    
    MatVecProdFD<dim,dim> *mvpfd = dynamic_cast<MatVecProdFD<dim,dim> *>(mvp);
    MatVecProdH1<dim,double,dim> *mvph1 = dynamic_cast<MatVecProdH1<dim,double,dim> *>(mvp);
    
    if (mvpfd) {

      this->spaceOp->computeJacobian(*this->X, *this->A, Q,this->distLSS, this->nodeTag, this->riemann,
                                     this->riemannNormal, this->Nsbar,this->ghostPoints, *_pc,this->timeState);
      this->timeState->addToJacobian(*this->A, *_pc, Q);
      this->spaceOp->applyBCsToJacobian(Q, *_pc);
    }
    else if (mvph1) {
      JacobiPrec<double,dim> *jac = dynamic_cast<JacobiPrec<double,dim> *>(pc);
      IluPrec<double,dim> *ilu = dynamic_cast<IluPrec<double,dim> *>(pc);
      
      if (jac)
	jac->getData(*mvph1);
      else if (ilu)
	ilu->getData(*mvph1);
    }
    
  }
  
  double t0 = this->timer->getTime();
  
  pc->setup();
  
  double t = this->timer->addPrecSetupTime(t0);
  
  this->com->printf(6, "Fluid preconditioner computation: %f s\n", t);
  
}

//------------------------------------------------------------------------------
template<int dim>
int ImplicitEmbeddedCoupledTsDesc<dim>::solveLinearSystem(int it, DistSVec<double,dim> &b,
				                   DistSVec<double,dim> &dQ)
{
 
  double t0 = this->timer->getTime();  

  this->embeddeddQ = 0.0;
  this->embeddedB.ghost() = 0.0;
  this->embeddedB.real() = b;
  
  if (this->modifiedGhidaglia)
    this->embeddedB.hh() = -1.0*(*this->hhResidual);
  
  ksp->setup(it, this->maxItsNewton, this->embeddedB);
  
  int lits = ksp->solve(this->embeddedB, this->embeddeddQ);

  if(this->data->checklinsolve && lits==ksp->maxits) this->data->badlinsolve=true;
 
  dQ = this->embeddeddQ.real();
  this->embeddedU.ghost() += this->embeddeddQ.ghost();
  if (this->modifiedGhidaglia) {
    this->embeddedU.hh() += this->embeddeddQ.hh();

    *this->bcData->getBoundaryStateHH() = this->embeddedU.hh();
  }

  this->timer->addKspTime(t0);
  
  return lits;
  
}

