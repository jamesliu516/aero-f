#include <ImplicitCoupledTsDesc.h>
#include <Domain.h>
#include <GeoSource.h>
#include <DistTimeState.h>
#include <MatVecProd.h>
#include <KspSolver.h>
#include <MemoryPool.h>
//#include <MultiGridPrec.h>

#ifdef TYPE_MAT
#define MatScalar TYPE_MAT
#else
#define MatScalar double
#endif

#ifdef TYPE_PREC
#define PrecScalar TYPE_PREC
#else
#define PrecScalar double
#endif

//------------------------------------------------------------------------------

template<int dim>
ImplicitCoupledTsDesc<dim>::
ImplicitCoupledTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitTsDesc<dim>(ioData, geoSource, dom)
{

  ImplicitData &implicitData = ioData.ts.implicit;

  if (implicitData.mvp == ImplicitData::FD || implicitData.mvp == ImplicitData::H1FD)
  {
    this->mvp = new MatVecProdFD<dim,dim>(implicitData, this->timeState, this->geoState, this->spaceOp, this->domain, ioData);
  }
  else if (implicitData.mvp == ImplicitData::H1)
  {
    mvp = new MatVecProdH1<dim,MatScalar,dim>(this->timeState, this->spaceOp, this->domain, ioData);
  }
  else if (implicitData.mvp == ImplicitData::H2)
  {
    mvp = new MatVecProdH2<dim,MatScalar,dim>(ioData, this->varFcn, this->timeState, this->spaceOp, this->domain, this->geoState);
  }

#ifdef MVP_CHECK
  ImplicitData fddata;
  fddata.mvp = ImplicitData::FD;
  mvpfd1 = new MatVecProdFD<dim,dim>(fddata, this->timeState, this->geoState, this->spaceOp, this->domain, ioData);
#endif

  pc = ImplicitTsDesc<dim>::template 
    createPreconditioner<PrecScalar,dim>(implicitData.newton.ksp.ns.pc, this->domain);

  ksp = this->createKrylovSolver(this->getVecInfo(), implicitData.newton.ksp.ns, mvp, pc, this->com);

  MemoryPool mp;

  mvp->exportMemory(&mp);
  pc->exportMemory(&mp);

  this->mmh = this->createMeshMotionHandler(ioData, geoSource, &mp);

}

//------------------------------------------------------------------------------

template<int dim>
ImplicitCoupledTsDesc<dim>::~ImplicitCoupledTsDesc()
{

  if (mvp) delete mvp;
  if (pc) delete pc;
  if (ksp) delete ksp;

#ifdef MVP_CHECK
  if (mvpfd1) delete mvpfd1;
#endif

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitCoupledTsDesc<dim>::computeJacobian(int it, DistSVec<double,dim> &Q,
						 DistSVec<double,dim> &F)
{

  if(this->wallRecType==BcsWallData::CONSTANT) {
    this->com->printf(6,"ImplicitCoupledTsDesc<dim>::computeJacobian 1\n");
    mvp->evaluate(it, *this->X, *this->A, Q, F);
  } else {
    this->com->printf(6,"ImplicitCoupledTsDesc<dim>::computeJacobian 2\n");
    mvp->evaluate(*this->riemann1, it, *this->X, *this->A, Q, F);
  }
#ifdef MVP_CHECK
  DistSVec<double,dim> p(this->getVecInfo());
  DistSVec<double,dim> prod(this->getVecInfo());

  p = 1.e-2;
  mvp->apply(p, prod);
  this->com->printf(6,"ImplicitCoupledTsDesc<dim>::computeJacobian 3\n");
  this->domain->checkMatVecProd(prod, "mvp");
  this->com->printf(6,"ImplicitCoupledTsDesc<dim>::computeJacobian 4\n");

//  ImplicitData fddata;
//  fddata.mvp = ImplicitData::FD;
  //fddata.type = ImplicitData::CRANK_NICOLSON;
//  MatVecProd<dim,dim> *mvpfd = new MatVecProdFD<dim>(fddata, timeState, geoState, spaceOp, this->domain);

  computeFunction(it, Q, F);
  this->com->printf(6,"ImplicitCoupledTsDesc<dim>::computeJacobian 5\n");
  mvpfd1->evaluate(it, *this->X, *this->A, Q, F);
  this->com->printf(6,"ImplicitCoupledTsDesc<dim>::computeJacobian 6\n");
  mvpfd1->apply(p, prod);
  this->com->printf(6,"ImplicitCoupledTsDesc<dim>::computeJacobian 7\n");
  this->domain->checkMatVecProd(prod, "mvpfd1");
  this->com->printf(6,"ImplicitCoupledTsDesc<dim>::computeJacobian 8\n");

  this->com->barrier();
  exit(1);
#endif

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitCoupledTsDesc<dim>::setOperators(DistSVec<double,dim> &Q)
{

  double t0 = this->timer->getTime();
  
  DistMat<PrecScalar,dim> *_pc = dynamic_cast<DistMat<PrecScalar,dim> *>(pc);
  DistMat<double,dim> *_pc2 = dynamic_cast<DistMat<double,dim> *>(pc);
  MultiGridPrec<PrecScalar,dim> *pmg = dynamic_cast<MultiGridPrec<PrecScalar,dim> *>(pc);

  if (_pc || _pc2) {

    MatVecProdFD<dim, dim> *mvpfd = dynamic_cast<MatVecProdFD<dim, dim> *>(mvp);
    MatVecProdH1<dim,MatScalar,dim> *mvph1 = dynamic_cast<MatVecProdH1<dim,MatScalar,dim> *>(mvp);
    MatVecProdH2<dim,MatScalar,dim> *mvph2 = dynamic_cast<MatVecProdH2<dim,MatScalar,dim> *>(mvp);

    if (mvpfd || mvph2) {
      if (_pc) {

        this->com->printf(6,"ImplicitCoupledTsDesc<dim>::setOperators 1\n");     
        this->spaceOp->computeJacobian(*this->X, *this->A, Q, *_pc, this->timeState);
        this->com->printf(6,"ImplicitCoupledTsDesc<dim>::setOperators 2\n");     
        this->timeState->addToJacobian(*this->A, *_pc, Q);
        this->com->printf(6,"ImplicitCoupledTsDesc<dim>::setOperators 3\n");     
        this->spaceOp->applyBCsToJacobian(Q, *_pc);
        this->com->printf(6,"ImplicitCoupledTsDesc<dim>::setOperators 4\n");     
      } else {
        this->com->printf(6,"ImplicitCoupledTsDesc<dim>::setOperators 5\n");     
        this->spaceOp->computeJacobian(*this->X, *this->A, Q, *_pc2, this->timeState);
        this->com->printf(6,"ImplicitCoupledTsDesc<dim>::setOperators 6\n");     
        this->timeState->addToJacobian(*this->A, *_pc2, Q);
        this->com->printf(6,"ImplicitCoupledTsDesc<dim>::setOperators 7\n");     
        this->spaceOp->applyBCsToJacobian(Q, *_pc2);
        this->com->printf(6,"ImplicitCoupledTsDesc<dim>::setOperators 8\n");     
        if (pmg) {
          if (!pmg->isInitialized())
            pmg->initialize();
          this->com->printf(6,"ImplicitCoupledTsDesc<dim>::setOperators 9\n");     
          pmg->getData(*_pc2);
          this->com->printf(6,"ImplicitCoupledTsDesc<dim>::setOperators 10\n");     
        }
      }
    }
    else if (mvph1) {
      JacobiPrec<PrecScalar,dim> *jac = dynamic_cast<JacobiPrec<PrecScalar,dim> *>(pc);
      IluPrec<PrecScalar,dim> *ilu = dynamic_cast<IluPrec<PrecScalar,dim> *>(pc);
      MultiGridPrec<PrecScalar,dim> *pmg = dynamic_cast<MultiGridPrec<PrecScalar,dim> *>(pc);
      
      if (jac) {
        this->com->printf(6,"ImplicitCoupledTsDesc<dim>::setOperators 11\n");     
        jac->getData(*mvph1);
      } else if (ilu) {
        this->com->printf(6,"ImplicitCoupledTsDesc<dim>::setOperators 12\n");     
        ilu->getData(*mvph1);
      } else if (pmg) {
        this->com->printf(6,"ImplicitCoupledTsDesc<dim>::setOperators 13\n");     
        if (!pmg->isInitialized())
          pmg->initialize();
        pmg->getData(*mvph1);
      }

    }

  }
  this->com->printf(6,"ImplicitCoupledTsDesc<dim>::setOperators 14\n");     
    
  
  //if (!pmg)    
    pc->setup();
  this->com->printf(6,"ImplicitCoupledTsDesc<dim>::setOperators 15\n");     
  //else
  //  pmg->setup(Q);
/*  
  MultiGridPrec<PrecScalar,dim> *pmg = dynamic_cast<MultiGridPrec<PrecScalar,dim> *>(pc);
  MatVecProdH1<dim,MatScalar,dim> *mvph1 = dynamic_cast<MatVecProdH1<dim,MatScalar,dim> *>(mvp);
  if (pmg && mvph1) { 
    //this->spaceOp->conservativeToPrimitive(Q);
    pmg->getData(*mvph1, *this->spaceOp->getCurrentPrimitiveVector(),
                 *this->spaceOp, this->timeState); 
  }
  //else
  //  fprintf(stderr,"Tried to getData but instead got %x, %x\n", pmg, mvph1);
  */

  double t = this->timer->addPrecSetupTime(t0);

  this->com->printf(6, "Fluid preconditioner computation: %f s\n", t);

}

//------------------------------------------------------------------------------

template<int dim>
int ImplicitCoupledTsDesc<dim>::solveLinearSystem(int it, DistSVec<double,dim> &b, 
						  DistSVec<double,dim> &dQ)
{

  double t0 = this->timer->getTime();

  dQ = 0.0;

  this->com->printf(6,"ImplicitCoupledTsDesc<dim>::solveLinearSystem 1\n");
  ksp->setup(it, this->maxItsNewton, b);
  this->com->printf(6,"ImplicitCoupledTsDesc<dim>::solveLinearSystem 2\n");

  int lits = ksp->solve(b, dQ);
  this->com->printf(6,"ImplicitCoupledTsDesc<dim>::solveLinearSystem 3\n");

  if (lits == ksp->maxits && this->data->checklinsolve) this->data->badlinsolve=true;
  this->com->printf(6,"ImplicitCoupledTsDesc<dim>::solveLinearSystem 4\n");

  this->timer->addKspTime(t0);
  this->com->printf(6,"ImplicitCoupledTsDesc<dim>::solveLinearSystem 5\n");

  return lits;

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void ImplicitCoupledTsDesc<dim>::rstVarImplicitCoupledTsDesc(IoData &ioData)
{

#ifdef MVP_CHECK
    mvpfd1->rstSpaceOp(ioData, this->varFcn, this->spaceOp, false);
#endif

    mvp->rstSpaceOp(ioData, this->varFcn, this->spaceOp, false);

}

//------------------------------------------------------------------------------
