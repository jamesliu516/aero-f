#include <ImplicitLevelSetTsDesc.h>

#include <GeoSource.h>
#include <DistTimeState.h>
#include <MatVecProd.h>
#include <KspSolver.h>
#include <SpaceOperator.h>
#include <Domain.h>
#include <NewtonSolver.h>


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
ImplicitLevelSetTsDesc<dim>::
ImplicitLevelSetTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  LevelSetTsDesc<dim>(ioData, geoSource, dom)
{
  tag = 0;
  ImplicitData &implicitData = ioData.ts.implicit;
  
  // NewtonSolver
  ns = new NewtonSolver<ImplicitLevelSetTsDesc<dim> >(this);
  failSafeNewton = implicitData.newton.failsafe;
  maxItsNewton = implicitData.newton.maxIts;
  epsNewton = implicitData.newton.eps;
  epsAbsResNewton = implicitData.newton.epsAbsRes;
  epsAbsIncNewton = implicitData.newton.epsAbsInc;

  // MatVecProd, Prec and Krylov solver for Euler equations
  if (implicitData.mvp == ImplicitData::FD || implicitData.mvp == ImplicitData::H1FD)
    mvp = new MatVecProdFD<dim, dim>(implicitData, this->timeState, this->geoState, this->spaceOp, this->domain, ioData);
  else if (implicitData.mvp == ImplicitData::H1)
    mvp = new MatVecProdH1<dim,MatScalar,dim>(this->timeState, this->spaceOp, this->domain, ioData);
  else if (implicitData.mvp == ImplicitData::H2)
    mvp = new MatVecProdH2<MatScalar,dim>(ioData, this->varFcn, this->timeState, this->spaceOp, this->domain);
  
  pc = createPreconditioner<PrecScalar,dim>(implicitData.newton.ksp.ns.pc, this->domain);
  
  ksp = createKrylovSolver(this->getVecInfo(), implicitData.newton.ksp.ns, mvp, pc, this->com);
  

  // MatVecProd, Prec and Krylov solver for LevelSet equation
  if (implicitData.mvp != ImplicitData::FD)
    this->com->printf(6, "H1 and H2 not implemented for Level-Set, FD only..... hum to check!!!\n");
  mvpLS  = new MatVecProdLS<MatScalar,dim,1>(ioData, this->varFcn, this->timeState, this->geoState, this->spaceOp, this->domain, this->LS);

  pcLS = createPreconditioner<PrecScalar,1>(implicitData.newton.ksp.lsi.pc, this->domain);
  kspLS = createKrylovSolverLS(this->getVecInfo(), implicitData.newton.ksp.lsi, mvpLS, pcLS, this->com);


  // meshmotion
  MemoryPool mp;
  
  mvp->exportMemory(&mp);
  pc->exportMemory(&mp);
  
  this->mmh = this->createMeshMotionHandler(ioData, geoSource, &mp);

}

//------------------------------------------------------------------------------

template<int dim>
ImplicitLevelSetTsDesc<dim>::~ImplicitLevelSetTsDesc()
{
  if (tag)   delete tag;
  if (mvp)   delete mvp;
  if (pc)    delete pc;
  if (ksp)   delete ksp;
  if (mvpLS) delete mvpLS;
  if (kspLS) delete kspLS;
  if (pcLS)  delete pcLS;
  if (ns)    delete ns;

}

//------------------------------------------------------------------------------
//  Internal routines to setup the class (called in constructor)
//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
KspPrec<neq> *ImplicitLevelSetTsDesc<dim>::createPreconditioner(PcData &pcdata, Domain *dom)
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
ImplicitLevelSetTsDesc<dim>::createKrylovSolver(const DistInfo &info, KspData &kspdata,
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
  
  return _ksp;
  
}

//------------------------------------------------------------------------------
template<int dim>
template<int neq>
KspSolver<DistVec<double>, MatVecProd<dim,neq>, KspPrec<neq>, Communicator, double> *
ImplicitLevelSetTsDesc<dim>::createKrylovSolverLS(const DistInfo &info, KspData &kspdata,
                                        MatVecProd<dim,neq> *_mvp, KspPrec<neq,double> *_pc,
                                        Communicator *_com)
{

  KspSolver<DistVec<double>, MatVecProd<dim,neq>,
    KspPrec<neq>, Communicator> *_ksp = 0;

  if (kspdata.type == KspData::RICHARDSON)
    _ksp = new RichardsonSolver<DistVec<double>, MatVecProd<dim,neq>,
                 KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  else if (kspdata.type == KspData::CG)
    _ksp = new CgSolver<DistVec<double>, MatVecProd<dim,neq>,
                 KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  else if (kspdata.type == KspData::GMRES)
    _ksp = new GmresSolver<DistVec<double>, MatVecProd<dim,neq>,
                 KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);

  return _ksp;
}


//------------------------------------------------------------------------------
// External routine to solve problem (called by TsSolver)
// It calls for the NewtonSolver ns, which in turn will
// call routines below from this same file or from LevelSetTsDesc
//------------------------------------------------------------------------------
template<int dim>
int ImplicitLevelSetTsDesc<dim>::solveNonLinearSystem(DistSVec<double,dim> &U)
{
  
  int its;

  double t0 = this->timer->getTime();
  its = this->ns->solve(U);
  this->timer->addFluidSolutionTime(t0);

  if(!(this->interfaceType==MultiFluidData::FSF)){
    this->varFcn->conservativeToPrimitive(U,this->V0,&(this->fluidSelector.fluidId));
    //this->spaceOp->storePreviousPrimitive(U, this->Vg, this->fluidSelector.fluidId, this->Vgf, 
    //                                      this->Vgfweight, *this->X);
    this->riemann->storePreviousPrimitive(this->V0, this->fluidSelector.fluidId, *this->X);

    double t1 = this->timer->getTime();
    int itsLS = this->ns->solveLS(this->Phi, U);
    avoidNewPhaseCreation(this->Phi);
    (this->fluidSelector).getFluidId(this->Phi);
    this->timer->addLevelSetSolutionTime(t1);

    //this->spaceOp->updatePhaseChange(this->Vg, U, this->fluidSelector.fluidId, this->fluidSelector.fluidIdn,
    //                                 this->Vgf, this->Vgfweight, 
    //                                 this->riemann);
    this->riemann->updatePhaseChange(this->V0, this->fluidSelector.fluidId, this->fluidSelector.fluidIdn);
    this->varFcn->primitiveToConservative(this->V0,U,&(this->fluidSelector.fluidId));
  }

  checkSolution(U);

  return its;
  
}
//------------------------------------------------------------------------------
// External routines to solve Euler equations implicitly (called by NewtonSolver)
//------------------------------------------------------------------------------

// this function evaluates (Aw),t + F(w,x,v)
template<int dim>
void ImplicitLevelSetTsDesc<dim>::computeFunction(int it, DistSVec<double,dim> &Q,
                                                  DistSVec<double,dim> &F)
{
  // phi is obtained once and for all for this iteration
  // no need to recompute it before computation of jacobian.
  this->LS->conservativeToPrimitive(this->Phi,this->PhiV,Q);
  this->spaceOp->computeResidual(*this->X, *this->A, Q, this->PhiV, this->fluidSelector.fluidId, F, this->riemann, it+1);
  this->timeState->add_dAW_dt(it, *this->geoState, *this->A, Q, F);
  this->spaceOp->applyBCsToResidual(Q, F);
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitLevelSetTsDesc<dim>::recomputeFunction(DistSVec<double,dim> &Q,
                                            DistSVec<double,dim> &rhs)
{
  this->spaceOp->recomputeRHS(*this->X, Q, this->fluidSelector.fluidId, rhs);
}

//------------------------------------------------------------------------------

template<int dim>
int ImplicitLevelSetTsDesc<dim>::checkFailSafe(DistSVec<double,dim>& U)
{
  fprintf(stdout, "CheckFailSafe for ImplicitLevelSetTsDesc to be rewritten\n");

  if (!this->failSafeNewton) return 0;

  if (!this->tag)
    this->tag = new DistSVec<bool,2>(this->getVecInfo());

  this->domain->checkFailSafe(this->varFcn, U, *this->tag, &(this->fluidSelector.fluidId));
  this->spaceOp->fix(*this->tag);

  return 1;

}

//------------------------------------------------------------------------------
template<int dim> 
void ImplicitLevelSetTsDesc<dim>::resetFixesTag()
{

  this->spaceOp->resetTag();

}

//------------------------------------------------------------------------------
template<int dim>
void ImplicitLevelSetTsDesc<dim>::computeJacobian(int it, DistSVec<double,dim> &Q,
                                                         DistSVec<double,dim> &F)
{
  this->mvp->evaluate(it, *this->X, *this->A, Q, this->PhiV, this->fluidSelector.fluidId, this->riemann, F);
}
//------------------------------------------------------------------------------
template<int dim>
void ImplicitLevelSetTsDesc<dim>::setOperators(DistSVec<double,dim> &Q)
{
  
  DistMat<PrecScalar,dim> *_pc = dynamic_cast<DistMat<PrecScalar,dim> *>(pc);
  
  if (_pc) {
    
    MatVecProdFD<dim, dim> *mvpfd = dynamic_cast<MatVecProdFD<dim, dim> *>(mvp);
    MatVecProdH1<dim,MatScalar,dim> *mvph1 = dynamic_cast<MatVecProdH1<dim,MatScalar,dim> *>(mvp);
    MatVecProdH2<MatScalar,dim> *mvph2 = dynamic_cast<MatVecProdH2<MatScalar,dim> *>(mvp);
    
    if (mvpfd || mvph2) {

      this->spaceOp->computeJacobian(*this->X, *this->A, Q, *_pc, this->fluidSelector.fluidId, this->riemann);
      this->timeState->addToJacobian(*this->A, *_pc, Q);
      this->spaceOp->applyBCsToJacobian(Q, *_pc);
    }
    else if (mvph1) {
      JacobiPrec<PrecScalar,dim> *jac = dynamic_cast<JacobiPrec<PrecScalar,dim> *>(pc);
      IluPrec<PrecScalar,dim> *ilu = dynamic_cast<IluPrec<PrecScalar,dim> *>(pc);
      
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
int ImplicitLevelSetTsDesc<dim>::solveLinearSystem(int it, DistSVec<double,dim> &b,
						  DistSVec<double,dim> &dQ)
{
  
  double t0 = this->timer->getTime();
  dQ = 0.0;
  
  ksp->setup(it, this->maxItsNewton, b);
  
  int lits = ksp->solve(b, dQ);
  
  this->timer->addKspTime(t0);
  
  return lits;
  
}

//------------------------------------------------------------------------------
// External routines to solve LevelSet equation implicitly (called by NewtonSolver)
//------------------------------------------------------------------------------
// this function evaluates (Aw),t + F(w,x,v)
template<int dim>
void ImplicitLevelSetTsDesc<dim>::computeFunctionLS(int it,
                                                    DistSVec<double,dim> &U,
                                                    DistVec<double> &PhiF)
{
  this->spaceOp->computeResidualLS(*this->X, *this->A, this->Phi, this->fluidSelector.fluidId, U, PhiF);

  this->timeState->add_dAW_dtLS(it, *this->geoState, *this->A, this->Phi,
			        this->LS->Phin, this->LS->Phinm1, 
      				this->LS->Phinm2, PhiF);

}

//------------------------------------------------------------------------------
template<int dim>
void ImplicitLevelSetTsDesc<dim>::computeJacobianLS(int it,
                                                    DistSVec<double,dim> &U,
                                                    DistVec<double> &PhiF)
{

  mvpLS->evaluateLS(it, *this->X, *this->A, this->Phi, this->LS->Phin, 
		    this->LS->Phinm1, this->LS->Phinm2, U, PhiF, this->fluidSelector.fluidId);

}
//------------------------------------------------------------------------------

template<int dim>
int ImplicitLevelSetTsDesc<dim>::solveLinearSystemLS(int it, DistVec<double> &b,
                                                  DistVec<double> &dQ)
{

  double t0 = this->timer->getTime();

  dQ = 0.0;

  kspLS->setup(it, this->maxItsNewton, b);

  int lits = kspLS->solveLS(b, dQ);

  this->timer->addLSKspTime(t0);

  return lits;

}

//------------------------------------------------------------------------------

