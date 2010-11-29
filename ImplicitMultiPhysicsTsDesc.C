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

template<int dim,int dimLS>
ImplicitMultiPhysicsTsDesc<dim,dimLS>::
ImplicitMultiPhysicsTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  MultiPhysicsTsDesc<dim,dimLS>(ioData, geoSource, dom)
{
  tag = 0;
  ImplicitData &implicitData = ioData.ts.implicit;
  
  // NewtonSolver
  ns = new NewtonSolver<ImplicitMultiPhysicsTsDesc<dim,dimLS> >(this);
  failSafeNewton = implicitData.newton.failsafe;
  maxItsNewton = implicitData.newton.maxIts;
  epsNewton = implicitData.newton.eps;
  epsAbsResNewton = implicitData.newton.epsAbsRes;
  epsAbsIncNewton = implicitData.newton.epsAbsInc;

  // MatVecProd, Prec and Krylov solver for Euler equations
  if (implicitData.mvp == ImplicitData::FD)
    mvp = new MatVecProdFDMultiPhase<dim,dimLS>(this->timeState, this->geoState,
                                      this->multiPhaseSpaceOp,this->riemann,&(this->fluidSelector),
                                      this->domain,ioData);
  else if (implicitData.mvp == ImplicitData::H1)
    mvp = new MatVecProdH1MultiPhase<dim,dimLS>(this->timeState, this->multiPhaseSpaceOp, this->riemann,&(this->fluidSelector),this->domain);
  else{
    this->com->fprintf(stdout, "*** Error: MatVecProdH2 is not available\n");
    exit(1);
  }
  
  pc = createPreconditioner<dim>(implicitData.newton.ksp.ns.pc, this->domain);
  
  ksp = createKrylovSolver(this->getVecInfo(), implicitData.newton.ksp.ns, mvp, pc, this->com);

  // MatVecProd, Prec and Krylov solver for LevelSet equation
  mvpLS  = new MatVecProdLS<dim,dimLS>(this->timeState, this->geoState, this->multiPhaseSpaceOp, this->domain, this->LS);

  pcLS = createPreconditioner<dimLS>(implicitData.newton.ksp.lsi.pc, this->domain);
  kspLS = createKrylovSolver(this->getVecInfo(), implicitData.newton.ksp.lsi, mvpLS, pcLS, this->com);
 
 
  //initialize mmh (EmbeddedMeshMotionHandler).
  if(this->dynNodalTransfer) 
    {
      /*
      MeshMotionHandler *_mmh = 0;
      _mmh = new EmbeddedMeshMotionHandler(ioData, dom, this->dynNodalTransfer, this->distLSS);
      this->mmh = _mmh;
      */
      this->mmh = new EmbeddedMeshMotionHandler(ioData, dom, this->dynNodalTransfer, this->distLSS);
    } 
  else
    { 
      this->mmh = 0;
    }

  typename MatVecProdMultiPhase<dim,dimLS>::_fsi fsi = {

    this->distLSS,
    this->fluidSelector.fluidId,
    this->riemann,
    this->linRecAtInterface,
    this->Nsbar,
    &this->Wtemp,
    this->riemannNormal,
    this->ghostPoints,
  };

  mvp->AttachStructure(fsi);
  
  this->existsWstarnm1 = false;
}

//------------------------------------------------------------------------------

template<int dim,int dimLS>
ImplicitMultiPhysicsTsDesc<dim,dimLS>::~ImplicitMultiPhysicsTsDesc()
{
  if (tag)   delete tag;
  if (mvp)   delete mvp;
  if (pc)    delete pc;
  if (ksp)   delete ksp;
  if (ns)    delete ns;

}

//------------------------------------------------------------------------------
//  Internal routines to setup the class (called in constructor)
//------------------------------------------------------------------------------

template<int dim,int dimLS>
template <int neq>
KspPrec<neq> *ImplicitMultiPhysicsTsDesc<dim,dimLS>::createPreconditioner(PcData &pcdata, Domain *dom)
{
  
  KspPrec<neq> *_pc = 0;
  
  if (pcdata.type == PcData::IDENTITY)
    _pc = new IdentityPrec<neq>();
  else if (pcdata.type == PcData::JACOBI)
    _pc = new JacobiPrec<double,neq>(DiagMat<double,neq>::DENSE, dom);
  else if (pcdata.type == PcData::AS ||
	   pcdata.type == PcData::RAS ||
	   pcdata.type == PcData::ASH ||
	   pcdata.type == PcData::AAS)
    _pc = new IluPrec<double,neq>(pcdata, dom);
  
  return _pc;
  
}

//------------------------------------------------------------------------------

template<int dim,int dimLS>
template<int neq, class MatVecProdOp>
KspSolver<DistSVec<double,neq>, MatVecProdOp, KspPrec<neq>, Communicator> *
ImplicitMultiPhysicsTsDesc<dim,dimLS>::createKrylovSolver(
                               const DistInfo &info, KspData &kspdata,
                               MatVecProdOp *_mvp, KspPrec<neq> *_pc,
                               Communicator *_com)
{
  
  KspSolver<DistSVec<double,neq>, MatVecProdOp, KspPrec<neq>, Communicator> *_ksp = 0;
  
  if (kspdata.type == KspData::RICHARDSON)
    _ksp = new RichardsonSolver<DistSVec<double,neq>, MatVecProdOp,
                 KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  else if (kspdata.type == KspData::CG)
    _ksp = new CgSolver<DistSVec<double,neq>, MatVecProdOp,
                 KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  else if (kspdata.type == KspData::GMRES)
    _ksp = new GmresSolver<DistSVec<double,neq>, MatVecProdOp,
                 KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  
  return _ksp;
  
}

template<int dim,int dimLS>
void ImplicitMultiPhysicsTsDesc<dim,dimLS>::commonPart(DistSVec<double,dim> &U)
{
  // Adam 04/06/10: Took everything in common in solveNLAllFE and solveNLAllRK2 and put it here. Added Ghost-Points treatment for viscous flows.

  if(this->mmh && !this->inSubCycling) {
    //get structure timestep dts
    this->dts = this->mmh->update(0, 0, 0, this->bcData->getVelocityVector(), *this->Xs);

    //recompute intersections
    double tw = this->timer->getTime();
    this->distLSS->recompute(this->dtf, this->dtfLeft, this->dts); 
    this->timer->addIntersectionTime(tw);
    this->com->barrier();
    this->timer->removeIntersAndPhaseChange(tw);
    if(this->riemannNormal==2)
      this->multiPhaseSpaceOp->computeCellAveragedStructNormal(*(this->Nsbar), this->distLSS);
   
    this->LS->conservativeToPrimitive(this->Phi, this->PhiV, U);
    this->multiPhaseSpaceOp->extrapolatePhiV(this->distLSS, this->PhiV);
    this->fluidSelector.updateFluidIdFS(this->distLSS, this->PhiV);
     this->PhiV = 0.0; //PhiV is no longer a distance function now. Only its sign (+/-)
                       //  is meaningful. We destroy it so people wouldn't use it
                       //  by mistake later on.

    //store previous states for phase-change update
    tw = this->timer->getTime();
    switch(this->phaseChangeChoice) {
      case 0:
        this->multiPhaseSpaceOp->computeWeightsForEmbeddedStruct(*this->X, U, this->Vtemp, *this->Weights,
                                                         *this->VWeights, this->Phi,this->PhiWeights,this->distLSS,
                                                         this->fluidSelector.fluidIdn, this->fluidSelector.fluidId);
        break;
      case 1:
        this->multiPhaseSpaceOp->computeRiemannWeightsForEmbeddedStruct(*this->X, U, this->Vtemp, *this->Wstarij,
                                                         *this->Wstarji, *this->Weights, *this->VWeights,
                                                         this->Phi, this->PhiWeights, this->distLSS,
                                                         this->fluidSelector.fluidIdn, this->fluidSelector.fluidId);
        break;
    }
    this->timer->addEmbedPhaseChangeTime(tw);
    this->com->barrier();
    this->timer->removeIntersAndPhaseChange(tw);

    //update phase-change
    tw = this->timer->getTime();
    this->multiPhaseSpaceOp->updatePhaseChange(this->Vtemp, U, this->Weights, this->VWeights, &(this->Phi), &(this->PhiWeights),
                                               this->distLSS, this->vfar, this->fluidSelector.fluidId);
   
    this->timeState->getUn() = U;

    // BDF update (Unm1)
    if (this->timeState->useNm1() && this->timeState->existsNm1()) {
      tw = this->timer->getTime();
      DistSVec<double,dim>& Unm1 = this->timeState->getUnm1();
  
      if (!this->existsWstarnm1) {
        fprintf(stderr,"*** Error: I ignored this case!\n");
        exit(1);  
      }

      switch(this->phaseChangeChoice) {
      case 0:
        this->multiPhaseSpaceOp->computeWeightsForEmbeddedStruct(*this->X, U, this->Vtemp, *this->Weights,
                                                         *this->VWeights, this->Phi,this->PhiWeights,this->distLSS,
                                                         this->fluidSelector.fluidIdn, this->fluidSelector.fluidId);
        break;
      case 1:
        this->multiPhaseSpaceOp->computeRiemannWeightsForEmbeddedStruct(*this->X, U, this->Vtemp, *this->Wstarij_nm1,
                                                         *this->Wstarji_nm1, *this->Weights, *this->VWeights,
                                                         this->Phi, this->PhiWeights, this->distLSS,
                                                         this->fluidSelector.fluidIdn, this->fluidSelector.fluidId);
        break;
    }
    this->timer->addEmbedPhaseChangeTime(tw);

      this->timer->addEmbedPhaseChangeTime(tw);
      this->com->barrier();
      this->timer->removeIntersAndPhaseChange(tw);

      //update phase-change
      tw = this->timer->getTime();
      this->multiPhaseSpaceOp->updatePhaseChange(this->Vtemp, Unm1, this->Weights, this->VWeights, &(this->Phi), &(this->PhiWeights),
                                                 this->distLSS, this->vfar, this->fluidSelector.fluidId);
      this->timer->addEmbedPhaseChangeTime(tw);
      this->com->barrier();
      this->timer->removeIntersAndPhaseChange(tw);
    }

    if (this->timeState->useNm1()) {
      *this->Wstarij_nm1 = *this->Wstarij;
      *this->Wstarji_nm1 = *this->Wstarji;
      this->existsWstarnm1 = true;
    }
  
  }

  // Ghost-Points Population
/*  if(this->eqsType == TsDesc<dim>::NAVIER_STOKES)
    {
      this->ghostPoints->deletePointers();
      this->mulitPhaseSpaceOp->populateGhostPoints(this->ghostPoints,U,this->varFcn,this->distLSS,this->nodeTag);
    }
*/
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// External routine to solve problem (called by TsSolver)
// It calls for the NewtonSolver ns, which in turn will
// call routines below from this same file or from LevelSetTsDesc
//------------------------------------------------------------------------------
template<int dim,int dimLS>
int ImplicitMultiPhysicsTsDesc<dim,dimLS>::solveNonLinearSystem(DistSVec<double,dim> &U)
{ 
  double t0 = this->timer->getTime();
  DistSVec<double,dim> Ubc(this->getVecInfo());

  commonPart(U);

  int its = this->ns->solve(U);
  
  this->timer->addFluidSolutionTime(t0);
   
  this->varFcn->conservativeToPrimitive(U,this->V0,this->fluidSelector.fluidId);
  this->riemann->storePreviousPrimitive(this->V0, *this->fluidSelector.fluidId, *this->X);
    
  double t1 = this->timer->getTime();
  int itsLS = this->ns->solveLS(this->Phi, U);
  this->riemann->storeOldV(U);
  this->riemann->avoidNewPhaseCreation(this->Phi, this->LS->Phin,this->distLSS);
//  this->fluidSelector.getFluidId(this->Phi,&(this->distLSS->getStatus()));
  DistVec<int> fluidId0(*this->fluidSelector.fluidId);
  this->fluidSelector.updateFluidIdFF(this->distLSS, this->Phi);
  this->timer->addLevelSetSolutionTime(t1);

  this->riemann->updatePhaseChange(this->V0, *this->fluidSelector.fluidId, fluidId0);
  this->varFcn->primitiveToConservative(this->V0,U,this->fluidSelector.fluidId);

  checkSolution(U);

  return its;
}
//------------------------------------------------------------------------------
// External routines to solve Euler equations implicitly (called by NewtonSolver)
//------------------------------------------------------------------------------

// this function evaluates (Aw),t + F(w,x,v)
template<int dim,int dimLS>
void ImplicitMultiPhysicsTsDesc<dim,dimLS>::computeFunction(int it, DistSVec<double,dim> &Q,
                                                  DistSVec<double,dim> &F)
{
  // phi is obtained once and for all for this iteration
  // no need to recompute it before computation of jacobian.
  
  this->LS->conservativeToPrimitive(this->Phi,this->PhiV,Q);
  this->multiPhaseSpaceOp->computeResidual(*this->X, *this->A, Q, *this->Wstarij, *this->Wstarji, this->distLSS,
                                 this->linRecAtInterface, this->riemann,  
                                 this->riemannNormal, this->Nsbar, this->PhiV, this->fluidSelector,F, 1, this->ghostPoints);
  this->timeState->add_dAW_dt(it, *this->geoState, *this->A, Q, F);
  this->multiPhaseSpaceOp->applyBCsToResidual(Q, F);
}

//------------------------------------------------------------------------------

template<int dim,int dimLS>
void ImplicitMultiPhysicsTsDesc<dim,dimLS>::recomputeFunction(DistSVec<double,dim> &Q,
                                            DistSVec<double,dim> &rhs)
{
  this->multiPhaseSpaceOp->recomputeRHS(*this->X, Q, rhs);
}

//------------------------------------------------------------------------------

template<int dim,int dimLS>
int ImplicitMultiPhysicsTsDesc<dim,dimLS>::checkFailSafe(DistSVec<double,dim>& U)
{
//  this->com->fprintf(stdout, "WARNING: At the moment CheckFailSafe is not supported by the embedded framework with an implicit time-integrator!\n");

  if (!this->failSafeNewton) return 0;

  if (!this->tag)
    this->tag = new DistSVec<bool,2>(this->getVecInfo());

  this->domain->checkFailSafe(this->varFcn, U, *this->tag, this->fluidSelector.fluidId);
  this->multiPhaseSpaceOp->fix(*this->tag);

  return 1;

}

//------------------------------------------------------------------------------
template<int dim,int dimLS> 
void ImplicitMultiPhysicsTsDesc<dim,dimLS>::resetFixesTag()
{

  this->multiPhaseSpaceOp->resetTag();

}

//------------------------------------------------------------------------------
template<int dim,int dimLS>
void ImplicitMultiPhysicsTsDesc<dim,dimLS>::computeJacobian(int it, DistSVec<double,dim> &Q,
							DistSVec<double,dim> &F)
{

  mvp->evaluate(it,*(this->X) ,*(this->A), Q,this->PhiV, F);

}
//------------------------------------------------------------------------------
template<int dim,int dimLS>
void ImplicitMultiPhysicsTsDesc<dim,dimLS>::setOperators(DistSVec<double,dim> &Q)
{
  
  DistMat<double,dim> *_pc = dynamic_cast<DistMat<double,dim> *>(pc);
  
  if (_pc) {
    
    MatVecProdFDMultiPhase<dim,dimLS> *mvpfd = dynamic_cast<MatVecProdFDMultiPhase<dim,dimLS> *>(mvp);
    MatVecProdH1MultiPhase<dim,dimLS> *mvph1 = dynamic_cast<MatVecProdH1MultiPhase<dim,dimLS> *>(mvp);
    
    if (mvpfd) {

      this->multiPhaseSpaceOp->computeJacobian(this->riemann, *this->X, Q,*this->A,this->distLSS,
                                   this->riemannNormal, this->Nsbar,(this->fluidSelector),*_pc,this->timeState);
      this->timeState->addToJacobian(*this->A, *_pc, Q);
      this->multiPhaseSpaceOp->applyBCsToJacobian(Q, *_pc);
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
template<int dim,int dimLS>
int ImplicitMultiPhysicsTsDesc<dim,dimLS>::solveLinearSystem(int it, DistSVec<double,dim> &b,
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
template<int dim, int dimLS>
void ImplicitMultiPhysicsTsDesc<dim,dimLS>::computeFunctionLS(int it,
                                                    DistSVec<double,dim> &U,
                                                    DistSVec<double,dimLS> &PhiF)
{
  this->multiPhaseSpaceOp->computeResidualLS(*this->X, *this->A, this->Phi, *this->fluidSelector.fluidId, U, PhiF,this->distLSS,this->linRecAtInterface);

  this->timeState->add_dAW_dtLS(it, *this->geoState, *this->A, this->Phi,
                                this->LS->Phin, this->LS->Phinm1,
                                this->LS->Phinm2, PhiF,this->requireSpecialBDF);

}

//------------------------------------------------------------------------------
template<int dim, int dimLS>
void ImplicitMultiPhysicsTsDesc<dim,dimLS>::computeJacobianLS(int it,
                                                    DistSVec<double,dim> &U,
                                                    DistSVec<double,dimLS> &PhiF)
{
  mvpLS->evaluate(it, *this->X, *this->A, this->Phi,
                  U,this->V0, PhiF, *this->fluidSelector.fluidId,this->requireSpecialBDF,this->distLSS);
}

//------------------------------------------------------------------------------
template<int dim, int dimLS>
void ImplicitMultiPhysicsTsDesc<dim,dimLS>::setOperatorsLS(DistSVec<double,dimLS> &Q)
{

  DistMat<double,dimLS> *_pc = dynamic_cast<DistMat<double,dimLS> *>(pcLS);

  if (_pc) {

      JacobiPrec<double,dimLS> *jac = dynamic_cast<JacobiPrec<double,dimLS> *>(pcLS);
      IluPrec<double,dimLS> *ilu = dynamic_cast<IluPrec<double,dimLS> *>(pcLS);

      if (jac)
        jac->getData(*mvpLS);
      else if (ilu)
        ilu->getData(*mvpLS);

  }

  double t0 = this->timer->getTime();

  pcLS->setup();

  double t = this->timer->addLSPrecSetupTime(t0);

  this->com->printf(6, "Fluid preconditioner computation: %f s\n", t);

}
//------------------------------------------------------------------------------

template<int dim, int dimLS>
int ImplicitMultiPhysicsTsDesc<dim,dimLS>::solveLinearSystemLS(int it, DistSVec<double,dimLS> &b,
                                                  DistSVec<double,dimLS> &dQ)
{

  double t0 = this->timer->getTime();

  dQ = 0.0;
  
  kspLS->setup(it, this->maxItsNewton, b);
  
  int lits = kspLS->solve(b, dQ);

  //mvpLS->apply(dQ, fnew);

  this->timer->addLSKspTime(t0);

  return lits;

}
  
//------------------------------------------------------------------------------

