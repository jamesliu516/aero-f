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
ImplicitEmbeddedTsDesc<dim>::
ImplicitEmbeddedTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  EmbeddedTsDesc<dim>(ioData, geoSource, dom)
{
  tag = 0;
  ImplicitData &implicitData = ioData.ts.implicit;
  
  // NewtonSolver
  ns = new NewtonSolver<ImplicitEmbeddedTsDesc<dim> >(this);
  failSafeNewton = implicitData.newton.failsafe;
  maxItsNewton = implicitData.newton.maxIts;
  epsNewton = implicitData.newton.eps;
  epsAbsResNewton = implicitData.newton.epsAbsRes;
  epsAbsIncNewton = implicitData.newton.epsAbsInc;

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
  
  pc = createPreconditioner<dim>(implicitData.newton.ksp.ns.pc, this->domain);
  
  ksp = createKrylovSolver(this->getVecInfo(), implicitData.newton.ksp.ns, mvp, pc, this->com);
  
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

  typename MatVecProd<dim,dim>::_fsi fsi = {

    this->distLSS,
    &this->nodeTag,
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

template<int dim>
ImplicitEmbeddedTsDesc<dim>::~ImplicitEmbeddedTsDesc()
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

template<int dim>
template <int neq>
KspPrec<neq> *ImplicitEmbeddedTsDesc<dim>::createPreconditioner(PcData &pcdata, Domain *dom)
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

template<int dim>
template<int neq, class MatVecProdOp>
KspSolver<DistSVec<double,neq>, MatVecProdOp, KspPrec<neq>, Communicator> *
ImplicitEmbeddedTsDesc<dim>::createKrylovSolver(
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

template<int dim>
void ImplicitEmbeddedTsDesc<dim>::commonPart(DistSVec<double,dim> &U)
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
      this->spaceOp->computeCellAveragedStructNormal(*(this->Nsbar), this->distLSS);

    //update nodeTags (only for numFluid>1)
 //   if(this->numFluid>1) {
      this->nodeTag0 = this->nodeTag;
      this->nodeTag = this->distLSS->getStatus();
 //   }

    //store previous states for phase-change update
    tw = this->timer->getTime();
    switch(this->phaseChangeChoice) {
      case 0:
        if(this->numFluid==1)
          this->spaceOp->computeWeightsForEmbeddedStruct(*this->X, U, this->Vtemp, *this->Weights,
                                                         *this->VWeights, this->distLSS);
        else //numFluid>1
          this->spaceOp->computeWeightsForEmbeddedStruct(*this->X, U, this->Vtemp, *this->Weights,
                                                         *this->VWeights, this->distLSS, &this->nodeTag0);
        break;
      case 1:
        if(this->numFluid==1)
          this->spaceOp->computeRiemannWeightsForEmbeddedStruct(*this->X, U, this->Vtemp, *this->Wstarij,
                                                         *this->Wstarji, *this->Weights, *this->VWeights,
                                                         this->distLSS);
        else //numFluid>1
          this->spaceOp->computeRiemannWeightsForEmbeddedStruct(*this->X, U, this->Vtemp, *this->Wstarij,
                                                         *this->Wstarji, *this->Weights, *this->VWeights,
                                                         this->distLSS, &this->nodeTag0);
        break;
    }
    this->timer->addEmbedPhaseChangeTime(tw);
    this->com->barrier();
    this->timer->removeIntersAndPhaseChange(tw);

    //update phase-change
    tw = this->timer->getTime();
    if(this->numFluid==1)
      this->spaceOp->updatePhaseChange(this->Vtemp, U, this->Weights, this->VWeights, this->distLSS, this->vfar);
    else //numFluid>1
      this->spaceOp->updatePhaseChange(this->Vtemp, U, this->Weights, this->VWeights, this->distLSS, this->vfar, &this->nodeTag);
   
    //this->timeState->update(U); 
    this->timeState->getUn() = U;

    // BDF update (Unm1)
    if (this->timeState->useNm1() && this->timeState->existsNm1()) {
      tw = this->timer->getTime();
      DistSVec<double,dim>& Unm1 = this->timeState->getUnm1();
  
      if (!this->existsWstarnm1) {
        
        this->spaceOp->computeResidual(*this->X, *this->A, Unm1, *this->Wstarij, *this->Wstarji, this->distLSS,
                                 this->linRecAtInterface, this->nodeTag, this->Vtemp, this->riemann, 
                                 this->riemannNormal, this->Nsbar, 1, this->ghostPoints);
      }

      switch(this->phaseChangeChoice) {
        case 0:
          if(this->numFluid==1)
            this->spaceOp->computeWeightsForEmbeddedStruct(*this->X, Unm1, this->Vtemp, *this->Weights,
                                                           *this->VWeights, this->distLSS);
          else //numFluid>1
            this->spaceOp->computeWeightsForEmbeddedStruct(*this->X, Unm1, this->Vtemp, *this->Weights,
                                                         *this->VWeights, this->distLSS, &this->nodeTag0);
          break;
      case 1:
        if(this->numFluid==1)
          this->spaceOp->computeRiemannWeightsForEmbeddedStruct(*this->X, Unm1, this->Vtemp, *this->Wstarij_nm1,
                                                         *this->Wstarji_nm1, *this->Weights, *this->VWeights,
                                                         this->distLSS);
        else //numFluid>1
          this->spaceOp->computeRiemannWeightsForEmbeddedStruct(*this->X, Unm1, this->Vtemp, *this->Wstarij_nm1,
                                                         *this->Wstarji_nm1, *this->Weights, *this->VWeights,
                                                         this->distLSS, &this->nodeTag0);
        break;
      }
      this->timer->addEmbedPhaseChangeTime(tw);
      this->com->barrier();
      this->timer->removeIntersAndPhaseChange(tw);

      //update phase-change
      tw = this->timer->getTime();
      if(this->numFluid==1)
        this->spaceOp->updatePhaseChange(this->Vtemp, Unm1, this->Weights, this->VWeights, this->distLSS, this->vfar);
      else //numFluid>1
        this->spaceOp->updatePhaseChange(this->Vtemp, Unm1, this->Weights, this->VWeights, this->distLSS, this->vfar, &this->nodeTag);
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
  if(this->eqsType == EmbeddedTsDesc<dim>::NAVIER_STOKES)
    {
      this->ghostPoints->deletePointers();
      this->spaceOp->populateGhostPoints(this->ghostPoints,U,this->varFcn,this->distLSS,this->nodeTag);
    }
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// External routine to solve problem (called by TsSolver)
// It calls for the NewtonSolver ns, which in turn will
// call routines below from this same file or from LevelSetTsDesc
//------------------------------------------------------------------------------
template<int dim>
int ImplicitEmbeddedTsDesc<dim>::solveNonLinearSystem(DistSVec<double,dim> &U)
{ 
  double t0 = this->timer->getTime();
  DistSVec<double,dim> Ubc(this->getVecInfo());

  commonPart(U);

  int its = this->ns->solve(U);
  
  this->timer->addFluidSolutionTime(t0);
   
  checkSolution(U);

  return its;
}
//------------------------------------------------------------------------------
// External routines to solve Euler equations implicitly (called by NewtonSolver)
//------------------------------------------------------------------------------

// this function evaluates (Aw),t + F(w,x,v)
template<int dim>
void ImplicitEmbeddedTsDesc<dim>::computeFunction(int it, DistSVec<double,dim> &Q,
                                                  DistSVec<double,dim> &F)
{
  // phi is obtained once and for all for this iteration
  // no need to recompute it before computation of jacobian.
  this->spaceOp->computeResidual(*this->X, *this->A, Q, *this->Wstarij, *this->Wstarji, this->distLSS,
                                 this->linRecAtInterface, this->nodeTag, F, this->riemann, 
                                 this->riemannNormal, this->Nsbar, 1, this->ghostPoints);
  this->timeState->add_dAW_dt(it, *this->geoState, *this->A, Q, F);
  this->spaceOp->applyBCsToResidual(Q, F);
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitEmbeddedTsDesc<dim>::recomputeFunction(DistSVec<double,dim> &Q,
                                            DistSVec<double,dim> &rhs)
{
  this->spaceOp->recomputeRHS(*this->X, Q, rhs);
}

//------------------------------------------------------------------------------

template<int dim>
int ImplicitEmbeddedTsDesc<dim>::checkFailSafe(DistSVec<double,dim>& U)
{
  this->com->fprintf(stdout, "WARNING: At the moment CheckFailSafe is not supported by the embedded framework with an implicit time-integrator!\n");
/*
  if (!this->failSafeNewton) return 0;

  if (!this->tag)
    this->tag = new DistSVec<bool,2>(this->getVecInfo());

  this->domain->checkFailSafe(this->varFcn, U, *this->tag, this->fluidSelector.fluidId);
  this->multiPhaseSpaceOp->fix(*this->tag);

  return 1;
*/
}

//------------------------------------------------------------------------------
template<int dim> 
void ImplicitEmbeddedTsDesc<dim>::resetFixesTag()
{

  this->spaceOp->resetTag();

}

//------------------------------------------------------------------------------
template<int dim>
void ImplicitEmbeddedTsDesc<dim>::computeJacobian(int it, DistSVec<double,dim> &Q,
							DistSVec<double,dim> &F)
{

  mvp->evaluate(it,*(this->X) ,*(this->A), Q, F);

}
//------------------------------------------------------------------------------
template<int dim>
void ImplicitEmbeddedTsDesc<dim>::setOperators(DistSVec<double,dim> &Q)
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
int ImplicitEmbeddedTsDesc<dim>::solveLinearSystem(int it, DistSVec<double,dim> &b,
				                   DistSVec<double,dim> &dQ)
{
  
  double t0 = this->timer->getTime();
  dQ = 0.0;
  
  ksp->setup(it, this->maxItsNewton, b);
  
  int lits = ksp->solve(b, dQ);
  
  this->timer->addKspTime(t0);
  
  return lits;
  
}
