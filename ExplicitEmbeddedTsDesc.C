
#include <GeoSource.h>
#include <DistTimeState.h>
#include <SpaceOperator.h>
#include <Domain.h>

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
ExplicitEmbeddedTsDesc<dim>::
ExplicitEmbeddedTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  EmbeddedTsDesc<dim>(ioData, geoSource, dom),
  k1(this->getVecInfo()), k2(this->getVecInfo()),
  k3(this->getVecInfo()), k4(this->getVecInfo()),
  p1(this->getVecInfo()), p2(this->getVecInfo()),
  p3(this->getVecInfo()), p4(this->getVecInfo()),
  U0(this->getVecInfo()), Phi0(this->getVecInfo())
{
  timeType = ioData.ts.expl.type;

  if(ioData.ts.expl.type == ExplicitData::RUNGE_KUTTA_4) RK4 = true;
  else {
    RK4 = false;
    if(ioData.ts.expl.type == ExplicitData::FORWARD_EULER) FE = true;
    else FE = false;
  } // if RK4 and FE are both false, do RK2.

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
}

//------------------------------------------------------------------------------

template<int dim>
ExplicitEmbeddedTsDesc<dim>::~ExplicitEmbeddedTsDesc()
{
}

//------------------------------------------------------------------------------

template<int dim>
int ExplicitEmbeddedTsDesc<dim>::solveNonLinearSystem(DistSVec<double,dim>& U,int)
{
  solveNLSystemOneBlock(U);
  return 1;
} //so far only consider one system to solve (there is no level-set)

//-----------------------------------------------------------------------------

template<int dim>
void ExplicitEmbeddedTsDesc<dim>::solveNLSystemOneBlock(DistSVec<double,dim> &U)
{
  double t0 = this->timer->getTime();
  DistSVec<double,dim> Ubc(this->getVecInfo());

  commonPart(U);
  if(RK4)     solveNLAllRK4(U,t0,Ubc);
  else if(FE) solveNLAllFE(U,t0,Ubc);
  else        solveNLAllRK2(U,t0,Ubc);

  this->updateBoundaryExternalState();
} 

//------------------------------------------------------------------------------

template<int dim>
void ExplicitEmbeddedTsDesc<dim>::commonPart(DistSVec<double,dim> &U)
{
  // Adam 04/06/10: Took everything in common in solveNLAllFE and solveNLAllRK2 and put it here. Added Ghost-Points treatment for viscous flows.

  if(this->mmh && !this->inSubCycling) {
    //get structure timestep dts
    this->dts = this->mmh->update(0, 0, 0, this->bcData->getVelocityVector(), *this->Xs);

    //recompute intersections
    double tw = this->timer->getTime();
    this->distLSS->recompute(this->dtf, this->dtfLeft, this->dts, true, TsDesc<dim>::failSafeFlag);
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
    this->spaceOp->updateSweptNodes(*this->X, this->phaseChangeChoice, this->phaseChangeAlg, U, this->Vtemp,
            *this->Weights, *this->VWeights, *this->Wstarij, *this->Wstarji,
            this->distLSS, (double*)this->vfar, (this->numFluid == 1 ? (DistVec<int>*)0 : &this->nodeTag));
    this->timer->addEmbedPhaseChangeTime(tw);
    this->timer->removeIntersAndPhaseChange(tw);

  }

  // Reset countWstar if second-order surrogate interface treatment is chosen
  if (this->interfaceAlg) {
    *this->countWstarij = 0;
    *this->countWstarji = 0;
    *this->Wstarij = 0.0;
    *this->Wstarji = 0.0;
  }

  // Ghost-Points Population
  if(this->eqsType == EmbeddedTsDesc<dim>::NAVIER_STOKES)
    {
      this->ghostPoints->deletePointers();
      this->spaceOp->populateGhostPoints(this->ghostPoints,*this->X,U,this->varFcn,this->distLSS,this->linRecAtInterface,this->nodeTag);
    }
}
//------------------------------------------------------------------------------

template<int dim>
void ExplicitEmbeddedTsDesc<dim>::solveNLAllFE(DistSVec<double,dim> &U, double t0, DistSVec<double,dim> &Ubc)
{
  computeRKUpdate(U, k1, 1);

  this->spaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - k1;
  this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  this->spaceOp->applyBCsToSolutionVector(U0,this->distLSS); //(?)for Navier-Stokes only

  U = U0;
  this->checkSolution(U);

  this->timer->addFluidSolutionTime(t0);
}

//------------------------------------------------------------------------------

template<int dim>
void ExplicitEmbeddedTsDesc<dim>::solveNLAllRK2(DistSVec<double,dim> &U, double t0, DistSVec<double,dim> &Ubc)
{
  computeRKUpdate(U, k1, 1);
  this->spaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - k1;
  this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  this->checkSolution(U0);

  // Ghost-Points Population
  if(this->eqsType == EmbeddedTsDesc<dim>::NAVIER_STOKES)
    {
      this->ghostPoints->deletePointers();
      this->spaceOp->populateGhostPoints(this->ghostPoints,*this->X,U0,this->varFcn,this->distLSS,this->linRecAtInterface,this->nodeTag);
    }
  computeRKUpdate(U0, k2, 1);
  this->spaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U = U - 1.0/2.0 * (k1 + k2);
  this->spaceOp->applyExtrapolationToSolutionVector(U, Ubc);

  this->spaceOp->applyBCsToSolutionVector(U,this->distLSS);

  this->checkSolution(U);

  this->timer->addFluidSolutionTime(t0);
}

//-----------------------------------------------------------------------------

template<int dim>
void ExplicitEmbeddedTsDesc<dim>::solveNLAllRK4(DistSVec<double,dim> &U, double t0, DistSVec<double,dim> &Ubc)
{ //TODO: no Ghost-Points Population ???
  computeRKUpdate(U, k1, 1);
  this->spaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - k1;
  this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  this->checkSolution(U0);

  computeRKUpdate(U0, k2, 1);
  this->spaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U0 = U - 0.5 * k2;
  this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  this->checkSolution(U0);

  computeRKUpdate(U0, k3, 1);
  this->spaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U0 = U - k3;
  this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  this->checkSolution(U0);

  computeRKUpdate(U0, k4, 1);
  this->spaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U = U - 1.0/6.0 * (k1 + 2.0 * (k2 + k3) + k4);
  this->spaceOp->applyExtrapolationToSolutionVector(U, Ubc);

  this->spaceOp->applyBCsToSolutionVector(U,this->distLSS);

  this->checkSolution(U);

  this->timer->addFluidSolutionTime(t0);
}

//-----------------------------------------------------------------------------

template<int dim>
void ExplicitEmbeddedTsDesc<dim>::computeRKUpdate(DistSVec<double,dim>& Ulocal,
                                  DistSVec<double,dim>& dU, int it)
//KW: 'it', positive or not, determines if Wstar, the Riemann solution (per edge), will be updated.
//    Only its sign (+ or 0) is used. Wstar is used for two purposes. 1) linear reconstruction at interface; 2) phase-change update
{
  this->spaceOp->applyBCsToSolutionVector(Ulocal,this->distLSS); //KW: (?)only for Navier-Stokes.
  if (this->interfaceAlg)
    this->spaceOp->computeResidual(*this->X, *this->A, Ulocal, *this->Wstarij, *this->Wstarji, *this->countWstarij, *this->countWstarji, this->distLSS,
                                   this->linRecAtInterface, this->nodeTag, dU, this->riemann, this->riemannNormal, this->Nsbar, this->timeState->getTime(), this->intersectAlpha, it, this->ghostPoints);
  else
    this->spaceOp->computeResidual(*this->X, *this->A, Ulocal, *this->Wstarij, *this->Wstarji, this->distLSS,
                                   this->linRecAtInterface, this->nodeTag, dU, this->riemann, this->riemannNormal, this->Nsbar, it, this->ghostPoints);

  this->timeState->multiplyByTimeStep(dU);
  
  if(this->numFluid==1&&!this->mmh)
    this->timeState->multiplyByPreconditioner(Ulocal,dU);
      //KW:This is the TEMPORAL low-mach precondition which is only for steady-state sims.
}


//------------------------------------------------------------------------------

