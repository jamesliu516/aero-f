#include <IoData.h>
#include <Domain.h>
#include <GeoSource.h>
#include <DistTimeState.h>
#include <SpaceOperator.h>
#include <LevelSet.h>
#include <FluidSelector.h>

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

template<int dim, int dimLS>
ExplicitMultiPhysicsTsDesc<dim,dimLS>::
ExplicitMultiPhysicsTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  MultiPhysicsTsDesc<dim,dimLS>(ioData, geoSource, dom),
  k1(this->getVecInfo()), k2(this->getVecInfo()),
  k3(this->getVecInfo()), k4(this->getVecInfo()), Ubc(this->getVecInfo()), 
  p1(this->getVecInfo()), p2(this->getVecInfo()),
  p3(this->getVecInfo()), p4(this->getVecInfo()), 
  U0(this->getVecInfo()), Phi0(this->getVecInfo()),
  fluidId0(this->getVecInfo())
{
  timeType = ioData.ts.expl.type;
  
  //initialize mmh (EmbeddedMeshMotionHandler).
  if(this->dynNodalTransfer) {
    MeshMotionHandler *_mmh = 0;
    _mmh = new EmbeddedMeshMotionHandler(ioData, dom, this->dynNodalTransfer, this->distLSS);
    this->mmh = _mmh;
  } else this->mmh = 0;
}

//------------------------------------------------------------------------------
  
template<int dim, int dimLS>
ExplicitMultiPhysicsTsDesc<dim,dimLS>::~ExplicitMultiPhysicsTsDesc()
{ /* nothing to destroy */ } 
  
//------------------------------------------------------------------------------

template<int dim, int dimLS>
int ExplicitMultiPhysicsTsDesc<dim,dimLS>::solveNonLinearSystem(DistSVec<double,dim> &U)
{
  solveNLSystemTwoBlocks(U);
  return 1;
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::solveNLSystemTwoBlocks(DistSVec<double,dim> &U)
{
  if(this->mmh && !this->inSubCycling) {
    // get structural time-step and recompute FS intersections.
    recomputeIntersections();
    // update fluidId.
    updateFluidIdFS(U);
    // update the phase-change (U & Phi) caused by the motion of FS interface
    updatePhaseChangeFS(U);
  }
  // populate ghost nodes (only for Navier-Stokes.)
  populateGhostPointsForNavierStokes(U);
  // evolve the fluid equation using FE, RK2, or RK4
  solveNLNavierStokes(U);
  // evolve the level-set equation using FE, RK2, or RK4.
  solveNLLevelSet(U);
  // update fluidId (fluidId0 = fluidId, fluidId = new).
  fluidId0 = *(this->fluidSelector.fluidId); // used in updatePhaseChangeFF
  this->fluidSelector.updateFluidIdFF(this->distLSS, this->Phi);
  // update the phase-change (only U) caused by the motion of FF interface
  updatePhaseChangeFF(U);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::recomputeIntersections()
{
  // just to get structural time-step (dts).
  this->dts = this->mmh->update(0, 0, 0, this->bcData->getVelocityVector(), *this->Xs);

  double tw = this->timer->getTime();
  this->distLSS->recompute(this->dtf, this->dtfLeft, this->dts);
  if(this->riemannNormal==2)
    this->multiPhaseSpaceOp->computeCellAveragedStructNormal(*(this->Nsbar), this->distLSS);
  this->timer->addIntersectionTime(tw);
  this->com->barrier();
  this->timer->removeIntersAndPhaseChange(tw);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::updatePhaseChangeFS(DistSVec<double,dim> &U)
{
  double tw = this->timer->getTime();
  switch(this->phaseChangeChoice) {
    case 0:
      this->multiPhaseSpaceOp->computeWeightsForEmbeddedStruct(*this->X, U, this->Vtemp, *this->Weights,
                                                     *this->VWeights, this->Phi, this->PhiWeights, this->distLSS, 
                                                     this->fluidSelector.fluidIdn, this->fluidSelector.fluidId);
      break;
    case 1:
      this->multiPhaseSpaceOp->computeRiemannWeightsForEmbeddedStruct(*this->X, U, this->Vtemp, *this->Wstarij,
                                                     *this->Wstarji, *this->Weights, *this->VWeights,
                                                     this->Phi, this->PhiWeights, this->distLSS, this->fluidSelector.fluidIdn,
                                                     this->fluidSelector.fluidId);
      break;
  }
  //update phase-change
  this->multiPhaseSpaceOp->updatePhaseChange(this->Vtemp, U, this->Weights, this->VWeights, &(this->Phi), &(this->PhiWeights), 
                                             this->distLSS, this->vfar, this->fluidSelector.fluidId);
                                            // Vtemp should have been filled in with primitive state
  this->timer->addEmbedPhaseChangeTime(tw);
  this->com->barrier();
  this->timer->removeIntersAndPhaseChange(tw);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::updateFluidIdFS(DistSVec<double,dim> &U)
{
  this->LS->conservativeToPrimitive(this->Phi, this->PhiV, U);
  this->multiPhaseSpaceOp->extrapolatePhiV(this->distLSS, this->PhiV);
  this->fluidSelector.updateFluidIdFS(this->distLSS, this->PhiV);
  this->PhiV = 0.0; //PhiV is no longer a distance function now. Only its sign (+/-)
                    //  is meaningful. We destroy it so people wouldn't use it
                    //  by mistake later on.
}

//------------------------------------------------------------------------------
//TODO: check with Adam to make sure this function works for Multi-Phase Fluid-Structure
template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::populateGhostPointsForNavierStokes(DistSVec<double,dim> &U)
{
  if(this->eqsType == MultiPhysicsTsDesc<dim,dimLS>::NAVIER_STOKES) {
    this->ghostPoints->deletePointers();
    this->multiPhaseSpaceOp->populateGhostPoints(this->ghostPoints,U,this->varFcn,this->distLSS,*(this->fluidSelector.fluidId));
  }
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::solveNLNavierStokes(DistSVec<double,dim> &U)
{
  double t0 = this->timer->getTime();
  switch (timeType) {
    case ExplicitData::FORWARD_EULER :
      solveNLNavierStokesFE(U); break;
    case ExplicitData::RUNGE_KUTTA_2 :
      solveNLNavierStokesRK2(U); break;
    case ExplicitData::RUNGE_KUTTA_4 :
      solveNLNavierStokesRK4(U); break;
    default:
      this->com->fprintf(stderr,"ERROR: Choose time-integrator from ForwardEuler, RungeKutta2, and RungeKutta4!\n");
  }

  // for FF phase-change update using extrapolation
  this->varFcn->conservativeToPrimitive(U,this->V0,this->fluidSelector.fluidId);
  this->riemann->storePreviousPrimitive(this->V0, *this->fluidSelector.fluidId, *this->X);

  this->timer->addFluidSolutionTime(t0);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::solveNLNavierStokesFE(DistSVec<double,dim> &U)
{
  this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);

  computeRKUpdate(U, k1, 1);
  this->multiPhaseSpaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - k1;
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  this->multiPhaseSpaceOp->applyBCsToSolutionVector(U0); //(?)for Navier-Stokes only

  U = U0;
  checkSolution(U);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::solveNLNavierStokesRK2(DistSVec<double,dim> &U)
{
  this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);

  computeRKUpdate(U, k1, 1);
  this->multiPhaseSpaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - k1;
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  checkSolution(U0);

  // Ghost-Points Population
  if(this->eqsType == MultiPhysicsTsDesc<dim,dimLS>::NAVIER_STOKES)
    {
      this->ghostPoints->deletePointers();
      this->multiPhaseSpaceOp->populateGhostPoints(this->ghostPoints,U,this->varFcn,this->distLSS, *this->fluidSelector.fluidId);
    }

  computeRKUpdate(U0, k2, 2);
  this->multiPhaseSpaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U = U - 1.0/2.0 * (k1 + k2);
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U, Ubc);
  this->multiPhaseSpaceOp->applyBCsToSolutionVector(U);
  checkSolution(U);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::solveNLNavierStokesRK4(DistSVec<double,dim> &U)
{ //TODO: no Ghost-Points Population ???
  this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);

  computeRKUpdate(U, k1, 1);
  this->multiPhaseSpaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - k1;
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  checkSolution(U0);

  computeRKUpdate(U0, k2, 2);
  this->multiPhaseSpaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U0 = U - 0.5 * k2;
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  checkSolution(U0);

  computeRKUpdate(U0, k3, 3);
  this->multiPhaseSpaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U0 = U - k3;
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  checkSolution(U0);

  computeRKUpdate(U0, k4, 4);
  this->multiPhaseSpaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U = U - 1.0/6.0 * (k1 + 2.0 * (k2 + k3) + k4);
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U, Ubc);
  this->multiPhaseSpaceOp->applyBCsToSolutionVector(U);
  checkSolution(U);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::computeRKUpdate(DistSVec<double,dim>& Ulocal,
                                                            DistSVec<double,dim>& dU, int it)
{
  this->multiPhaseSpaceOp->applyBCsToSolutionVector(Ulocal);
  this->multiPhaseSpaceOp->computeResidual(*this->X, *this->A, Ulocal, *this->Wstarij, *this->Wstarji,
                                           this->distLSS, this->linRecAtInterface, this->riemann, 
                                           this->riemannNormal, this->Nsbar, this->PhiV, this->fluidSelector,
                                           dU, it, this->ghostPoints);
                                           //Q: why send PhiV?
                                           //A: Riemann solver needs gradPhi.
                                           //Note: PhiV should be pre-computed. 
  this->timeState->multiplyByTimeStep(dU);

  if(this->numFluid==1&&!this->mmh)
    this->timeState->multiplyByPreconditioner(Ulocal,dU); //low-mach preconditioner
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::solveNLLevelSet(DistSVec<double,dim>& U)
{
  double t0 = this->timer->getTime();
  switch (timeType) {
    case ExplicitData::FORWARD_EULER :
      solveNLLevelSetFE(U); break;
    case ExplicitData::RUNGE_KUTTA_2 :
      solveNLLevelSetRK2(U); break;
    case ExplicitData::RUNGE_KUTTA_4 :
      solveNLLevelSetRK4(U); break;
    default:
      this->com->fprintf(stderr,"ERROR: Choose time-integrator from ForwardEuler, RungeKutta2, and RungeKutta4!\n");
  }
  this->timer->addLevelSetSolutionTime(t0);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::solveNLLevelSetFE(DistSVec<double,dim> &U)
{
  computeRKUpdateLS(this->Phi, *this->fluidSelector.fluidId, p1, U);
  Phi0 = this->Phi - p1;
  this->riemann->avoidNewPhaseCreation(this->Phi, this->LS->Phin, this->distLSS);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::solveNLLevelSetRK2(DistSVec<double,dim> &U)
{
  computeRKUpdateLS(this->Phi, *this->fluidSelector.fluidId, p1, U);
  Phi0 = this->Phi - p1;
//  this->riemann->avoidNewPhaseCreation(this->Phi0, this->LS->Phin, this->distLSS);
  this->fluidSelector.getFluidId(fluidId0,Phi0,&(this->distLSS->getStatus()));

  computeRKUpdateLS(Phi0, fluidId0, p2, U);
  this->Phi = this->Phi - 1.0/2.0 * (p1+p2);
  this->riemann->avoidNewPhaseCreation(this->Phi, this->LS->Phin, this->distLSS);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::solveNLLevelSetRK4(DistSVec<double,dim> &U)
{
  computeRKUpdateLS(this->Phi, *this->fluidSelector.fluidId, p1, U);
  Phi0 = this->Phi - 0.5 * p1;
//  this->riemann->avoidNewPhaseCreation(this->Phi0, this->LS->Phin, this->distLSS);
  this->fluidSelector.getFluidId(fluidId0,Phi0,&(this->distLSS->getStatus()));

  computeRKUpdateLS(Phi0, fluidId0, p2, U);
  Phi0 = this->Phi - 0.5 * p2;
//  this->riemann->avoidNewPhaseCreation(this->Phi0, this->LS->Phin, this->distLSS);
  this->fluidSelector.getFluidId(fluidId0,Phi0,&(this->distLSS->getStatus()));

  computeRKUpdateLS(Phi0, fluidId0, p3, U);
  Phi0 = this->Phi - p3;
//  this->riemann->avoidNewPhaseCreation(this->Phi0, this->LS->Phin, this->distLSS);
  this->fluidSelector.getFluidId(fluidId0,Phi0,&(this->distLSS->getStatus()));

  computeRKUpdateLS(Phi0, fluidId0, p4, U);
  this->Phi -= 1.0/6.0 * (p1 + 2.0 * (p2 + p3) + p4);
  this->riemann->avoidNewPhaseCreation(this->Phi, this->LS->Phin, this->distLSS);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::computeRKUpdateLS(DistSVec<double,dimLS> &Philocal,
                                            DistVec<int> &localFluidId,
                                            DistSVec<double,dimLS> &dPhi, DistSVec<double,dim> &U)
{
  this->multiPhaseSpaceOp->computeResidualLS(*this->X, *this->A, Philocal, localFluidId, U, dPhi, this->distLSS, this->linRecAtInterface);
  this->timeState->multiplyByTimeStep(dPhi);
  this->LS->checkTrueLevelSetUpdate(dPhi);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::updatePhaseChangeFF(DistSVec<double,dim> &U)
{
  this->riemann->updatePhaseChange(this->V0, *this->fluidSelector.fluidId, fluidId0);
  this->varFcn->primitiveToConservative(this->V0,U,this->fluidSelector.fluidId);
  checkSolution(U);
}

//------------------------------------------------------------------------------








