#include <ExplicitLevelSetTsDesc.h>

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
ExplicitLevelSetTsDesc<dim>::
ExplicitLevelSetTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  LevelSetTsDesc<dim>(ioData, geoSource, dom), 
  k1(this->getVecInfo()), k2(this->getVecInfo()), 
  k3(this->getVecInfo()), k4(this->getVecInfo()), 
  p1(this->getVecInfo()), p2(this->getVecInfo()), 
  p3(this->getVecInfo()), p4(this->getVecInfo()), 
  U0(this->getVecInfo()), Phi0(this->getVecInfo())
{
  this->mmh = this->createMeshMotionHandler(ioData, geoSource, 0);

  timeType = ioData.ts.expl.type;
}

//------------------------------------------------------------------------------

template<int dim>
ExplicitLevelSetTsDesc<dim>::~ExplicitLevelSetTsDesc()
{
}

//------------------------------------------------------------------------------

template<int dim>
int ExplicitLevelSetTsDesc<dim>::solveNonLinearSystem(DistSVec<double,dim> &U)
{

  if(timeType == ExplicitData::FORWARD_EULER || 
     timeType == ExplicitData::ONE_BLOCK_RK2 )
    solveNLSystemOneBlock(U);
  else if(timeType == ExplicitData::RUNGE_KUTTA_2 ||
          timeType == ExplicitData::RUNGE_KUTTA_4 )
    solveNLSystemTwoBlocks(U);

  return 1;

}

//------------------------------------------------------------------------------

template<int dim>
void ExplicitLevelSetTsDesc<dim>::solveNLSystemOneBlock(DistSVec<double,dim> &U)
{

  if(timeType == ExplicitData::FORWARD_EULER)
    solveNLAllFE(U);
  else if(timeType == ExplicitData::ONE_BLOCK_RK2)
    solveNLAllRK2(U);

}
//------------------------------------------------------------------------------

template<int dim>
void ExplicitLevelSetTsDesc<dim>::solveNLAllFE(DistSVec<double,dim> &U)
{

  double t0 = this->timer->getTime();

  DistSVec<double,dim> Ubc(this->getVecInfo());
  this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);
  computeRKUpdate(U, k1,1);
  this->spaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - k1;
  this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);

  this->timer->addFluidSolutionTime(t0);

  if(!(this->interfaceType==MultiFluidData::FSF)){
    this->spaceOp->storePreviousPrimitive(U0, this->Vg, this->Phi,
                                          this->Vgf, this->Vgfweight);
    t0 = this->timer->getTime();

    computeRKUpdateLS(this->Phi, p1, U);
    this->Phi = this->Phi - p1;

    this->timer->addLevelSetSolutionTime(t0);

    U = U0;

    // Riemann overwrite using the value of Phi_{n+1}
    this->spaceOp->updatePhaseChange(this->Vg, U, this->Phi,
                                     this->LS->Phin, this->Vgf,
                                     this->Vgfweight, this->riemann);
  }else U = U0;

  checkSolution(U);
  this->boundaryFlux = this->tmpDistSVec;

}

//------------------------------------------------------------------------------

template<int dim>
void ExplicitLevelSetTsDesc<dim>::solveNLAllRK2(DistSVec<double,dim> &U)
{

  double t0 = this->timer->getTime();

  DistSVec<double,dim> Ubc(this->getVecInfo());
  this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);
  // *** prediction step ***
  computeRKUpdate(U, k1,1);
  this->spaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - k1;
  this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);

  this->timer->addFluidSolutionTime(t0);

  if(!(this->interfaceType==MultiFluidData::FSF)){
    this->spaceOp->storePreviousPrimitive(U0, this->Vg, this->Phi,
                                          this->Vgf, this->Vgfweight);

    t0 = this->timer->getTime();

    computeRKUpdateLS(this->Phi, p1, U);
    Phi0 = this->Phi - p1;

    this->timer->addLevelSetSolutionTime(t0);

    // Riemann overwrite on U0 using the value of Phi0
    this->spaceOp->updatePhaseChange(this->Vg, U0, Phi0,
                                     this->LS->Phin, this->Vgf,
                                     this->Vgfweight, this->riemann);
  }
  this->boundaryFlux = this->tmpDistSVec;

  t0 = this->timer->getTime();

  // *** corrector step ***
  computeRKUpdate(U0, k2,2);
  this->spaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U -= 0.5 * (k1 + k2);
  this->spaceOp->applyExtrapolationToSolutionVector(U, Ubc);
  this->spaceOp->applyBCsToSolutionVector(U);

  this->timer->addFluidSolutionTime(t0);

  if(!(this->interfaceType==MultiFluidData::FSF)){
    this->spaceOp->storePreviousPrimitive(U, this->Vg, this->Phi,
                                          this->Vgf, this->Vgfweight);

    t0 = this->timer->getTime();

    computeRKUpdateLS(Phi0, p2, U0);
    this->Phi = this->Phi - 0.5 * (p1 + p2);

    this->timer->addLevelSetSolutionTime(t0);

    // Riemann overwrite on U_{n+1} using the value of Phi_{n+1}
    this->spaceOp->updatePhaseChange(this->Vg, U, this->Phi,
                                     this->LS->Phin, this->Vgf,
                                     this->Vgfweight, this->riemann);
  }
  checkSolution(U);
  this->boundaryFlux += this->tmpDistSVec;
  this->boundaryFlux *= 0.5;


}

//------------------------------------------------------------------------------

template<int dim>
void ExplicitLevelSetTsDesc<dim>::solveNLSystemTwoBlocks(DistSVec<double,dim> &U)
{
 
  solveNLEuler(U);

  if(!(this->interfaceType==MultiFluidData::FSF)){

    this->spaceOp->storePreviousPrimitive(U, this->Vg, this->Phi, 
                                          this->Vgf, this->Vgfweight);

    solveNLLevelSet(U);

    this->spaceOp->updatePhaseChange(this->Vg, U, this->Phi, 
                                     this->LS->Phin, this->Vgf, 
                                     this->Vgfweight, this->riemann);
  }

  checkSolution(U);

}

//------------------------------------------------------------------------------

template<int dim>
void ExplicitLevelSetTsDesc<dim>::solveNLEuler(DistSVec<double,dim> &U)
{
  double t0 = this->timer->getTime();

  if (timeType == ExplicitData::RUNGE_KUTTA_4) solveNLEulerRK4(U);
  else                                         solveNLEulerRK2(U);
  
  this->timer->addFluidSolutionTime(t0);

}
//------------------------------------------------------------------------------

template<int dim>
void ExplicitLevelSetTsDesc<dim>::solveNLEulerRK2(DistSVec<double,dim> &U)
{
  DistSVec<double,dim> Ubc(this->getVecInfo());
  this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);

  computeRKUpdate(U, k1,1);
  this->spaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - k1;
  this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  checkSolution(U0);
  this->boundaryFlux = this->tmpDistSVec;


  computeRKUpdate(U0, k2,2);
  this->spaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U -= 1.0/2.0 * (k1 + k2);
  this->spaceOp->applyExtrapolationToSolutionVector(U, Ubc);
  this->spaceOp->applyBCsToSolutionVector(U);
  checkSolution(U);
  this->boundaryFlux += this->tmpDistSVec;
  this->boundaryFlux *= 0.5;

}

//------------------------------------------------------------------------------

template<int dim>
void ExplicitLevelSetTsDesc<dim>::solveNLEulerRK4(DistSVec<double,dim> &U)
{

  DistSVec<double,dim> Ubc(this->getVecInfo());
  this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);

  computeRKUpdate(U, k1, 1);
  this->spaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - 0.5 * k1;
  this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  checkSolution(this->U0);
  this->boundaryFlux = this->tmpDistSVec;


  computeRKUpdate(U0, k2, 2);
  this->spaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U0 = U - 0.5 * k2;
  this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  checkSolution(U0);
  this->boundaryFlux += 2.0*this->tmpDistSVec;


  computeRKUpdate(U0, k3, 3);
  this->spaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U0 = U - k3;
  this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  checkSolution(U0);
  this->boundaryFlux += 2.0*this->tmpDistSVec;


  computeRKUpdate(U0, k4, 4);
  this->spaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U -= 1.0/6.0 * (k1 + 2.0 * (k2 + k3) + k4);
  this->spaceOp->applyExtrapolationToSolutionVector(U, Ubc);
  this->spaceOp->applyBCsToSolutionVector(U);
  checkSolution(U);
  this->boundaryFlux += this->tmpDistSVec;
  this->boundaryFlux *= 1.0/6.0;


}

//------------------------------------------------------------------------------
template<int dim>
void ExplicitLevelSetTsDesc<dim>::solveNLLevelSet(DistSVec<double,dim> &U)
{

  double t0 = this->timer->getTime();

  if (timeType == ExplicitData::RUNGE_KUTTA_4) solveNLLevelSetRK4(U);
  else                                         solveNLLevelSetRK2(U);

  this->timer->addLevelSetSolutionTime(t0);

}
//------------------------------------------------------------------------------
template<int dim>
void ExplicitLevelSetTsDesc<dim>::solveNLLevelSetRK2(DistSVec<double,dim> &U)
{

  computeRKUpdateLS(this->Phi, p1, U);
  Phi0 = this->Phi - p1;

  computeRKUpdateLS(Phi0, p2, U);
  this->Phi -= 1.0/2.0 * (p1+p2);

}
//------------------------------------------------------------------------------
template<int dim>
void ExplicitLevelSetTsDesc<dim>::solveNLLevelSetRK4(DistSVec<double,dim> &U)
{

  computeRKUpdateLS(this->Phi, p1, U);
  Phi0 = this->Phi - 0.5 * p1;

  computeRKUpdateLS(Phi0, p2, U);
  Phi0 = this->Phi - 0.5 * p2;

  computeRKUpdateLS(Phi0, p3, U);
  Phi0 = this->Phi - p3;

  computeRKUpdateLS(Phi0, p4, U);
  this->Phi -= 1.0/6.0 * (p1 + 2.0 * (p2 + p3) + p4);

}
//------------------------------------------------------------------------------

template<int dim>
void ExplicitLevelSetTsDesc<dim>::computeRKUpdate(DistSVec<double,dim>& Ulocal,
                                  DistSVec<double,dim>& dU, int it)
{
  this->tmpDistSVec = 0.0;
  this->spaceOp->applyBCsToSolutionVector(Ulocal);
  this->spaceOp->computeResidual(*this->X, *this->A, Ulocal, this->PhiV, 
                                 dU, this->riemann,it, &this->tmpDistSVec);
  // for RK2 on moving grids
  this->domain->computeVolumeChangeTerm(*this->A, *this->geoState, Ulocal, dU);
  this->timeState->multiplyByTimeStep(dU);
}

//------------------------------------------------------------------------------
template<int dim>
void ExplicitLevelSetTsDesc<dim>::computeRKUpdateLS(DistVec<double> &Philocal,
                                  DistVec<double> &dPhi, DistSVec<double,dim> &U)
{

  this->spaceOp->computeResidualLS(*this->X, *this->A, Philocal, U, dPhi);
  // for RK2 on moving grids
  this->domain->computeVolumeChangeTerm(*this->A, *this->geoState, Philocal, dPhi);
  this->timeState->multiplyByTimeStep(dPhi);

}
