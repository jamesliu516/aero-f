#include <ExplicitLevelSetTsDesc.h>

#include "IoData.h"
#include "Domain.h"
#include <GeoSource.h>
#include <DistTimeState.h>
#include <SpaceOperator.h>
#include "LevelSet.h"

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
ExplicitLevelSetTsDesc<dim,dimLS>::
ExplicitLevelSetTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  LevelSetTsDesc<dim,dimLS>(ioData, geoSource, dom), 
  k1(this->getVecInfo()), k2(this->getVecInfo()), 
  k3(this->getVecInfo()), k4(this->getVecInfo()), 
  p1(this->getVecInfo()), p2(this->getVecInfo()), 
  p3(this->getVecInfo()), p4(this->getVecInfo()), 
  U0(this->getVecInfo()), Phi0(this->getVecInfo()),
  ratioTimesU(this->getVecInfo()), 
  ratioTimesPhi(this->getVecInfo()),
  fluidId0(this->getVecInfo())
{
  this->mmh = this->createMeshMotionHandler(ioData, geoSource, 0);

  timeType = ioData.ts.expl.type;
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
ExplicitLevelSetTsDesc<dim,dimLS>::~ExplicitLevelSetTsDesc()
{
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
int ExplicitLevelSetTsDesc<dim,dimLS>::solveNonLinearSystem(DistSVec<double,dim> &U, int)
{

  if(timeType == ExplicitData::FORWARD_EULER || 
     timeType == ExplicitData::ONE_BLOCK_RK2 ||
     timeType == ExplicitData::ONE_BLOCK_RK2bis )
    solveNLSystemOneBlock(U);
  else if(timeType == ExplicitData::RUNGE_KUTTA_2 ||
          timeType == ExplicitData::RUNGE_KUTTA_4 )
    solveNLSystemTwoBlocks(U);

  return 1;

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLSystemOneBlock(DistSVec<double,dim> &U)
{

  if(timeType == ExplicitData::FORWARD_EULER)
    solveNLAllFE(U);
  else if(timeType == ExplicitData::ONE_BLOCK_RK2)
    solveNLAllRK2(U);
  else if(timeType == ExplicitData::ONE_BLOCK_RK2bis)
    solveNLAllRK2bis(U);

}
//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLAllFE(DistSVec<double,dim> &U)
{

  double t0 = this->timer->getTime();

  DistSVec<double,dim> Ubc(this->getVecInfo());
  if (this->lsMethod == 0) {
    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);
  }

  computeRKUpdate(U, k1,1);
  this->multiPhaseSpaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - k1;
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U0, Ubc);

  this->timer->addFluidSolutionTime(t0);

  if(!(this->interfaceType==MultiFluidData::FSF)){
    this->varFcn->conservativeToPrimitive(U0,this->V0,this->fluidSelector.fluidId);
    this->riemann->storePreviousPrimitive(this->V0, *this->fluidSelector.fluidId, *this->X);
    t0 = this->timer->getTime();

    computeRKUpdateLS(this->Phi, *this->fluidSelector.fluidId, p1, U);
    this->Phi = this->Phi - p1;
    this->riemann->avoidNewPhaseCreation(this->Phi, this->LS->Phin);
    (this->fluidSelector).getFluidId(this->Phi); //update fluidId accordingly

    this->timer->addLevelSetSolutionTime(t0);

    // Riemann overwrite using the value of Phi_{n+1}
    this->riemann->updatePhaseChange(this->V0, *this->fluidSelector.fluidId, *this->fluidSelector.fluidIdn);
    this->varFcn->primitiveToConservative(this->V0,U,this->fluidSelector.fluidId);
  }else U = U0;

  checkSolution(U);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLAllRK2(DistSVec<double,dim> &U)
{

  double t0 = this->timer->getTime();
  this->domain->computePrdtWCtrlVolRatio(ratioTimesU, U, *this->A, *this->geoState);
  this->domain->computePrdtPhiCtrlVolRatio(ratioTimesPhi, this->Phi, *this->A, *this->geoState);

  DistSVec<double,dim> Ubc(this->getVecInfo());
  if (this->lsMethod == 0) {
    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);
  }
  // *** prediction step ***
  computeRKUpdate(U, k1,1);
  this->multiPhaseSpaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = ratioTimesU - k1;
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U0, Ubc);

  this->timer->addFluidSolutionTime(t0);

  if(!(this->interfaceType==MultiFluidData::FSF)){
    this->varFcn->conservativeToPrimitive(U0,this->V0,this->fluidSelector.fluidId);
    this->riemann->storePreviousPrimitive(this->V0, *this->fluidSelector.fluidId, *this->X);

    t0 = this->timer->getTime();

    computeRKUpdateLS(this->Phi, *this->fluidSelector.fluidId, p1, U);
    Phi0 = ratioTimesPhi - p1;
    this->fluidSelector.getFluidId(fluidId0,Phi0);

    this->timer->addLevelSetSolutionTime(t0);

    // Riemann overwrite on U0 using the value of Phi0
  }

  t0 = this->timer->getTime();

  // *** corrector step ***
  computeRKUpdate(U0, k2,2);
  this->multiPhaseSpaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U = ratioTimesU - 0.5 * (k1 + k2);
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U, Ubc);
  this->multiPhaseSpaceOp->applyBCsToSolutionVector(U);

  this->timer->addFluidSolutionTime(t0);

  if(!(this->interfaceType==MultiFluidData::FSF)){

    t0 = this->timer->getTime();

    computeRKUpdateLS(Phi0, fluidId0, p2, U0);
    this->Phi = ratioTimesPhi - 0.5 * (p1 + p2);
    this->riemann->avoidNewPhaseCreation(this->Phi, this->LS->Phin);
    (this->fluidSelector).getFluidId(this->Phi);

    this->timer->addLevelSetSolutionTime(t0);

    // Riemann overwrite on U_{n+1} using the value of Phi_{n+1}
    this->riemann->updatePhaseChange(this->V0, *this->fluidSelector.fluidId, *this->fluidSelector.fluidIdn);
    this->varFcn->primitiveToConservative(this->V0,U,this->fluidSelector.fluidId);
  }
  checkSolution(U);

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLAllRK2bis(DistSVec<double,dim> &U)
{

  double t0 = this->timer->getTime();

  DistSVec<double,dim> Ubc(this->getVecInfo());
  if (this->lsMethod == 0) {
    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);
  }
  // *** prediction step ***
  computeRKUpdate(U, k1,1);
  this->multiPhaseSpaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - k1;
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U0, Ubc);

  this->timer->addFluidSolutionTime(t0);

  if(!(this->interfaceType==MultiFluidData::FSF)){
    this->varFcn->conservativeToPrimitive(U0,this->V0,this->fluidSelector.fluidId);
    this->riemann->storePreviousPrimitive(this->V0, *this->fluidSelector.fluidId, *this->X);

    t0 = this->timer->getTime();

    computeRKUpdateLS(this->Phi, *this->fluidSelector.fluidId, p1, U0);
    Phi0 = this->Phi - p1;
    
    this->fluidSelector.getFluidId(fluidId0,Phi0);

    this->timer->addLevelSetSolutionTime(t0);

    // Riemann overwrite on U0 using the value of Phi0
    //this->multiPhaseSpaceOp->updatePhaseChange(this->Vg, U0, fluidId0,
    //                                 this->fluidSelector.fluidIdn, this->Vgf,
    //                                 this->Vgfweight, this->riemann);
    //this->riemann->updatePhaseChange();
  }

  t0 = this->timer->getTime();

  // *** corrector step ***
  computeRKUpdate(U0, k2,2);
  this->multiPhaseSpaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U -= 0.5 * (k1 + k2);
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U, Ubc);
  this->multiPhaseSpaceOp->applyBCsToSolutionVector(U);

  this->timer->addFluidSolutionTime(t0);

  if(!(this->interfaceType==MultiFluidData::FSF)){
    //this->multiPhaseSpaceOp->storePreviousPrimitive(U, this->Vg, this->Phi,
    //                                      this->Vgf, this->Vgfweight);
    //this->riemann->storePreviousPrimitive();

    t0 = this->timer->getTime();

    computeRKUpdateLS(Phi0, fluidId0, p2, U);
    this->Phi = this->Phi - 0.5 * (p1 + p2);
    this->riemann->avoidNewPhaseCreation(this->Phi, this->LS->Phin);
    (this->fluidSelector).getFluidId(this->Phi);

    this->timer->addLevelSetSolutionTime(t0);

    // Riemann overwrite on U_{n+1} using the value of Phi_{n+1}
    this->riemann->updatePhaseChange(this->V0, *this->fluidSelector.fluidId, *this->fluidSelector.fluidIdn);
    this->varFcn->primitiveToConservative(this->V0,U,this->fluidSelector.fluidId);
  }
  checkSolution(U);

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLSystemTwoBlocks(DistSVec<double,dim> &U)
{
 
  solveNLEuler(U);

  if(!(this->interfaceType==MultiFluidData::FSF)){

    this->varFcn->conservativeToPrimitive(U,this->V0,this->fluidSelector.fluidId);
    this->riemann->storePreviousPrimitive(this->V0, *this->fluidSelector.fluidId, *this->X);

    solveNLLevelSet(U);

    this->riemann->updatePhaseChange(this->V0, *this->fluidSelector.fluidId, *this->fluidSelector.fluidIdn);
    this->varFcn->primitiveToConservative(this->V0,U,this->fluidSelector.fluidId);
  }

  checkSolution(U);

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLEuler(DistSVec<double,dim> &U)
{
  double t0 = this->timer->getTime();

  if (timeType == ExplicitData::RUNGE_KUTTA_4) solveNLEulerRK4(U);
  else                                         solveNLEulerRK2(U);
  
  this->timer->addFluidSolutionTime(t0);

}
//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLEulerRK2(DistSVec<double,dim> &U)
{
  this->domain->computePrdtWCtrlVolRatio(ratioTimesU, U, *this->A, *this->geoState);
  DistSVec<double,dim> Ubc(this->getVecInfo());
  if (this->lsMethod == 0) {
    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);
  }

  computeRKUpdate(U, k1,1);
  this->multiPhaseSpaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = ratioTimesU - k1;
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  checkSolution(U0);

  computeRKUpdate(U0, k2,2);
  this->multiPhaseSpaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U = ratioTimesU - 1.0/2.0 * (k1 + k2);
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U, Ubc);
  this->multiPhaseSpaceOp->applyBCsToSolutionVector(U);
  checkSolution(U);
}
//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLEulerRK4(DistSVec<double,dim> &U)
{

  DistSVec<double,dim> Ubc(this->getVecInfo());
  if (this->lsMethod == 0) {
    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);
  }

  computeRKUpdate(U, k1, 1);
  this->multiPhaseSpaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - 0.5 * k1;
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  checkSolution(this->U0);

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
  U -= 1.0/6.0 * (k1 + 2.0 * (k2 + k3) + k4);
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U, Ubc);
  this->multiPhaseSpaceOp->applyBCsToSolutionVector(U);
  checkSolution(U);

}

//------------------------------------------------------------------------------
template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLLevelSet(DistSVec<double,dim> &U)
{

  double t0 = this->timer->getTime();

  if (timeType == ExplicitData::RUNGE_KUTTA_4) solveNLLevelSetRK4(U);
  else                                         solveNLLevelSetRK2(U);

  this->timer->addLevelSetSolutionTime(t0);

}
//------------------------------------------------------------------------------
template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLLevelSetRK2(DistSVec<double,dim> &U)
{

  this->domain->computePrdtPhiCtrlVolRatio(ratioTimesPhi, this->Phi, *this->A, *this->geoState);
  computeRKUpdateLS(this->Phi, *this->fluidSelector.fluidId, p1, U);
  Phi0 = ratioTimesPhi - p1;
  this->riemann->avoidNewPhaseCreation(this->Phi0, this->LS->Phin);
  this->fluidSelector.getFluidId(fluidId0,Phi0);

  computeRKUpdateLS(Phi0, fluidId0, p2, U);
  this->Phi = ratioTimesPhi - 1.0/2.0 * (p1+p2);
  this->riemann->avoidNewPhaseCreation(this->Phi, this->LS->Phin);
  (this->fluidSelector).getFluidId(this->Phi);
}
//------------------------------------------------------------------------------
template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLLevelSetRK4(DistSVec<double,dim> &U)
{
  computeRKUpdateLS(this->Phi, *this->fluidSelector.fluidId, p1, U);
  Phi0 = this->Phi - 0.5 * p1;
  this->riemann->avoidNewPhaseCreation(this->Phi0, this->LS->Phin);
  this->fluidSelector.getFluidId(fluidId0,Phi0);

  computeRKUpdateLS(Phi0, fluidId0, p2, U);
  Phi0 = this->Phi - 0.5 * p2;
  this->riemann->avoidNewPhaseCreation(this->Phi0, this->LS->Phin);
  this->fluidSelector.getFluidId(fluidId0,Phi0);

  computeRKUpdateLS(Phi0, fluidId0, p3, U);
  Phi0 = this->Phi - p3;
  this->riemann->avoidNewPhaseCreation(this->Phi0, this->LS->Phin);
  this->fluidSelector.getFluidId(fluidId0,Phi0);

  computeRKUpdateLS(Phi0, fluidId0, p4, U);
  this->Phi -= 1.0/6.0 * (p1 + 2.0 * (p2 + p3) + p4);
  this->riemann->avoidNewPhaseCreation(this->Phi, this->LS->Phin);
  (this->fluidSelector).getFluidId(this->Phi);
}
//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::computeRKUpdate(DistSVec<double,dim>& Ulocal,
							DistSVec<double,dim>& dU, int it)
{
  this->multiPhaseSpaceOp->applyBCsToSolutionVector(Ulocal);

  DistSVec<double,dimLS>& locphi = this->Phi;
  if (this->lsMethod == 0)
    locphi = this->PhiV;

  if (this->interfaceOrder == 2) {
    this->multiPhaseSpaceOp->findCutCells(locphi,
					  this->cutCellStatus,
					  *this->fluidSelector.fluidId,
					  Ulocal,
					  *this->X);
  }

  this->multiPhaseSpaceOp->computeResidual(*this->X, *this->A, Ulocal, locphi, this->fluidSelector, 
					   dU, this->riemann,it);
                                 //Q: why send PhiV?
                                 //A: Riemann solver needs gradPhi.
  // for RK2 on moving grids
  this->timeState->multiplyByTimeStep(dU);
}

//------------------------------------------------------------------------------
template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::computeRKUpdateLS(DistSVec<double,dimLS> &Philocal,
                                  DistVec<int> &localFluidId,
                                  DistSVec<double,dimLS> &dPhi, DistSVec<double,dim> &U)
{

  this->multiPhaseSpaceOp->computeResidualLS(*this->X, *this->A, Philocal, localFluidId, U, dPhi);
  // for RK2 on moving grids
  this->timeState->multiplyByTimeStep(dPhi);
  this->LS->checkTrueLevelSetUpdate(dPhi);

}
