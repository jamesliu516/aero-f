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

  //if(ioData.ts.expl.type == ExplicitData::RUNGE_KUTTA_4) RK4 = true;
  //else RK4 = false;
  timeIntegrationType = ioData.ts.expl.type;
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

  if(timeIntegrationType==ExplicitData::RUNGE_KUTTA_4 ||
     timeIntegrationType==ExplicitData::RUNGE_KUTTA_2 ){

  /* resolution of the non linear system (U,rho*phi)
  ** using a staggered scheme, first we solve U,
  ** then using the new value of U, we solve rho*phi.
  **
  ** We use a GFMP scheme, ie the following steps are done:
  ** 1- solve for U -> Utilde using rho*phi
  ** 2- Utilde -> V (conservative to primitive) using rho*phi
  ** 3- solve for rho*phi using Utilde
  ** 4- V -> U using the new rho*phi
  */

    solveNonLinearSystemEuler(U);

    if(!(this->interfaceType==MultiFluidData::FSF)){
      this->spaceOp->storePreviousPrimitive(U, this->Vg, this->Phi, 
                                            this->Vgf, this->Vgfweight);

      solveNonLinearSystemLevelSet(U);

      this->spaceOp->updatePhaseChange(this->Vg, U, this->Phi, 
                                       this->LS->Phin, this->Vgf, 
                                       this->Vgfweight, this->riemann);
    }
  }else if(timeIntegrationType==ExplicitData::FORWARD_EULER){

    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);
    computeRKUpdate(U, k1,1);
    U0 = U - k1;

    if(!(this->interfaceType==MultiFluidData::FSF)){
      this->spaceOp->storePreviousPrimitive(U0, this->Vg, this->Phi,
                                            this->Vgf, this->Vgfweight);

      computeRKUpdateLS(this->Phi, p1, U);
      this->Phi = this->Phi - p1;
      U = U0;

      // Riemann overwrite using the value of Phi_{n+1}
      this->spaceOp->updatePhaseChange(this->Vg, U, this->Phi,
                                       this->LS->Phin, this->Vgf,
                                       this->Vgfweight, this->riemann);
    }

  }else if(timeIntegrationType==ExplicitData::RK2_TALLEC){
    DistSVec<double,dim> Ubc(this->getVecInfo());
    // prediction step
    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);
    computeRKUpdate(U, k1,1);
    this->spaceOp->getExtrapolationValue(U, Ubc, *this->X);
    U0 = U - k1;
    this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);

    if(!(this->interfaceType==MultiFluidData::FSF)){
      this->spaceOp->storePreviousPrimitive(U0, this->Vg, this->Phi,
                                            this->Vgf, this->Vgfweight);

      computeRKUpdateLS(this->Phi, p1, U);
      Phi0 = this->Phi - p1;

      // Riemann overwrite using the value of Phi0
      this->spaceOp->updatePhaseChange(this->Vg, U0, Phi0,
                                       this->LS->Phin, this->Vgf,
                                       this->Vgfweight, this->riemann);
    }

    checkSolution(U0);

    // correction step
    this->LS->conservativeToPrimitive(Phi0,this->PhiV,U0);
    computeRKUpdate(U0, k2,2);
    this->spaceOp->getExtrapolationValue(U0, Ubc, *this->X);
    U -= 1.0/2.0 * (k1 + k2);
    this->spaceOp->applyExtrapolationToSolutionVector(U, Ubc);
    this->spaceOp->applyBCsToSolutionVector(U);

    if(!(this->interfaceType==MultiFluidData::FSF)){
      this->spaceOp->storePreviousPrimitive(U, this->Vg, Phi0,
                                            this->Vgf, this->Vgfweight);

      computeRKUpdateLS(Phi0, p2, U0);
      this->Phi -= 1.0/2.0 * (p1+p2);

      // Riemann overwrite using the value of Phi (at n+1)
      this->spaceOp->updatePhaseChange(this->Vg, U, this->Phi,
                                       this->LS->Phin, this->Vgf,
                                       this->Vgfweight, this->riemann);
      // or should it be the following?
      //this->spaceOp->updatePhaseChange(this->Vg, U, this->Phi,
      //                                 Phi0, this->Vgf,
      //                                 this->Vgfweight, this->riemann);
    }
  }

  checkSolution(U);

  return 1;
}

//------------------------------------------------------------------------------

template<int dim>
void ExplicitLevelSetTsDesc<dim>::solveNonLinearSystemEuler(DistSVec<double,dim> &U)
{
  double t0 = this->timer->getTime();

  if (timeIntegrationType==ExplicitData::RUNGE_KUTTA_4){
    DistSVec<double,dim> Ubc(this->getVecInfo());
    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);

    computeRKUpdate(U, k1, 1);
    this->spaceOp->getExtrapolationValue(U, Ubc, *this->X);
    U0 = U - 0.5 * k1;
    this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
    checkSolution(this->U0);


    computeRKUpdate(U0, k2, 2);
    this->spaceOp->getExtrapolationValue(U0, Ubc, *this->X);
    U0 = U - 0.5 * k2;
    this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
    checkSolution(U0);


    computeRKUpdate(U0, k3, 3);
    this->spaceOp->getExtrapolationValue(U0, Ubc, *this->X);
    U0 = U - k3;
    this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
    checkSolution(U0);


    computeRKUpdate(U0, k4, 4);
    this->spaceOp->getExtrapolationValue(U0, Ubc, *this->X);
    U -= 1.0/6.0 * (k1 + 2.0 * (k2 + k3) + k4);
    this->spaceOp->applyExtrapolationToSolutionVector(U, Ubc);
    checkSolution(U);


    this->spaceOp->applyBCsToSolutionVector(U);

  } else{

    DistSVec<double,dim> Ubc(this->getVecInfo());
    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);

    computeRKUpdate(U, k1,1);
    this->spaceOp->getExtrapolationValue(U, Ubc, *this->X);

    U0 = U - k1;
    this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
    checkSolution(U0);
    computeRKUpdate(U0, k2,2);
    this->spaceOp->getExtrapolationValue(U0, Ubc, *this->X);

    U -= 1.0/2.0 * (k1 + k2);
    this->spaceOp->applyExtrapolationToSolutionVector(U, Ubc);
    this->spaceOp->applyBCsToSolutionVector(U);
    checkSolution(U);
  
  }
  this->timer->addFluidSolutionTime(t0);

}
//------------------------------------------------------------------------------
template<int dim>
void ExplicitLevelSetTsDesc<dim>::solveNonLinearSystemLevelSet(DistSVec<double,dim> &U)
{

  double t0 = this->timer->getTime();

  if (timeIntegrationType==ExplicitData::RUNGE_KUTTA_4){
    computeRKUpdateLS(this->Phi, p1, U);
    Phi0 = this->Phi - 0.5 * p1;

    computeRKUpdateLS(Phi0, p2, U);
    Phi0 = this->Phi - 0.5 * p2;

    computeRKUpdateLS(Phi0, p3, U);
    Phi0 = this->Phi - p3;

    computeRKUpdateLS(Phi0, p4, U);
    this->Phi -= 1.0/6.0 * (p1 + 2.0 * (p2 + p3) + p4);

  }else{

    computeRKUpdateLS(this->Phi, p1, U);
    Phi0 = this->Phi - p1;

    computeRKUpdateLS(Phi0, p2, U);
    this->Phi -= 1.0/2.0 * (p1+p2);

  }

  this->timer->addLevelSetSolutionTime(t0);

}
//------------------------------------------------------------------------------

template<int dim>
void ExplicitLevelSetTsDesc<dim>::computeRKUpdate(DistSVec<double,dim>& Ulocal,
				DistSVec<double,dim>& dU, int it)
{
  this->spaceOp->applyBCsToSolutionVector(Ulocal);
  this->spaceOp->computeResidual(*this->X, *this->A, Ulocal, this->PhiV, dU, this->riemann,it);
  this->domain->computeVolumeChangeTerm(*this->A, *this->geoState, Ulocal, dU);
  this->timeState->multiplyByTimeStep(dU);
}

//------------------------------------------------------------------------------
template<int dim>
void ExplicitLevelSetTsDesc<dim>::computeRKUpdateLS(DistVec<double> &Philocal,
 				    DistVec<double> &dPhi, DistSVec<double,dim> &U)
{

  this->spaceOp->computeResidualLS(*this->X, *this->A, Philocal, U, dPhi);
  this->domain->computeVolumeChangeTerm(*this->A, *this->geoState, Philocal, dPhi);
  this->timeState->multiplyByTimeStep(dPhi);

}
