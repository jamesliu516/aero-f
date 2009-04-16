
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
ExplicitStructLevelSetTsDesc<dim>::
ExplicitStructLevelSetTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  StructLevelSetTsDesc<dim>(ioData, geoSource, dom),
  k1(this->getVecInfo()), k2(this->getVecInfo()),
  k3(this->getVecInfo()), k4(this->getVecInfo()),
//  p1(this->getVecInfo()), p2(this->getVecInfo()),
//  p3(this->getVecInfo()), p4(this->getVecInfo()),
  U0(this->getVecInfo())//, Phi0(this->getVecInfo())
{
  timeType = ioData.ts.expl.type;
  this->mmh = 0;
}

//------------------------------------------------------------------------------

template<int dim>
ExplicitStructLevelSetTsDesc<dim>::~ExplicitStructLevelSetTsDesc()
{
}

//------------------------------------------------------------------------------

template<int dim>
int ExplicitStructLevelSetTsDesc<dim>::solveNonLinearSystem(DistSVec<double,dim>& U)
// the key player.
{solveNLSystemOneBlock(U);} //let's just do a Forward Euler for now.

//-----------------------------------------------------------------------------

template<int dim>
void ExplicitStructLevelSetTsDesc<dim>::solveNLSystemOneBlock(DistSVec<double,dim> &U)
{solveNLAllFE(U);} // Forward Euler.
//------------------------------------------------------------------------------

template<int dim>
void ExplicitStructLevelSetTsDesc<dim>::solveNLAllFE(DistSVec<double,dim> &U)
{

  double t0 = this->timer->getTime();

  DistSVec<double,dim> Ubc(this->getVecInfo());
  computeRKUpdate(U, k1, 1);  //TODO: changes happen here.
  this->spaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - k1;
  this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);

  this->timer->addFluidSolutionTime(t0);
  U = U0;
  checkSolution(U);
}


template<int dim>
void ExplicitStructLevelSetTsDesc<dim>::computeRKUpdate(DistSVec<double,dim>& Ulocal,
                                  DistSVec<double,dim>& dU, int it)
{
  this->spaceOp->applyBCsToSolutionVector(Ulocal);
  this->distLSS->clearTotalForce();

  this->spaceOp->computeResidual(*this->X, *this->A, Ulocal, *this->Wstarij, *this->Wstarji, this->distLSS,
                                 dU, this->riemann,it);
  this->distLSS->getTotalForce(this->pressureRef);
  // for RK2 on moving grids
//  this->domain->computeVolumeChangeTerm(*this->A, *this->geoState, Ulocal, dU);
  this->timeState->multiplyByTimeStep(dU);
}

//------------------------------------------------------------------------------
