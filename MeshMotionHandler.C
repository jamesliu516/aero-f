#include <stdlib.h>
#include <math.h>

#include <MeshMotionHandler.h>

#include <Domain.h>
#include <PostOperator.h>
#include <StructExc.h>
#include <MeshMotionSolver.h>

//------------------------------------------------------------------------------
// TODO variable time-step

template<int dim>
void MeshMotionHandler::computeInterfaceWork(double dt, PostOperator<dim>* postOp, 
					     DistSVec<double,3>& Xn, DistSVec<double,dim>& Un, 
					     DistSVec<double,3>& X, DistSVec<double,dim>& U, 
					     double* Wnp1)
{

  Wnp1[0] = 1./3.*Wn + 2./3.*dt*postOp->computeInterfaceWork(X, U, Pin);
  Wn = Wnp1[0];

  postOp->computeNodalForce(Xn, Un, Pin, dX);
  postOp->computeNodalForce(X, U, Pin, F);
  F = 0.5 * (dX + F);
  // F needs to be assembled in order to get the correct scalar product
  domain->assemble(domain->getVec3DPat(), F);
  Wnp1[1] = F * (Xn - X);

}

//------------------------------------------------------------------------------
// note: this initialization is not correct in case of a restart with an averaged force

template<int dim>
void AeroMeshMotionHandler::setup(int *rstrt, double *maxTime, PostOperator<dim>* postOp,
				  DistSVec<double,3> &X, DistSVec<double,dim> &U)
{

  if (strExc->getAlgorithmNumber() == 4)
    postOp->computeNodalForce(X, U, Pin, F);
  
  *rstrt = strExc->getRestartFrequency();
  *maxTime = strExc->getMaxTime();

  mms->setup(X);

}

//------------------------------------------------------------------------------

template<int dim>
void AeroMeshMotionHandler::resetOutputToStructure(PostOperator<dim>* postOp,
						   DistSVec<double,3> &X, 
						   DistSVec<double,dim> &U)
{

  if (forceComputation == AeroelasticData::AVERAGED) {
    *Favg = 0.0;
    if (timeIntegrator == IMPLICIT_SECOND_ORDER) 
      postOp->computeNodalForce(X, U, Pin, *Fn);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void AeroMeshMotionHandler::updateOutputToStructure(double dt, double dtLeft,
						    PostOperator<dim>* postOp,
						    DistSVec<double,3> &X, 
						    DistSVec<double,dim> &U)
{

  if (forceComputation == AeroelasticData::AVERAGED) {
    switch (timeIntegrator) {
    case IMPLICIT_FIRST_ORDER : {
      postOp->computeNodalForce(X, U, Pin, F);
      break;
    }
    case IMPLICIT_SECOND_ORDER : {
      postOp->computeNodalForce(X, U, Pin, *Fnp1);
      F = 0.5 * (*Fn + *Fnp1);
      *Fn = *Fnp1;
      break;
    }
    }
    *Favg += dt * F;
  }

  if (dtLeft == 0.0) {
    if (forceComputation == AeroelasticData::AVERAGED)
      F = (1.0 / strExc->getTimeStep()) * (*Favg);
    else if (forceComputation == AeroelasticData::LAST)
      postOp->computeNodalForce(X, U, Pin, F);
    else if (forceComputation == AeroelasticData::LAST_KRIS) {
      postOp->computeNodalForce(X, U, Pin, F);
      F = 1./3.*(*Fn) + 2./3.*F;
      *Fn = F;
    }
  }

}

//------------------------------------------------------------------------------
