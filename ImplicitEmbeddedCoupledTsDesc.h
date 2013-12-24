#pragma once

#include <ImplicitEmbeddedTsDesc.h>

#include "IoData.h"
#include "GeoSource.h"
#include "Domain.h"
#include "LevelSet.h"
#include "DistTimeState.h"

#include <MatVecProd.h>
#include <KspSolver.h>
#include <SpaceOperator.h>
#include <NewtonSolver.h>
#include <DistEmbeddedVector.h>

#ifdef TYPE_PREC
#define PrecScalar TYPE_PREC
#else
#define PrecScalar double
#endif

//------------------------------------------------------------------------------

template<int dim>
class ImplicitEmbeddedCoupledTsDesc : public ImplicitEmbeddedTsDesc<dim> {

  MatVecProd<dim,dim> *mvp;
  KspPrec<dim> *pc;
  KspSolver<DistEmbeddedVec<double,dim>, MatVecProd<dim,dim>, KspPrec<dim>, Communicator> *ksp;
	
public:
  
  ImplicitEmbeddedCoupledTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom);

  ~ImplicitEmbeddedCoupledTsDesc();

  void computeJacobian(int it, DistSVec<double,dim> &Q,
		       DistSVec<double,dim> &F);

  void setOperators(DistSVec<double,dim> &Q);

  int solveLinearSystem(int it, DistSVec<double,dim> &b,
			DistSVec<double,dim> &dQ);
};


#ifdef TEMPLATE_FIX
#include <ImplicitEmbeddedCoupledTsDesc.C>
#endif
