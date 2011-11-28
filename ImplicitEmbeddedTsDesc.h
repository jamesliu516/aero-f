#pragma once

#include <EmbeddedTsDesc.h>

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
class ImplicitEmbeddedTsDesc : public EmbeddedTsDesc<dim> {

  DistSVec<bool,2> *tag;

  MatVecProd<dim,dim> *mvp;
  KspPrec<dim> *pc;
  KspSolver<DistEmbeddedVec<double,dim>, MatVecProd<dim,dim>, KspPrec<dim>, Communicator> *ksp;
	
  NewtonSolver<ImplicitEmbeddedTsDesc<dim> > *ns;

  int failSafeNewton;
  int maxItsNewton;
  double epsNewton;
  double epsAbsResNewton, epsAbsIncNewton;

  template <int neq>
  KspPrec<neq> *createPreconditioner(PcData &pcdata, Domain *dom);
 
  DistEmbeddedVec<double,dim> embeddedU,embeddedB,embeddeddQ;

public:
  
  ImplicitEmbeddedTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom);

  ~ImplicitEmbeddedTsDesc();

  template<int neq, class MatVecProdOp>
  KspSolver<DistEmbeddedVec<double,neq>, MatVecProdOp, KspPrec<neq>, Communicator> *
  createKrylovSolver(const DistInfo &info, KspData &kspdata,
                     MatVecProdOp *_mvp, KspPrec<neq> *_pc,
                     Communicator *_com);

  int commonPart(DistSVec<double,dim> &U);

  int solveNonLinearSystem(DistSVec<double,dim> &U, int);
  
  void computeFunction(int it, DistSVec<double,dim> &Q,
                       DistSVec<double,dim> &F);

  void recomputeFunction(DistSVec<double,dim> &Q,
                         DistSVec<double,dim> &rhs);

  void doErrorEstimation(DistSVec<double,dim> &U);

  int checkFailSafe(DistSVec<double,dim>& U);

  void resetFixesTag();

  void computeJacobian(int it, DistSVec<double,dim> &Q,
		       DistSVec<double,dim> &F);

  void setOperators(DistSVec<double,dim> &Q);

  int solveLinearSystem(int it, DistSVec<double,dim> &b,
			DistSVec<double,dim> &dQ);

  int getMaxItsNewton() const { return maxItsNewton; }
  double getEpsNewton() const { return epsNewton; }
  double getEpsAbsResNewton() const { return epsAbsResNewton; }
  double getEpsAbsIncNewton() const { return epsAbsIncNewton; }
};


#ifdef TEMPLATE_FIX
#include <ImplicitEmbeddedTsDesc.C>
#endif
