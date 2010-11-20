#ifndef _IMPLICIT_COUPLED_TS_DESC_H_
#define _IMPLICIT_COUPLED_TS_DESC_H_

#include <ImplicitTsDesc.h>
#include <KspPrec.h>
#include <IoData.h>
#include <Domain.h>

template<int dim, int neq> class MatVecProd;

#ifndef _KSPSLVR_TMPL_
#define _KSPSLVR_TMPL_
template<class VecType, class MvpOp, class PrecOp, class IoOp, class ScalarT = double> class KspSolver;
#endif

//------------------------------------------------------------------------------

template<int dim>
class ImplicitCoupledTsDesc : public ImplicitTsDesc<dim> {

protected:
  MatVecProd<dim,dim> *mvp;
  KspPrec<dim> *pc;
  KspSolver<DistSVec<double,dim>, MatVecProd<dim,dim>, KspPrec<dim>, Communicator> *ksp;
  
#ifdef MVP_CHECK
  MatVecProd<dim,dim> *mvpfd1;
#endif

public:

  ImplicitCoupledTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitCoupledTsDesc();

  void computeJacobian(int, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void setOperators(DistSVec<double,dim> &);
  int solveLinearSystem(int, DistSVec<double,dim> &, DistSVec<double,dim> &);
  
// Included (MB)
  void rstVarImplicitCoupledTsDesc(IoData &);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitCoupledTsDesc.C>
#endif

#endif
