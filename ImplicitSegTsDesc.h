#ifndef _IMPLICIT_SEG_TS_DESC_H_
#define _IMPLICIT_SEG_TS_DESC_H_

#include <ImplicitTsDesc.h>
#include <KspPrec.h>

template<int dim, int neq> class MatVecProd;

#ifndef _KSPSLVR_TMPL_
#define _KSPSLVR_TMPL_
template<class VecType, class MvpOp, class PrecOp, class IoOp, class ScalarT = double> class KspSolver;
#endif


//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
class ImplicitSegTsDesc : public ImplicitTsDesc<dim> {

  DistSVec<double,neq1> b1, dQ1;
  DistSVec<double,neq2> b2, dQ2;

  SpaceOperator<dim> *spaceOp1;
  SpaceOperator<dim> *spaceOp2;

  MatVecProd<dim,neq1> *mvp1;
  MatVecProd<dim,neq2> *mvp2;
  MatVecProd<dim,dim> *mvpfd;

  KspPrec<neq1> *pc1;
  KspPrec<neq2> *pc2;

  KspSolver<DistSVec<double,neq1>, MatVecProd<dim,neq1>, KspPrec<neq1>, Communicator> *ksp1;
  KspSolver<DistSVec<double,neq2>, MatVecProd<dim,neq2>, KspPrec<neq2>, Communicator> *ksp2;

private:

  SpaceOperator<dim> *createSpaceOperator1(IoData &, SpaceOperator<dim> *);
  SpaceOperator<dim> *createSpaceOperator2(IoData &, SpaceOperator<dim> *);

  template<int neq>
  void setOperator(MatVecProd<dim,neq> *, KspPrec<neq> *, DistSVec<double,dim> &, SpaceOperator<dim> *);

public:

  ImplicitSegTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitSegTsDesc();

  void computeJacobian(int, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void computeJacobianLS(int, DistVec<double> &, DistVec<double> &,  DistVec<double> &, DistVec<double> &) {};
  void setOperators(DistSVec<double,dim> &);
  int solveLinearSystem(int, DistSVec<double,dim> &, DistSVec<double,dim> &);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitSegTsDesc.C>
#endif

#endif
