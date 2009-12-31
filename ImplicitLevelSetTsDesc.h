#ifndef _IMPLICIT_LEVELSET_TS_DESC_H_
#define _IMPLICIT_LEVELSET_TS_DESC_H_

#include <LevelSetTsDesc.h>

#include <IoData.h>
#include <KspPrec.h>
#include <Domain.h>

struct DistInfo;

class GeoSource;
class LevelSet; 
template<class Scalar, int dim> class DistSVec;
template<int dim, int neq> class MatVecProd;

                                                                                                                
#ifndef _KSPSLVR_TMPL_
#define _KSPSLVR_TMPL_
template<class VecType, class MvpOp, class PrecOp, class IoOp, class ScalarT = double> class KspSolver;
#endif


//------------------------------------------------------------------------

template<int dim>
class ImplicitLevelSetTsDesc : public LevelSetTsDesc<dim> {

 protected: 
  DistSVec<bool,2> *tag;

  MatVecProd<dim,dim> *mvp;
  KspPrec<dim> *pc;
  KspSolver<DistSVec<double,dim>, MatVecProd<dim,dim>, KspPrec<dim>, Communicator> *ksp;
	
  MatVecProd<dim,1> *mvpLS;
  KspPrec<1> *pcLS;
  KspSolver<DistVec<double>, MatVecProd<dim,1>, KspPrec<1>, Communicator> *kspLS;

  NewtonSolver<ImplicitLevelSetTsDesc<dim> > *ns;

  int failSafeNewton;
  int maxItsNewton;
  double epsNewton;
  double epsAbsResNewton, epsAbsIncNewton;


 public:
  ImplicitLevelSetTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitLevelSetTsDesc();
  
  int getMaxItsNewton() const { return maxItsNewton; }
  double getEpsNewton() const { return epsNewton; }
  double getEpsAbsResNewton() const { return epsAbsResNewton; }
  double getEpsAbsIncNewton() const { return epsAbsIncNewton; }

  //-- functions for solving Euler equations
  int solveNonLinearSystem(DistSVec<double,dim> &);
  void computeFunction(int, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void recomputeFunction(DistSVec<double,dim> &, DistSVec<double,dim> &);
  
  int checkFailSafe(DistSVec<double,dim>&);
  void resetFixesTag();

  void computeJacobian(int, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void setOperators(DistSVec<double,dim> &);
  int solveLinearSystem(int, DistSVec<double,dim> &, DistSVec<double,dim> &);

  //-- new functions for solving LevelSet equation
  void computeFunctionLS(int, DistSVec<double,dim> &,DistVec<double> &);
  void computeJacobianLS(int, DistSVec<double,dim> &,DistVec<double> &);
  int solveLinearSystemLS(int, DistVec<double> &, DistVec<double> &);


 protected:
  template<class Scalar, int neq>
  KspPrec<neq> *createPreconditioner(PcData &, Domain *);

  template<int neq>
  KspSolver<DistSVec<double,neq>, MatVecProd<dim,neq>, KspPrec<neq>,
  Communicator> *createKrylovSolver(const DistInfo &, KspData &, MatVecProd<dim,neq> *,
                                    KspPrec<neq> *, Communicator *);


  template<int neq>
  KspSolver<DistVec<double>, MatVecProd<dim,neq>, KspPrec<neq,double>,
  Communicator> *createKrylovSolverLS(const DistInfo &, KspData &, MatVecProd<dim,neq> *,
                                        KspPrec<neq,double> *, Communicator *);
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitLevelSetTsDesc.C>
#endif

#endif
