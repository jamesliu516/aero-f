#pragma once

#include <IoData.h>
#include <TsDesc.h>
#include <KspPrec.h>

struct DistInfo;

class GeoSource;
class Domain;
class Communicator;

template<class Scalar, int dim> class DistSVec;
template<int dim, int neq> class MatVecProd;
template<class ProblemDescriptor> class NewtonSolver;

//------------------------------------------------------------------------------

template<int dim>
class MultiGridTsDesc : public TsDesc<dim> {

protected:

  int failSafeNewton;
  int maxItsNewton;
  double epsNewton;
  double epsAbsResNewton, epsAbsIncNewton;
  int lineSearch;
  int maxItsLS;
  double contractionLS, sufficDecreaseLS;

  DistSVec<bool,2> *tag;

  NewtonSolver<ImplicitTsDesc<dim> > *ns;

  IoData* myIoDataPtr;

protected:

  template<class Scalar, int neq>
  KspPrec<neq> *createPreconditioner(PcData &, Domain *);

  template<int neq>
  KspSolver<DistSVec<double,neq>, MatVecProd<dim,neq>, KspPrec<neq>, 
    Communicator> *createKrylovSolver(const DistInfo &, KspData &, MatVecProd<dim,neq> *, 
				      KspPrec<neq> *, Communicator *);

public:
  
  MultiGridTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitTsDesc();

  virtual void computeJacobian(int, DistSVec<double,dim> &, DistSVec<double,dim> &) = 0;
  virtual void setOperators(DistSVec<double,dim> &) = 0;
  virtual int solveLinearSystem(int, DistSVec<double,dim> &, DistSVec<double,dim> &) = 0;
  int solveNonLinearSystem(DistSVec<double,dim> &, const int);
  void computeFunction(int, DistSVec<double,dim> &, DistSVec<double,dim> &);  
  void recomputeFunction(DistSVec<double,dim> &, DistSVec<double,dim> &);
  void doErrorEstimation(DistSVec<double,dim> &);

  int checkFailSafe(DistSVec<double,dim>&);
  void resetFixesTag();

  int getMaxItsNewton() const { return maxItsNewton; }
  double getEpsNewton() const { return epsNewton; }
  double getEpsAbsResNewton() const { return epsAbsResNewton; }
  double getEpsAbsIncNewton() const { return epsAbsIncNewton; }
  int getLineSearch() const { return lineSearch; }
  int getMaxItsLineSearch() const { return maxItsLS; }
  double getContractionLineSearch() const { return contractionLS; }
  double getSufficientDecreaseLineSearch() const { return sufficDecreaseLS; }
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitTsDesc.C>
#endif

#endif
