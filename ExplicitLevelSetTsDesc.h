#ifndef _EXPLICIT_LEVELSET_TS_DESC_H_
#define _EXPLICIT_LEVELSET_TS_DESC_H_

#include <ExplicitTsDesc.h>

#include <IoData.h>
#include <KspPrec.h>
#include <Domain.h>

struct DistInfo;

class GeoSource;
class LevelSet; 
template<class Scalar, int dim> class DistSVec;
//template<int dim, int neq> class MatVecProd;

//#ifndef _KSPSLVR_TMPL_
//#define _KSPSLVR_TMPL_
//template<class VecType, class MvpOp, class PrecOp, class IoOp, class ScalarT = double> class KspSolver;
//#endif

//------------------------------------------------------------------------

template<int dim>
class ExplicitLevelSetTsDesc : public LevelSetTsDesc<dim> {

 private:
  //bool RK4;
  ExplicitData::Type timeIntegrationType;

  DistSVec<double,dim> U0;
  DistSVec<double,dim> k1;
  DistSVec<double,dim> k2;
  DistSVec<double,dim> k3;
  DistSVec<double,dim> k4;

  DistVec<double> Phi0;
  DistVec<double> p1;
  DistVec<double> p2;
  DistVec<double> p3;
  DistVec<double> p4;

 public:
  ExplicitLevelSetTsDesc(IoData &, GeoSource &, Domain *);
  ~ExplicitLevelSetTsDesc();

  int solveNonLinearSystem(DistSVec<double,dim> &);

 private:
  void solveNonLinearSystemEuler(DistSVec<double,dim> &);
  void solveNonLinearSystemLevelSet(DistSVec<double,dim> &);
  void computeRKUpdate(DistSVec<double,dim>&, DistSVec<double,dim>&, int);
  void computeRKUpdateLS(DistVec<double>&, DistVec<double>&, DistSVec<double,dim>&);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ExplicitLevelSetTsDesc.C>
#endif

#endif
