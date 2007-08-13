#ifndef _EXPLICIT_LEVELSET_TS_DESC_H_
#define _EXPLICIT_LEVELSET_TS_DESC_H_

#include <ExplicitTsDesc.h>

#include <IoData.h>
#include <KspPrec.h>
#include <Domain.h>
#include <TsOutput.h>
#include <TsInput.h>
#include <TsRestart.h>
#include <TsParameters.h>

struct DistInfo;

class GeoSource;
class DistGeoState;

class LevelSet; 
template<class Scalar, int dim> class DistSVec;
template<int dim, int neq> class MatVecProd;

                                                                                                                
#ifndef _KSPSLVR_TMPL_
#define _KSPSLVR_TMPL_
template<class VecType, class MvpOp, class PrecOp, class IoOp, class ScalarT = double> class KspSolver;
#endif


//------------------------------------------------------------------------

template<int dim>
class ExplicitLevelSetTsDesc : public ExplicitTsDesc<dim> {


  LevelSet *LS;
  DistSVec<double,dim> *Vgf;
  DistSVec<double,dim> Vg;
  NewtonSolver<ExplicitTsDesc<dim> > *ns;
  MatVecProd<dim,1> *mvpLS;
  KspPrec<dim> *pc;
  KspSolver<DistVec<double>, MatVecProd<dim,1>, KspPrec<dim>, Communicator> *kspLS;

 public:
 
  ExplicitLevelSetTsDesc(IoData &, GeoSource &, Domain *);
  ~ExplicitLevelSetTsDesc();


  template<int neq>
  KspSolver<DistVec<double>, MatVecProd<dim,1>, KspPrec<neq,double>,
  Communicator, double> *createKrylovSolverLS(const DistInfo &, KspData &, MatVecProd<dim,1> *,
                              KspPrec<neq,double> *, Communicator *);
  
  LevelSet *getLevelSet() {return LS;};
  
  //-- overrides the functions implemented in TsDesc.
  virtual void setupTimeStepping(DistSVec<double,dim> *, IoData &);
  virtual double computeTimeStep(int, double *, DistSVec<double,dim> &);
  virtual void updateStateVectors(DistSVec<double,dim> &);
  virtual int checkSolution(DistSVec<double,dim> &);
  virtual void setupOutputToDisk(IoData &, bool *, int, double, DistSVec<double,dim> &);
  virtual void outputToDisk(IoData &, bool*, int, int, int, double, double, 
		    DistSVec<double,dim> &);
  virtual void resetOutputToStructure(DistSVec<double,dim> &);
  virtual void updateOutputToStructure(double, double, DistSVec<double,dim> &);

  //-- overrides the functions implemented in ExplicitTsDesc.
  virtual int solveNonLinearSystem(DistSVec<double,dim> &);
  virtual void computeRKUpdate(DistSVec<double,dim>&, 
			DistSVec<double,dim>&, DistSVec<double,dim>&);

  virtual void computeFunctionLS(int, DistVec<double> &, DistVec<double> &,  DistVec<double> &,
                         DistSVec<double,dim> &,DistVec<double> &);

  //-- overrides the functions implemented in ExplicitCoupledTsDesc.
  virtual void computeJacobianLS(int, DistVec<double> &, DistVec<double> &,
                         DistVec<double> &, DistSVec<double,dim> &,DistVec<double> &);
  virtual int solveLinearSystemLS(int, DistVec<double> &, DistVec<double> &);



};


//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ExplicitLevelSetTsDesc.C>
#endif

#endif
