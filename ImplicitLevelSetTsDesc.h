#ifndef _IMPLICIT_LEVELSET_TS_DESC_H_
#define _IMPLICIT_LEVELSET_TS_DESC_H_

#include <ImplicitCoupledTsDesc.h>

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
class ImplicitLevelSetTsDesc : public ImplicitCoupledTsDesc<dim> {


  LevelSet *LS;
  MatVecProd<dim,1> *mvpLS;
  KspSolver<DistVec<double>, MatVecProd<dim,1>, KspPrec<dim>, Communicator> *kspLS;
  DistSVec<double,dim> Vg;

 public:
 
  ImplicitLevelSetTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitLevelSetTsDesc();

  
  LevelSet *getLevelSet() {return LS;};
  
  //-- overrides the functions implemented in TsDesc.
  virtual void setupTimeStepping(DistSVec<double,dim> *, IoData &);
  virtual double computeTimeStep(int, double *, DistSVec<double,dim> &);
  virtual void updateStateVectors(DistSVec<double,dim> &);
  virtual int checkSolution(DistSVec<double,dim> &);
  virtual void setupOutputToDisk(IoData &, bool *, int, double, 
                        DistSVec<double,dim> &);
  virtual void outputToDisk(IoData &, bool*, int, int, int, double, double, 
		    	DistSVec<double,dim> &);
  virtual void resetOutputToStructure(DistSVec<double,dim> &);
  virtual void updateOutputToStructure(double, double, DistSVec<double,dim> &);

  //-- overrides the functions implemented in ImplicitTsDesc.
  virtual int solveNonLinearSystem(DistSVec<double,dim> &);
  virtual void computeFunction(int, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void computeFunctionLS(int, DistVec<double> &, DistVec<double> &,  DistVec<double> &,
                         DistSVec<double,dim> &,DistVec<double> &);
  virtual void recomputeFunction(DistSVec<double,dim> &, DistSVec<double,dim> &);
  
  virtual int checkFailSafe(DistSVec<double,dim>&);
  virtual void resetFixesTag();

  //-- overrides the functions implemented in ImplicitCoupledTsDesc.
  virtual void computeJacobian(int, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void computeJacobianLS(int, DistVec<double> &, DistVec<double> &,
                         DistVec<double> &, DistSVec<double,dim> &,DistVec<double> &);
  int solveLinearSystemLS(int, DistVec<double> &, DistVec<double> &);



};


//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitLevelSetTsDesc.C>
#endif

#endif
