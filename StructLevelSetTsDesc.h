#ifndef _STRUCT_LEVELSET_TS_DESC_H_
#define _STRUCT_LEVELSET_TS_DESC_H_

#include <TsDesc.h>

#include <IoData.h>
#include <Domain.h>
// MLX TODO REMOVE #include "Ghost/DistEulerStructGhostFluid.h"

struct DistInfo;

class GeoSource;
template<class Scalar, int dim> class DistSVec;
template<int dim> class DistExactRiemannSolver;




//------------------------------------------------------------------------

template<int dim>
class StructLevelSetTsDesc : public TsDesc<dim> {

 protected:
  DistExactRiemannSolver<dim> *riemann;
  DistLevelSetStructure *distLSS;
  DistSVec<double,dim> *Wstarij;  // stores the FS Riemann solution (i->j) along edges 
  DistSVec<double,dim> *Wstarji;  // stores the FS Riemann solution (j->i) along edges 

 public:
  StructLevelSetTsDesc(IoData &, GeoSource &, Domain *);
  ~StructLevelSetTsDesc();


  //-- overrides the functions implemented in TsDesc.
  void setupTimeStepping(DistSVec<double,dim> *, IoData &);
  double computeTimeStep(int, double *, DistSVec<double,dim> &);
  void updateStateVectors(DistSVec<double,dim> &, int = 0);
  int checkSolution(DistSVec<double,dim> &);
  void setupOutputToDisk(IoData &, bool *, int, double,
                        DistSVec<double,dim> &);
  void outputToDisk(IoData &, bool*, int, int, int, double, double,
                        DistSVec<double,dim> &);
  void outputForces(IoData &, bool*, int, int, int, double, double,
                    DistSVec<double,dim> &);
  void outputPositionVectorToDisk();
  void resetOutputToStructure(DistSVec<double,dim> &);
  void updateOutputToStructure(double, double, DistSVec<double,dim> &);
  double computeResidualNorm(DistSVec<double,dim>& );
  void monitorInitialState(int, DistSVec<double,dim>& );
  virtual int solveNonLinearSystem(DistSVec<double,dim> &)=0;

};








#ifdef TEMPLATE_FIX
#include <StructLevelSetTsDesc.C>
#endif

#endif

