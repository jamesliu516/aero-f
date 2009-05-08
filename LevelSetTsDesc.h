#ifndef _LEVELSET_TS_DESC_H_
#define _LEVELSET_TS_DESC_H_

#include <TsDesc.h>

#include <IoData.h>
#include <Domain.h>
#include <LevelSet.h>

struct DistInfo;

class GeoSource;
template<class Scalar, int dim> class DistSVec;
template<int dim> class DistExactRiemannSolver;

                                                                                                                

//------------------------------------------------------------------------

template<int dim>
class LevelSetTsDesc : public TsDesc<dim> {

 protected:
  LevelSet *LS;
  DistExactRiemannSolver<dim> *riemann;
  DistVec<double> Phi;           //conservative variables
  DistVec<double> PhiV;          //primitive variables
  DistSVec<double,dim> Vg;       //primitive V for GFMP
  DistSVec<double,dim> *Vgf;     //primitive V storage for phase change (if extrapolation)
  DistVec<double> *Vgfweight;

  // multiphase conservation check
  DistSVec<double,dim> boundaryFlux;
  DistSVec<double,dim> interfaceFlux;
  DistSVec<double,dim> computedQty;
  DistSVec<double,dim> *tmpDistSVec;
  DistSVec<double,dim> *tmpDistSVec2;
  double expectedTot[dim];
  double expectedF1[dim];
  double expectedF2[dim];
  double computedTot[dim];
  double computedF1[dim];
  double computedF2[dim];

  // frequency for reinitialization of level set
  int frequencyLS;

  MultiFluidData::InterfaceType interfaceType; //to advance levelset or not
  
  //buckling cylinder parameters
  // pressure is increased in the fluid at rate Prate from
  // initial pressure Pinit until it reaches the pressure
  // given by boundary conditions which happens at tmax.
  double tmax;
  double Prate;
  double Pinit;

 public:
  LevelSetTsDesc(IoData &, GeoSource &, Domain *);
  ~LevelSetTsDesc();

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
  void resetOutputToStructure(DistSVec<double,dim> &);
  void updateOutputToStructure(double, double, DistSVec<double,dim> &);

  void conservationErrors(DistSVec<double,dim> &U, int it);


  bool IncreasePressure(double dt, double t, DistSVec<double,dim> &U);
  virtual int solveNonLinearSystem(DistSVec<double,dim> &)=0;

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <LevelSetTsDesc.C>
#endif

#endif
