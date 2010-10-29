#ifndef _LEVELSET_TS_DESC_H_
#define _LEVELSET_TS_DESC_H_

#include <TsDesc.h>

class IoData;
class GeoSource;
class Domain;
class FluidSelector;
template<int dimLS> class LevelSet;
template<int dim> class DistExactRiemannSolver;

struct DistInfo;
template<class Scalar, int dim> class DistSVec;

//------------------------------------------------------------------------

template<int dim, int dimLS>
class LevelSetTsDesc : public TsDesc<dim> {

 protected:
  MultiPhaseSpaceOperator<dim,dimLS> *multiPhaseSpaceOp;
  FluidSelector fluidSelector;
  LevelSet<dimLS> *LS;
  DistExactRiemannSolver<dim> *riemann;
  DistSVec<double,dimLS> Phi;           //conservative variables
  DistSVec<double,dimLS> PhiV;          //primitive variables
  DistSVec<double,dim> V0;

  DistVec<double> umax;
  // multiphase conservation check
  DistSVec<double,dim> boundaryFlux;
  DistSVec<double,dim> interfaceFlux;
  DistSVec<double,dim> computedQty;
  DistSVec<double,dim> *tmpDistSVec;
  DistSVec<double,dim> *tmpDistSVec2;
  double **expected;//[dimLS+2][dim]; // dimLS+1 different fluids and the total
  double **computed;//[dimLS+2][dim];

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
  
  void fixSolution(DistSVec<double,dim>& U, DistSVec<double,dim>& dU);

  virtual int solveNonLinearSystem(DistSVec<double,dim> &)=0;

 protected:
  void avoidNewPhaseCreation(DistSVec<double,dimLS> &localPhi);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <LevelSetTsDesc.C>
#endif

#endif
