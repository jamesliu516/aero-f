#ifndef _LEVEL_SET_H_
#define _LEVEL_SET_H_

#include <IoData.h>
#include <Domain.h>
#include <TimeData.h>


class LevelSet {

 private:
  Communicator *com;
  int numLocSub;
  Domain *domain;
  TimeData *data;

  // for reinitialization from Input file
  MultiFluidData::InterfaceTracking typeTracking;  // how to find interface location
  int bandlevel;                                   // number of node layers to reinitialize 
  int localtime;                                   // time stepping is local or global
  int subIt;                                       // max number of subiterations
  double cfl_psi;                                  // cfl for fictitious time step

  double invertGasLiquid; //to run Liquid in Gas Simulation....

  // for reinitialization
  DistVec<double> Phi0;
  DistSVec<double,1> Psi;		// the steady state solution of Psi will reinitialize the level set
  DistVec<double> dt;			// pseudo-time steps
  DistSVec<double,1> PsiRes;		// residual of the reinitialization equation
  DistVec<double> w;			// weights in epxlicit positive coefficient scheme (cf Barth and Sethian)
  DistNodalGrad<1, double> *lsgrad;	// least-square gradient to compute normal at each node
  DistVec<int> Tag;			// node tags for reinitialization in a narrow band

 public:

  DistVec<double> Phin;
  DistVec<double> Phinm1;
  DistVec<double> Phinm2;

  LevelSet(IoData &iod, Domain *domain);
  ~LevelSet();

  template<int dim>
  void setup(char *name, DistSVec<double,3> &X, DistVec<double> &Phi,
             DistSVec<double,dim> &U, IoData &iod);
  void update(DistVec<double> &Phi);
  void writeToDisk(char *name);

  template<int dim>
  void reinitializeLevelSet(DistGeoState &geoState,
                            DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                            DistSVec<double,dim> &U, DistVec<double> &Phi);
  template<int dim>
  void computeSteadyState(DistGeoState &geoState,
                          DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                          DistSVec<double,dim> &U, DistVec<double> &Phi);

  bool checkConvergencePsi(int iteration, double &res0);


  template<int dim>
  void conservativeToPrimitive(DistVec<double> &Cons, DistVec<double> &Prim,
                               DistSVec<double,dim> &U);

  template<int dim>
  void primitiveToConservative(DistVec<double> &Prim, DistVec<double> &Cons,
	                       DistSVec<double,dim> &U);
};

#ifdef TEMPLATE_FIX
#include <LevelSet.C>
#endif

#endif


