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
  MultiFluidData::CopyCloseNodes copy;             // true if nodes close to interface are unchanged
  int bandlevel;                                   // number of node layers to reinitialize 
  int localtime;                                   // time stepping is local or global
  int subIt;                                       // max number of subiterations
  double cfl_psi;                                  // cfl for fictitious time step
  double conv_eps;
  bool diff;

  double invertGasLiquid; //to run Liquid in Gas Simulation....

  // for reinitialization
  DistSVec<double,1> Phi0;
  DistSVec<double,1> Psi;		// the steady state solution of Psi will reinitialize the level set
  DistVec<double> dt;			// pseudo-time steps
  DistSVec<double,1> PsiRes;		// residual of the reinitialization equation
  DistVec<double> w;			// weights in epxlicit positive coefficient scheme (cf Barth and Sethian)
  DistNodalGrad<1, double> *lsgrad;	// least-square gradient to compute normal at each node
  DistVec<int> Tag;			// node tags for reinitialization in a narrow band

 public:

  DistSVec<double,1> Phin;
  DistSVec<double,1> Phinm1;
  DistSVec<double,1> Phinm2;

  LevelSet(IoData &iod, Domain *domain);
  ~LevelSet();

  template<int dim>
  void setup(const char *name, DistSVec<double,3> &X, DistSVec<double,dim> &U,
             DistSVec<double,1> &Phi, IoData &iod);
  void setupPhiVolumesInitialConditions(IoData &iod, DistSVec<double,1> &Phi);
  void setupPhiMultiFluidInitialConditions(IoData &iod, 
                              DistSVec<double,3> &X, DistSVec<double,1> &Phi);
  void update(DistSVec<double,1> &Phi);
  void writeToDisk(char *name);

  template<int dim>
  void reinitializeLevelSet(DistGeoState &geoState,
                            DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                            DistSVec<double,dim> &U, DistSVec<double,1> &Phi);
  template<int dim>
  void reinitializeLevelSetPDE(DistGeoState &geoState,
                            DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                            DistSVec<double,dim> &U, DistSVec<double,1> &Phi);
  template<int dim>
  void computeSteadyState(DistGeoState &geoState,
                          DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                          DistSVec<double,dim> &U, DistSVec<double,1> &Phi);
  template<int dim>
  void reinitializeLevelSetFM(DistGeoState &geoState,
                            DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                            DistSVec<double,dim> &U, DistSVec<double,1> &Phi);

  bool checkConvergencePsi(int iteration, double &res0);


  template<int dim>
  void conservativeToPrimitive(DistSVec<double,1> &Cons, DistSVec<double,1> &Prim,
                               DistSVec<double,dim> &U);

  template<int dim>
  void primitiveToConservative(DistSVec<double,1> &Prim, DistSVec<double,1> &Cons,
	                       DistSVec<double,dim> &U);
};

#ifdef TEMPLATE_FIX
#include <LevelSet.C>
#endif

#endif


