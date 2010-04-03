#ifndef _LEVEL_SET_H_
#define _LEVEL_SET_H_

#include "DistVector.h"

class IoData;
class Domain;
class TimeData;
class Communicator;

#ifndef _DNDGRAD_TMPL_
#define _DNDGRAD_TMPL_
template<int dim, class Scalar = double> class DistNodalGrad;
#endif

// This LevelSet class stores the level-set at the previous time-steps
// and helps compute the next-step level-set as well as the
// reinitialization of the level-set
// Note that the current level-set is not stored here, but is available
// in LevelSetTsDesc.

template<int dimLS>
class LevelSet {

 private:
  Communicator *com;
  int numLocSub;
  Domain *domain;
  TimeData *data;

  bool trueLevelSet[dimLS]; // some level-sets take only values -1 and +1 (volumeID initialization)
                            // and do not need to be advected, since they don't represent a
                            // fluid-fluid interface. However, other level-sets are distance functions
                            // to a fluid-fluid interface and must therefore be updated and reinitialized.

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
  DistSVec<double,dimLS> Phi0;
  DistSVec<double,1> Psi;		// the steady state solution of Psi will reinitialize the level set
  DistVec<double> dt;			// pseudo-time steps
  DistSVec<double,dimLS> PsiRes;        // residual of the reinitialization equation
  DistVec<double> w;			// weights in epxlicit positive coefficient scheme (cf Barth and Sethian)
  DistNodalGrad<dimLS, double> *lsgrad;	// least-square gradient to compute normal at each node
  DistVec<int> Tag;			// node tags for reinitialization in a narrow band

 public:

  DistSVec<double,dimLS> Phin;
  DistSVec<double,dimLS> Phinm1;
  DistSVec<double,dimLS> Phinm2;

  LevelSet(IoData &iod, Domain *domain);
  ~LevelSet();

  // initialization routines
  template<int dim>
  void setup(const char *name, DistSVec<double,3> &X, DistSVec<double,dim> &U,
             DistSVec<double,dimLS> &Phi, IoData &iod);
  void setupPhiVolumesInitialConditions(IoData &iod, DistSVec<double,dimLS> &Phi);
  void setupPhiOneDimensionalSolution(IoData &iod, DistSVec<double,3> &X,
                                      DistSVec<double,dimLS> &Phi);
  void setupPhiMultiFluidInitialConditions(IoData &iod, 
                              DistSVec<double,3> &X, DistSVec<double,dimLS> &Phi);

  // update the members of the class and writing them to the disk
  void checkTrueLevelSetUpdate(DistSVec<double,dimLS> &dPhi);
  void update(DistSVec<double,dimLS> &Phi);
  void writeToDisk(char *name);

  // reinitialization routines
  template<int dim>
  void reinitializeLevelSet(DistGeoState &geoState,
                            DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                            DistSVec<double,dim> &U, DistSVec<double,dimLS> &Phi);
  template<int dim>
  void reinitializeLevelSetFM(DistGeoState &geoState,
                            DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                            DistSVec<double,dim> &U, DistSVec<double,dimLS> &Phi);

  // switching from conservative (rho*phi) to primitive (phi) and vice-versa
  template<int dim>
  void conservativeToPrimitive(DistSVec<double,dimLS> &Cons, DistSVec<double,dimLS> &Prim,
                               DistSVec<double,dim> &U);

  template<int dim>
  void primitiveToConservative(DistSVec<double,dimLS> &Prim, DistSVec<double,dimLS> &Cons,
	                       DistSVec<double,dim> &U);
};

#ifdef TEMPLATE_FIX
#include <LevelSet.C>
#endif

#endif


