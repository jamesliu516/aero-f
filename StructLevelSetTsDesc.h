#ifndef _STRUCT_LEVELSET_TS_DESC_H_
#define _STRUCT_LEVELSET_TS_DESC_H_

#include <TsDesc.h>

#include <IoData.h>
#include <Domain.h>
#include <LevelSet.h>
#include <PostOperator.h>
// MLX TODO REMOVE #include "Ghost/DistEulerStructGhostFluid.h"

struct DistInfo;
class DynamicNodalTransfer;
class GeoSource;
template<class Scalar, int dim> class DistSVec;
template<int dim> class DistExactRiemannSolver;




//------------------------------------------------------------------------

template<int dim>
class StructLevelSetTsDesc : public TsDesc<dim> , ForceGenerator<dim> {

 protected:
  DistExactRiemannSolver<dim> *riemann; //Riemann solver -- used at both FF and FS interfaces

  DistVec<int> nodeTag; // = 1 for fluid #1; = -1 for fluid #2.
  DistVec<int> nodeTag0; // node tag for the previous time-step.

  double timeStep;

  double (*Fs)[3]; //force distribution on the structure surfac3
  int numStructNodes;
  int numStructElems;

  // coefficients for piston simulations.
  Vec3D fsiPosition;
  Vec3D fsiNormal;
  double fsiVelocity;
  double pressureRef;

  // ----------- time steps -----------------------------------------------------------
  double dtf;     //<! fluid time-step
  double dtfLeft; //<! time until next structure time-step is reached.
  double dts;     //<! structure time-step
  // ----------------------------------------------------------------------------------

  // ----------- components for Fluid-Structure interface. -----------------------------
  DistLevelSetStructure *distLSS; //<! tool for FS tracking (not necessarily a  "levelset solver".)
  DistSVec<double,dim> *Wstarij;  //<! stores the FS Riemann solution (i->j) along edges
  DistSVec<double,dim> *Wstarji;  //<! stores the FS Riemann solution (j->i) along edges
  DistSVec<double,dim> Vtemp;     //<! the primitive variables.
  DistSVec<double,dim> *VWeights; //<! stores U*Weights for each node. Used in updating phase change.
  DistVec<double> *Weights;       //<! weights for each node. Used in updating phase change.
  // ------------------------------------------------------------------------------------

  // ----------- components for Fluid-Fluid interface -----------------------------------
  LevelSet *LS;
  DistVec<double> Phi;            //<! conservative variables
  DistVec<double> PhiV;           //<! primitive variables
  DistSVec<double,dim> Vg;        //<! primitive V for GFMP
  DistSVec<double,dim> *Vgf;      //<! primitive V storage for phase change (if extrapolation)
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

  MultiFluidData::InterfaceType interfaceTypeFF; //to advance levelset or not
  // --------------------------------------------------------------------------------------

  DynamicNodalTransfer *dynNodalTransfer;

 public:
  int orderOfAccuracy; // consistent with the reconstruction type for space
  int forceApp; // now have four options.
                // = 1 : on GammaF, formula 1;
                // = 2 : on GammaF, formula 2;
                // = 3 : on Gamma*, formula 1;
                // = 4 : on Gamma*, formula 2;
  const int TYPE;   //TYPE = 1 (for fluid-fullbody)
                    //     = 2 (for fluid-shell-fluid)

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
  /** Override the TsDesc routine because forces are sent to the structure
   * in a different way than the general case */
  void updateOutputToStructure(double, double, DistSVec<double,dim> &);
//  double computePositionVector(bool *lastIt, int it, double t);

  double computeResidualNorm(DistSVec<double,dim>& );
  void monitorInitialState(int, DistSVec<double,dim>& );
  void conservationErrors(DistSVec<double,dim> &U, int it);

  void updateFSInterface();
  void updateNodeTag();

  void computeForceLoad(); //compute Fs.

  virtual int solveNonLinearSystem(DistSVec<double,dim> &)=0;

  void getForcesAndMoments(DistSVec<double,dim> &U, DistSVec<double,3> &X,
                                           double F[3], double M[3]);

};








#ifdef TEMPLATE_FIX
#include <StructLevelSetTsDesc.C>
#endif

#endif

