#ifndef _EMBEDDED_TS_DESC_H_
#define _EMBEDDED_TS_DESC_H_

#include <TsDesc.h>

#include <IoData.h>
#include <PostOperator.h>
#include <GhostPoint.h>

struct DistInfo;
class DynamicNodalTransfer;
class GeoSource;
class Domain;
template<class Scalar, int dim> class DistSVec;
template<int dim> class DistExactRiemannSolver;

//------------------------------------------------------------------------

template<int dim>
class EmbeddedTsDesc : public TsDesc<dim> , ForceGenerator<dim> {

 protected:
  DistExactRiemannSolver<dim> *riemann; //Riemann solver -- used at both FF and FS interfaces
  double vfar[dim]; //farfield state

  DistVec<int> nodeTag; // = 1 for fluid #1; = -1 for fluid #2.
  DistVec<int> nodeTag0; // node tag for the previous time-step.

  double timeStep;

  bool withCracking;

  double (*Fs)[3]; //force distribution on the structure surfac3
  bool FsComputed; //whether Fs has been computed for this (fluid-)time step.
  int numStructNodes;
  int totStructNodes;
  bool linRecAtInterface;
  int simType;        // 0: steady-state    1: unsteady
  int riemannNormal;  // 0: struct normal;  1: fluid normal (w.r.t. control volume face)
                      // 2: averaged structure normal;

  bool increasingPressure;
 
  // ----------- time steps -----------------------------------------------------------
  double dtf;     //<! fluid time-step
  double dtfLeft; //<! time until next structure time-step is reached.
  double dts;     //<! structure time-step
  int globIt;         //<! current global(i.e. structure) iteration
  bool inSubCycling;  //<! is it in subcyling (i.e. itSc>1)
  // ----------------------------------------------------------------------------------

  bool existsWstarnm1;

  // ----------- components for Fluid-Structure interface. -----------------------------
  DistLevelSetStructure *distLSS; //<! tool for FS tracking (not necessarily a  "levelset solver".)
  DistSVec<double,dim> *Wstarij,*Wstarij_nm1;  //<! stores the FS Riemann solution (i->j) along edges
  DistSVec<double,dim> *Wstarji,*Wstarji_nm1;  //<! stores the FS Riemann solution (j->i) along edges
  DistSVec<double,dim> Vtemp;     //<! the primitive variables.
  DistSVec<double,dim> *VWeights; //<! stores U*Weights for each node. Used in updating phase change.
  DistVec<double> *Weights;       //<! weights for each node. Used in updating phase change.

  DistSVec<double,3> *Nsbar;      //<! cell-averaged structure normal (optional)
  // ------------------------------------------------------------------------------------

  // Copies for fail safe ----- -----------------------------
  DistSVec<double,dim> *WstarijCopy,*Wstarij_nm1Copy;  //<! stores the FS Riemann solution (i->j) along edges
  DistSVec<double,dim> *WstarjiCopy,*Wstarji_nm1Copy;  //<! stores the FS Riemann solution (j->i) along edges
  DistVec<int> *nodeTagCopy; // = 1 for fluid #1; = -1 for fluid #2.
  DistVec<int> *nodeTag0Copy; // node tag for the previous time-step.
  DistSVec<double,dim> *UCopy;     //<! the primitive variables.
  // ------------------------------------------------------------------------------------

  DistSVec<double,dim> Wtemp;
  DynamicNodalTransfer *dynNodalTransfer;

  //buckling cylinder parameters
  // pressure is increased in the fluid at rate Prate from
  // initial pressure Pinit until it reaches the pressure
  // given by boundary conditions which happens at tmax.
  double tmax;
  double Prate;
  double Pinit;
  double Pscale;

  // Adam 04/06/2010
  enum Type {EULER = 0, NAVIER_STOKES = 1} eqsType;

  DistVec<GhostPoint<dim>*> *ghostPoints;
  // ghostPoints is a pointer on a vector of pointer on GhostPoint<dim> objects. Each subdomain can be accessed separatly.

 public:
  int orderOfAccuracy; // consistent with the reconstruction type for space
  int forceApp; // now have four options.
                // = 1 : on GammaF, formula 1;
                // = 2 : on GammaF, formula 2;
                // = 3 : on Gamma*, formula 1;
                // = 4 : on Gamma*, formula 2;
  int phaseChangeChoice; // = 0. use nodal values.
                         // = 1. use solutions of Riemann problems.
  const int numFluid;   //numFluid = 1 (for fluid-fullbody)
                            //     = 2 (for fluid-shell-fluid)

  EmbeddedTsDesc(IoData &, GeoSource &, Domain *);
  ~EmbeddedTsDesc();


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
  void outputPositionVectorToDisk(DistSVec<double,dim>&);
  void resetOutputToStructure(DistSVec<double,dim> &);
  /** Override the TsDesc routine because forces are sent to the structure
   * in a different way than the general case */
  void updateOutputToStructure(double, double, DistSVec<double,dim> &);
//  double computePositionVector(bool *lastIt, int it, double t);

  double computeResidualNorm(DistSVec<double,dim>& );
  void monitorInitialState(int, DistSVec<double,dim>& );

  void computeForceLoad(DistSVec<double,dim> *Wij, DistSVec<double,dim> *Wji);
  /** computes the force load. Wij and Wji must be edge-based primitive state vectors. */ 

  virtual int solveNonLinearSystem(DistSVec<double,dim> &, int)=0;

  void getForcesAndMoments(DistSVec<double,dim> &U, DistSVec<double,3> &X,
                                           double F[3], double M[3]);

  bool IncreasePressure(double dt, double t, DistSVec<double,dim> &U);

  void fixSolution(DistSVec<double,dim>& U,DistSVec<double,dim>& dU);
};


#ifdef TEMPLATE_FIX
#include <EmbeddedTsDesc.C>
#endif

#endif

