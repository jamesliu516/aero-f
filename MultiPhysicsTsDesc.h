#ifndef _MULTIPHYSICS_TS_DESC_H_
#define _MULTIPHYSICS_TS_DESC_H_

#include <GhostPoint.h>
#include <TsDesc.h>

struct DistInfo;
class Domain;
class DynamicNodalTransfer;
class GeoSource;
class IoData;

template<int dim> class DistExactRiemannSolver;
template<class Scalar, int dim> class DistSVec;
template<int dimLS> class LevelSet;

//------------------------------------------------------------------------

template<int dim, int dimLS>
class MultiPhysicsTsDesc : public TsDesc<dim> , ForceGenerator<dim> {

 protected:

  double vfar[dim]; //farfield state
  DistExactRiemannSolver<dim> *riemann; //Riemann solver -- used at both FF and FS interfaces
  enum Type {EULER = 0, NAVIER_STOKES = 1} eqsType;
  const int numFluid;

  // time-step and related.
  double dtf;     //<! fluid time-step
  double dtfLeft; //<! time until next structure time-step is reached.
  double dts;     //<! structure time-step
  int globIt;         //<! current global(i.e. structure) iteration
  bool inSubCycling;  //<! is it in subcyling (i.e. itSc>1)

  bool requireSpecialBDF;

  //------------------------------------------------------------------------
  // EulerFSI: basic parameters
  int orderOfAccuracy; // consistent with the reconstruction type for space
  bool linRecAtInterface;
  int simType;        // 0: steady-state    1: unsteady
  int riemannNormal;  // 0: struct normal;  1: fluid normal (w.r.t. control volume face)
                      // 2: averaged structure normal;
  int numStructNodes;
  int numStructElems;
  int forceApp; // now have four options.
                // = 1 : on GammaF, formula 1;
                // = 2 : on GammaF, formula 2; (not used)
                // = 3 : on Gamma*, formula 1;
                // = 4 : on Gamma*, formula 2; (not used)
  int phaseChangeChoice; // = 0. use nodal values.
                         // = 1. use solutions of Riemann problems.

  // EulerFSI: force calculation
  double (*Fs)[3]; //force distribution on the structure surfac3
  bool FsComputed; //whether Fs has been computed for this (fluid-)time step.

  bool existsWstarnm1;

  // EulerFSI: interface tracking
  DistLevelSetStructure *distLSS; //<! tool for FS tracking (not necessarily a  "levelset solver".)

  DistSVec<double,dim> *Wstarij,*Wstarij_nm1;  //<! stores the FS Riemann solution (i->j) along edges
  DistSVec<double,dim> *Wstarji,*Wstarji_nm1;  //<! stores the FS Riemann solution (j->i) along edges
  DistSVec<double,dim> Vtemp;     //<! the primitive variables.
  DistSVec<double,dim> *VWeights; //<! stores U*Weights for each node. Used in updating phase change.
  DistVec<double> *Weights;       //<! weights for each node. Used in updating phase change.
  DistVec<GhostPoint<dim>*> *ghostPoints;

  DistSVec<double,3> *Nsbar;      //<! cell-averaged structure normal (optional)

  // EulerFSI: FS communication
  DistSVec<double,dim> Wtemp;
  DynamicNodalTransfer *dynNodalTransfer;

  //------------------------------------------------------------------------
  // MultiPhaseFlow: basics
  MultiPhaseSpaceOperator<dim,dimLS> *multiPhaseSpaceOp;
  int frequencyLS; // frequency for reinitialization of level set

  // MultiPhaseFlow: level-sets and fluidIds
  LevelSet<dimLS> *LS;
  FluidSelector fluidSelector;
  DistSVec<double,dimLS> Phi;           //conservative variables
  DistSVec<double,dimLS> PhiWeights;    //<! stores Phi*Weights for each ndoe. Used in updating phase change
  DistSVec<double,dimLS> PhiV;          //primitive variables
  DistSVec<double,dim> V0;

  //------------------------------------------------------------------------
  // buckling cylinder parameters
  // pressure is increased in the fluid at rate Prate from
  // initial pressure Pinit until it reaches the pressure
  // given by boundary conditions which happens at tmax.
  double tmax;
  double Prate;
  double Pinit;
  double Pscale;

 protected:
  void setupEmbeddedFSISolver(IoData &ioData);
  void setupMultiPhaseFlowSolver(IoData &ioData);
  /** computes the force load. Wij and Wji must be edge-based primitive state vectors. */ 
  void computeForceLoad(DistSVec<double,dim> *Wij, DistSVec<double,dim> *Wji);

 public:
  MultiPhysicsTsDesc(IoData &, GeoSource &, Domain *);
  ~MultiPhysicsTsDesc();

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

  void getForcesAndMoments(DistSVec<double,dim> &U, DistSVec<double,3> &X,
                                           double F[3], double M[3]);

  bool IncreasePressure(double dt, double t, DistSVec<double,dim> &U);

  virtual int solveNonLinearSystem(DistSVec<double,dim> &)=0;

};

#ifdef TEMPLATE_FIX
#include <MultiPhysicsTsDesc.C>
#endif

#endif
