#ifndef _TS_DESC_H_
#define _TS_DESC_H_

#include <IoData.h>
#include <TsInput.h>
#include <TsOutput.h>
#include <TsRestart.h>
#include <TsParameters.h>
#include <Domain.h>
#include <DistVector.h>
//#include <DistEulerStructGhostFluid.h>

class RefVal;
class VarFcn;
class GeoSource;
class DistGeoState;
class MeshMotionHandler;
class HeatTransferHandler;
class MemoryPool;
class Timer;
//class DistEulerStructGhostFluid;

template<int dim> class DistBcData;
template<int dim> class DistTimeState;
template<int dim> class SpaceOperator;
template<int dim> class PostOperator;
template<int dim> class DistExactRiemannSolver;

//------------------------------------------------------------------------------
/** Base class for time integration descriptors */
template<int dim>
class TsDesc {

public:

  typedef DistSVec<double,dim> SolVecType;
  typedef DistSVec<double,3> PosVecType;
  typedef DistVec<double> VolVecType;
  typedef DistSVec<double,1> PhiVecType;

  TsParameters *data;
  virtual double computeResidualNorm(DistSVec<double,dim>&);
  virtual void monitorInitialState(int, DistSVec<double,dim> &);

  Timer *timer;
protected:
  DistExactRiemannSolver<dim> *riemann1;

  PosVecType *X;
  VolVecType *A;
  PosVecType *Xs;

  bool* problemType;
  int numPhase;
  TsData::Clipping clippingType;
  BcsWallData::Integration wallType;
  BcsWallData::Reconstruction wallRecType;
  TsData::TimeStepCalculation timeStepCalculation;

//  TsParameters *data;
  TsInput *input;
  TsOutput<dim> *output;
  TsRestart *restart;

  DistSVec<double,dim> *V;
  DistSVec<double,dim> *R;
  DistSVec<double,dim> *Rinlet;
  DistSVec<double,dim> *Rreal;
  RefVal *refVal;
  VarFcn *varFcn;

  DistVec<int> fluidIdDummy;
  DistTimeState<dim> *timeState;
  DistBcData<dim> *bcData;
  DistGeoState *geoState;
  SpaceOperator<dim> *spaceOp;
  PostOperator<dim> *postOp;

  MeshMotionHandler* mmh;
  HeatTransferHandler* hth;
//  DistEulerStructGhostFluid* eulerFSI;

  Domain *domain;

//  Timer *timer;
  Communicator *com;

// Included (MB)
  int fixSol;
  int iForce;
  int iTotal;

  double forceNorm;
  double *forceNorms;

  bool failSafeFlag;

protected:

//  void monitorInitialState(int, DistSVec<double,dim> &);
  bool monitorConvergence(int, DistSVec<double,dim> &);

// Included (MB)
  bool monitorForceConvergence(IoData &, int, DistSVec<double,dim> &);
  bool monitorAvgForceConvergence(IoData &, int, DistSVec<double,dim> &);

public:

  TsDesc(IoData &, GeoSource &, Domain *);
  virtual ~TsDesc();

  void printf(int, const char *, ...);
  VarFcn *createVarFcn(IoData &);
  DistBcData<dim> *createBcData(IoData &);
  MeshMotionHandler *createMeshMotionHandler(IoData &, GeoSource &, MemoryPool *);

  HeatTransferHandler* createHeatTransferHandler(IoData&, GeoSource&);

  double recomputeResidual(DistSVec<double,dim> &, DistSVec<double,dim> &);
  virtual void setupTimeStepping(DistSVec<double,dim> *, IoData &);
  virtual double computeTimeStep(int, double *, DistSVec<double,dim> &);
  virtual double computePositionVector(bool *, int, double);
  void interpolatePositionVector(double, double);
  void computeMeshMetrics(int it = -1);
  virtual void updateStateVectors(DistSVec<double,dim> &, int = 0);
  bool checkForLastIteration(int, double, double, DistSVec<double,dim> &); //KW: not used?

  void setFailSafe(bool flag){ failSafeFlag = flag; }

// Modified (MB)
  bool checkForLastIteration(IoData &, int, double, double, DistSVec<double,dim> &);

  virtual void setupOutputToDisk(IoData &, bool *, int, double,
			 	DistSVec<double,dim> &);
  virtual void outputToDisk(IoData &, bool*, int, int, int, double, double,
				DistSVec<double,dim> &);
	virtual void writeStateRomToDisk(int it, double cpu) {};
	virtual void writeErrorToDisk(int it, double cpu) {};

  virtual void outputForces(IoData &, bool*, int, int, int, double, double,
		    DistSVec<double,dim> &);

  void outputPositionVectorToDisk(DistSVec<double,dim> &U);
  virtual void resetOutputToStructure(DistSVec<double,dim> &);
  virtual void updateOutputToStructure(double, double, DistSVec<double,dim> &);

  virtual bool IncreasePressure(double dt, double t, DistSVec<double,dim> &U){return true;}
  virtual int solveNonLinearSystem(DistSVec<double,dim> &U, int) { return 0; }
  virtual int checkSolution(DistSVec<double,dim> &);

  int getInitialIteration() const { return restart->iteration; }
  double getInitialTime() const { return restart->etime; }
  DistInfo &getVecInfo() const { return domain->getNodeDistInfo(); }
  DistInfo &getInletVecInfo() const {return domain->getInletNodeDistInfo(); }
  int getNumPhase() { return numPhase; }

// Included (MB)
  virtual void fixSolution(DistSVec<double,dim> &, DistSVec<double,dim> &);

  virtual void setCurrentTime(double t,DistSVec<double,dim>& U) { }

  virtual void writeBinaryVectorsToDiskRom(bool, int, double, DistSVec<double,dim> *, DistSVec<double,dim> *, VecSet<DistSVec<double,dim> > *);
  void updateGhostFluid(DistSVec<double,dim> &, Vec3D&, double);

  void printNodalDebug(int globNodeId, int identifier, DistSVec<double,dim> *U, DistVec<int> *Id=0, DistVec<int> *Id0=0);
};

//------------------------------------------------------------------------------


#ifdef TEMPLATE_FIX
#include <TsDesc.C>
#endif

#endif
