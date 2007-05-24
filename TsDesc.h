#ifndef _TS_DESC_H_
#define _TS_DESC_H_

#include <IoData.h>
#include <TsInput.h>
#include <TsOutput.h>
#include <TsRestart.h>
#include <TsParameters.h>
#include <Domain.h>
#include <DistVector.h>

class RefVal;
class VarFcn;
class GeoSource;
class DistGeoState;
class MeshMotionHandler;
class HeatTransferHandler;
class MemoryPool;
class Timer;

template<int dim> class DistBcData;
template<int dim> class DistTimeState;
template<int dim> class SpaceOperator;
template<int dim> class PostOperator;

//------------------------------------------------------------------------------

template<int dim>
class TsDesc {

public:

  typedef DistSVec<double,dim> SolVecType;
  typedef DistSVec<double,3> PosVecType;
  typedef DistVec<double> VolVecType;

protected:

  PosVecType *X;
  VolVecType *A;
  PosVecType *Xs;

  bool* problemType;
  int numPhase;
  TsData::Clipping clippingType;
  BcsWallData::Integration wallType;

  TsParameters *data;
  TsInput *input;
  TsOutput<dim> *output;
  TsRestart *restart;

  DistSVec<double,dim> *V;
  DistSVec<double,dim> *R;
  DistSVec<double,dim> *Rinlet;

  RefVal *refVal;
  VarFcn *varFcn;

  DistTimeState<dim> *timeState;
  DistBcData<dim> *bcData;
  DistGeoState *geoState;
  SpaceOperator<dim> *spaceOp;
  PostOperator<dim> *postOp;

  MeshMotionHandler* mmh;
  HeatTransferHandler* hth;

  Domain *domain;

  Timer *timer;
  Communicator *com;

protected:

  double computeResidualNorm(DistSVec<double,dim>&);
  void monitorInitialState(int, DistSVec<double,dim> &);
  bool monitorConvergence(int, DistSVec<double,dim> &);

public:

  TsDesc(IoData &, GeoSource &, Domain *);
  ~TsDesc();

  void printf(int, const char *, ...);
  VarFcn *createVarFcn(IoData &);
  DistBcData<dim> *createBcData(IoData &);
  MeshMotionHandler *createMeshMotionHandler(IoData &, GeoSource &, MemoryPool *);
  HeatTransferHandler* createHeatTransferHandler(IoData&, GeoSource&);

  double recomputeResidual(DistSVec<double,dim> &, DistSVec<double,dim> &);
  virtual void setupTimeStepping(DistSVec<double,dim> *, IoData &);
  virtual double computeTimeStep(int, double *, DistSVec<double,dim> &);
  double computePositionVector(bool *, int, double);
  void interpolatePositionVector(double, double);
  void computeMeshMetrics();
  virtual void updateStateVectors(DistSVec<double,dim> &, int);
  bool checkForLastIteration(int, double, double, DistSVec<double,dim> &);

  virtual void setupOutputToDisk(IoData &, bool *, int, double, 
			 	DistSVec<double,dim> &);
  virtual void outputToDisk(IoData &, bool*, int, int, int, double, double, 
				DistSVec<double,dim> &);

  virtual void outputForces(IoData &, bool*, int, int, int, double, double, 
		    DistSVec<double,dim> &);

  void outputPositionVectorToDisk();
  void resetOutputToStructure(DistSVec<double,dim> &);
  void updateOutputToStructure(double, double, DistSVec<double,dim> &);

  virtual int solveNonLinearSystem(DistSVec<double,dim> &U) { return 0; }
  virtual int checkSolution(DistSVec<double,dim> &);

  int getInitialIteration() const { return restart->iteration; }
  double getInitialTime() const { return restart->etime; }
  DistInfo &getVecInfo() const { return domain->getNodeDistInfo(); }
  DistInfo &getInletVecInfo() const {return domain->getInletNodeDistInfo(); }

};

//------------------------------------------------------------------------------


#ifdef TEMPLATE_FIX
#include <TsDesc.C>
#endif

#endif
