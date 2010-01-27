#ifndef _TIMER_H_
#define _TIMER_H_

#include <IoData.h>

class Communicator;
class IoData;
//------------------------------------------------------------------------------

class Timer {

  enum TimerIndex {
    setup, run, total,
    fluid, nodalWeights, nodalGrad, fvTerm, feTerm, 
    fvJac, feJac, vms, dvms, h2Assembly, fluidPrecSetup, fluidKsp, meshMetrics, structUpd, 
    mesh, meshAssembly, meshPrecSetup, meshKsp,
    podConstr, snapsLinSolv, padeReconstr, correlMatrix, eigSolv, gramSchmidt,
    romSol, romConstr, romTimeInteg,
    comm, localCom, globalCom, interCom,
    io, binread, binwrite,
    levelSet, lsNodalWeightsAndGrad, lsFvTerm, lsKsp,
    timeStep, intersect, embedPhaseChange, eulerFSI
  };

  int numTimings;

  double initialTime;

  int *counter;
  double *data;
  
  IoData *ioData;
  
  Communicator *com;

public:

  Timer(Communicator *);
  ~Timer();

  double getTime();
  double getTimeSyncro();
  double getRunTime();

  void setIoData(IoData &_ioData);
  void setSetupTime();
  void setRunTime();

  double addTimeStepTime(double);
  double addNodalWeightsTime(double);
  double addNodalGradTime(double);
  double addFiniteVolumeTermTime(double);
  double addFiniteElementTermTime(double);
  double addFiniteVolumeJacTime(double);
  double addFiniteElementJacTime(double);
  double addVMSLESTime(double);
  double addDynamicVMSLESTime(double);
  double addH2SetupTime(double);
  double addPrecSetupTime(double);
  double addKspTime(double);
  double addMeshMetricsTime(double);
  double addStructUpdTime(double);
  double addMeshSolutionTime(double);
  double addMeshAssemblyTime(double);
  double addMeshPrecSetupTime(double);
  double addMeshKspTime(double);
  double removeForceAndDispComm(double);
  double addPodConstrTime(double);
  double addSnapsLinSolvTime(double);
  double addPadeReconstrTime(double);
  double addCorrelMatrixTime(double);
  double addEigSolvTime(double);
  double addGramSchmidtTime(double);
  double addRomSolTime(double);
  double addRomConstrTime(double);
  double addRomTimeIntegTime(double);
  double addLocalComTime(double);
  double addGlobalComTime(double);
  double addInterComTime(double);
  double addBinaryReadTime(double);
  double addBinaryWriteTime(double);
  double addFluidSolutionTime(double);

  // Level-Set Timer Functions
  double addLevelSetSolutionTime(double);
  double addLSNodalWeightsAndGradTime(double);
  double addLSFiniteVolumeTermTime(double);
  double addLSKspTime(double);

  // Embedded FSI Timer Functions
  double addIntersectionTime(double);
  double addEmbedPhaseChangeTime(double);
//  double addEmbedComTime(double);
  double removeIntersAndPhaseChange(double);

  void print(Timer *, FILE * = stdout);

};

//------------------------------------------------------------------------------

#endif
