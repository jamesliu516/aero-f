#include <stdio.h>
#include <stdlib.h>
#include <alloca.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <Timer.h>
#include <Communicator.h>

//------------------------------------------------------------------------------

Timer::Timer(Communicator *communicator) : com(communicator)
{

  ioData = 0;
  initialTime = getTime();

  numTimings = 46;
  
  counter = new int[numTimings];
  data = new double[numTimings];

  for (int i=0; i<numTimings; ++i) {
    counter[i] = 0;
    data[i] = 0.0;
  }

}

//------------------------------------------------------------------------------

Timer::~Timer()
{

  if (counter) delete [] counter;
  if (data) delete [] data;

}

//------------------------------------------------------------------------------

double Timer::getTime()
{

  static double micro = 1.e-6;

  timeval tp;
  struct timezone tz;
  
  gettimeofday(&tp, &tz);

  // return 1000.0*tp.tv_sec + tp.tv_usec/1000.0;
  return double(tp.tv_sec) + double(tp.tv_usec) * micro;

}

//------------------------------------------------------------------------------

double Timer::getTimeSyncro()
{

  if (com) com->barrier();

  return getTime();

}

//------------------------------------------------------------------------------

double Timer::getRunTime()
{

  return getTime() - data[setup] - initialTime;

}

//------------------------------------------------------------------------------
/*
double Timer::getCpuTime()
{

  static struct rusage r;
  static double micro = (1./1.e6);

  double t;

  getrusage(RUSAGE_SELF,&r);

  // maximum resident set size utilized (in kilobytes)
  // size = r.ru_maxrss;

  // amount of time spent executing in user mode (in seconds)
  t = double(r.ru_utime.tv_sec) + double(r.ru_utime.tv_usec) * micro;

  // amount of time spent in the system executing on behalf of the process(es)
  t += double(r.ru_stime.tv_sec) + double(r.ru_stime.tv_usec) * micro;

  return t;

}
*/

//------------------------------------------------------------------------------
void Timer::setIoData(IoData &_ioData)
{

ioData = &_ioData;

/*if (ioData->problem.alltype == ProblemData::_POD_CONSTRUCTION_)
  numTimings += 4;
else if (ioData->problem.alltype == ProblemData::_ROM_AEROELASTIC_)
  numTimings += 3;  
*/
}

//------------------------------------------------------------------------------

void Timer::setSetupTime() 
{ 

  counter[setup]++;
  data[setup] = getTime() - initialTime; 

}

//------------------------------------------------------------------------------

void Timer::setRunTime() 
{ 

  counter[run]++;
  data[run] = getRunTime(); 

}

//------------------------------------------------------------------------------

double Timer::addTimeStepTime(double t0) 
{ 

  double t = getTime() - t0;
  
  counter[timeStep]++;
  data[timeStep] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addNodalWeightsTime(double t0) 
{ 

  double t = getTime() - t0;
  
  counter[nodalWeights]++;
  data[nodalWeights] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addNodalGradTime(double t0) 
{ 
  
  double t = getTime() - t0;

  counter[nodalGrad]++;
  data[nodalGrad] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addFiniteVolumeTermTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[fvTerm]++;
  data[fvTerm] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addFiniteElementTermTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[feTerm]++;
  data[feTerm] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addFiniteVolumeJacTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[fvJac]++;
  data[fvJac] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addFiniteElementJacTime(double t0) 
{ 
  
  double t = getTime() - t0;

  counter[feJac]++;
  data[feJac] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addVMSLESTime(double t0)
{

  double t = getTime() - t0;

  counter[vms]++;
  data[vms] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addDynamicVMSLESTime(double t0)
{
  double t = getTime() - t0;
  counter[dvms]++;
  data[dvms] += t;
  return t;
}

//------------------------------------------------------------------------------
	   
double Timer::addH2SetupTime(double t0) 
{ 
  
  double t = getTime() - t0;

  counter[h2Assembly]++;
  data[h2Assembly] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addPrecSetupTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[fluidPrecSetup]++;
  data[fluidPrecSetup] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addKspTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[fluidKsp]++;
  data[fluidKsp] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addMeshMetricsTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[meshMetrics]++;
  data[meshMetrics] += t; 

  return t;

}

//------------------------------------------------------------------------------
                                                                                                                             
double Timer::addStructUpdTime(double t0)
{

  double t = getTime() - t0;

  counter[structUpd]++;
  data[structUpd] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addFluidSolutionTime(double t0)
{

  double t = getTime() - t0;

  counter[fluid]++;
  data[fluid] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addMeshSolutionTime(double t0)
{

  double t = getTime() - t0;

  counter[mesh]++;
  data[mesh] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addMeshAssemblyTime(double t0)
{

  double t = getTime() - t0;

  counter[meshAssembly]++;
  data[meshAssembly] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addMeshPrecSetupTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[meshPrecSetup]++;
  data[meshPrecSetup] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addMeshKspTime(double t0)
{

  double t = getTime() - t0;

  counter[meshKsp]++;
  data[meshKsp] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::removeForceAndDispComm(double t0)
{

  double t = getTime() - t0;

  data[mesh] -= t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addPodConstrTime(double t0)
{

  double t = getTime() - t0;

  counter[podConstr]++;
  data[podConstr] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addSnapsLinSolvTime(double t0)
{

  double t = getTime() - t0;

  counter[snapsLinSolv]++;
  data[snapsLinSolv] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addPadeReconstrTime(double t0)
{

  double t = getTime() - t0;

  counter[padeReconstr]++;
  data[padeReconstr] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addCorrelMatrixTime(double t0)
{

  double t = getTime() - t0;

  counter[correlMatrix]++;
  data[correlMatrix] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addEigSolvTime(double t0)
{

  double t = getTime() - t0;

  counter[eigSolv]++;
  data[eigSolv] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addGramSchmidtTime(double t0)
{

  double t = getTime() - t0;

  counter[gramSchmidt]++;
  data[gramSchmidt] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addRomSolTime(double t0)
{

  double t = getTime() - t0;

  counter[romSol]++;
  data[romSol] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addRomConstrTime(double t0)
{

  double t = getTime() - t0;

  counter[romConstr]++;
  data[romConstr] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addRomTimeIntegTime(double t0)
{

  double t = getTime() - t0;

  counter[romTimeInteg]++;
  data[romTimeInteg] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addLocalComTime(double t0) 
{ 

  double t = getTime() - t0;

  //counter[localCom]++;
  data[localCom] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addGlobalComTime(double t0) 
{ 

  double t = getTime() - t0;

  //counter[globalCom]++;
  data[globalCom] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addInterComTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[interCom]++;
  data[interCom] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addBinaryReadTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[binread]++;
  data[binread] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addBinaryWriteTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[binwrite]++;
  data[binwrite] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addLevelSetSolutionTime(double t0)
{

  double t = getTime() - t0;

  counter[levelSet]++;
  data[levelSet] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addLSNodalWeightsAndGradTime(double t0)  {

  double t = getTime() - t0;

  counter[lsNodalWeightsAndGrad]++;
  data[lsNodalWeightsAndGrad] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addLSFiniteVolumeTermTime(double t0)
{

  double t = getTime() - t0;

  counter[lsFvTerm]++;
  data[lsFvTerm] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addLSKspTime(double t0)
{

  double t = getTime() - t0;

  counter[lsKsp]++;
  data[lsKsp] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addWaitAndReceiveDisp(double t0)
{

  double t = getTime() - t0;

  counter[waitrec]++;
  data[waitrec] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addIntersectionTime(double t0)
{
  double t = getTime() - t0;
  
  counter[intersect]++;
  data[intersect] += t;
  data[eulerFSI] += t;

  return t;
}

//------------------------------------------------------------------------------

double Timer::addEmbedPhaseChangeTime(double t0)
{
  double t = getTime() - t0;

  counter[embedPhaseChange]++;
  data[embedPhaseChange] += t;
  data[eulerFSI] += t;

  return t;
}

//------------------------------------------------------------------------------
/*
double Timer::addEmbedComTime(double t0)
{
  double t = getTime() - t0;

  counter[embedCom]++;
  data[embedCom] += t;

  return t;
}
*/
//------------------------------------------------------------------------------

double Timer::removeIntersAndPhaseChange(double t0)
{

  double t = getTime() - t0;

  data[fluid] -= t;

  return t;

}

//------------------------------------------------------------------------------
// note: the timings of both fluid and mesh parts contain their communication
void Timer::print(Timer *str, FILE *fp)
{

  if (!com) return;

  double *tmin = reinterpret_cast<double *>(alloca(numTimings * sizeof(double)));
  double *tmax = reinterpret_cast<double *>(alloca(numTimings * sizeof(double)));
  double *tavg = reinterpret_cast<double *>(alloca(numTimings * sizeof(double)));

  if (str) {
    counter[interCom] += str->counter[interCom];
    data[interCom] += str->data[interCom];
  }

  data[total] = data[setup] + data[run];


  data[comm] = data[localCom] + data[globalCom] + data[interCom];
  data[io] = data[binread] + data[binwrite];
  
  if (ioData->problem.alltype == ProblemData::_POD_CONSTRUCTION_)
    data[podConstr] -= data[io];

  int i;

  for (i=0; i<numTimings ; ++i) {
    tmin[i] = data[i];
    tmax[i] = data[i];
    tavg[i] = data[i];
  }

  com->globalMin(numTimings, tmin);
  com->globalMax(numTimings, tmax);
  com->globalSum(numTimings, tavg);

  int numCPU = com->size();

  for (i=0; i<numTimings ; ++i)
    tavg[i] /= numCPU;

  com->fprintf(fp, "\n");
  com->fprintf(fp, "----------------------------------------------------------------------\n");
  com->fprintf(fp, "Elapsed Time Report (s)       :        Min        Max        Avg   # Calls\n");
  com->fprintf(fp, "\n");
  com->fprintf(fp, "Problem Setup                 : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[setup], tmax[setup], tavg[setup], 
	       counter[setup]);
  com->fprintf(fp, "\n");
  com->fprintf(fp, "Fluid Solution                : %10.2f %10.2f %10.2f         -\n", 
	       tmin[fluid], tmax[fluid], tavg[fluid]);
  com->fprintf(fp, "  Time Steps                  : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[timeStep], tmax[timeStep], tavg[timeStep], 
	       counter[timeStep]);
  com->fprintf(fp, "  Nodal Weights and Gradients : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[nodalGrad], tmax[nodalGrad], tavg[nodalGrad], 
	       counter[nodalGrad]);
  com->fprintf(fp, "  FV Fluxes                   : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[fvTerm], tmax[fvTerm], tavg[fvTerm], 
	       counter[fvTerm]);
  com->fprintf(fp, "  FE Fluxes                   : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[feTerm], tmax[feTerm], tavg[feTerm], 
	       counter[feTerm]);
  com->fprintf(fp, "  FV Jacobian                 : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[fvJac], tmax[fvJac], tavg[fvJac], 
	       counter[fvJac]);
  com->fprintf(fp, "  FE Jacobian                 : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[feJac], tmax[feJac], tavg[feJac], 
	       counter[feJac]);
  com->fprintf(fp, "  VMS-LES Modeling            : %10.2f %10.2f %10.2f %9d\n",
               tmin[vms], tmax[vms], tavg[vms],
               counter[vms]);
  com->fprintf(fp, "  Dynamic VMS-LES Modeling    : %10.2f %10.2f %10.2f %9d\n",
               tmin[dvms], tmax[dvms], tavg[dvms],
               counter[dvms]);
  com->fprintf(fp, "  H2 Matrix Assembly          : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[h2Assembly], tmax[h2Assembly], tavg[h2Assembly], 
	       counter[h2Assembly]);
  com->fprintf(fp, "  Preconditioner Setup        : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[fluidPrecSetup], tmax[fluidPrecSetup], tavg[fluidPrecSetup], 
	       counter[fluidPrecSetup]);
  com->fprintf(fp, "  Linear Solver               : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[fluidKsp], tmax[fluidKsp], tavg[fluidKsp], 
	       counter[fluidKsp]);
  com->fprintf(fp, "  Mesh Metrics Update         : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[meshMetrics], tmax[meshMetrics], tavg[meshMetrics], 
	       counter[meshMetrics]);
  if (ioData->problem.alltype == ProblemData::_UNSTEADY_LINEARIZED_AEROELASTIC_)  {
    com->fprintf(fp, "  Structural Update           : %10.2f %10.2f %10.2f %9d\n",
               tmin[structUpd], tmax[structUpd], tavg[structUpd],
               counter[structUpd]);
  }
  com->fprintf(fp, "\n");

  // Output Mesh solution time (except for Euler FSI)
  if(ioData->strucIntersect.intersectorName == 0) {
    com->fprintf(fp, "Mesh Solution                 : %10.2f %10.2f %10.2f         -\n", 
                 tmin[mesh], tmax[mesh], tavg[mesh]);
    com->fprintf(fp, "  K Matrix Assembly           : %10.2f %10.2f %10.2f %9d\n", 
                 tmin[meshAssembly], tmax[meshAssembly], tavg[meshAssembly], 
                 counter[meshAssembly]);
    com->fprintf(fp, "  Preconditioner Setup        : %10.2f %10.2f %10.2f %9d\n", 
                 tmin[meshPrecSetup], tmax[meshPrecSetup], tavg[meshPrecSetup], 
	         counter[meshPrecSetup]);
    com->fprintf(fp, "  Linear Solver               : %10.2f %10.2f %10.2f %9d\n", 
                 tmin[meshKsp], tmax[meshKsp], tavg[meshKsp], 
	         counter[meshKsp]);
    com->fprintf(fp, "\n");
  }

  // Output Level-Set Timers
  if (ioData->eqs.numPhase > 1) {
    com->fprintf(fp, "LevelSet Solution             : %10.2f %10.2f %10.2f         -\n", tmin[levelSet], tmax[levelSet], tavg[levelSet]);

    com->fprintf(fp, "  Nodal Weights and Grad      : %10.2f %10.2f %10.2f %9d\n",
               tmin[lsNodalWeightsAndGrad], tmax[lsNodalWeightsAndGrad], tavg[lsNodalWeightsAndGrad],
               counter[lsNodalWeightsAndGrad]);
    com->fprintf(fp, "  FV Fluxes                   : %10.2f %10.2f %10.2f %9d\n",
               tmin[lsFvTerm], tmax[lsFvTerm], tavg[lsFvTerm], counter[lsFvTerm]);
    com->fprintf(fp, "  Linear Solver               : %10.2f %10.2f %10.2f %9d\n", tmin[lsKsp], tmax[lsKsp], tavg[lsKsp], counter[lsKsp]);
    com->fprintf(fp, "\n");
  }


  // Output POD Timers
  if (ioData->problem.alltype == ProblemData::_POD_CONSTRUCTION_) {
    com->fprintf(fp, "POD Basis Construction        : %10.2f %10.2f %10.2f         -\n",
               tmin[podConstr], tmax[podConstr], tavg[podConstr]);
    com->fprintf(fp, "  Snapshot Linear Solver      : %10.2f %10.2f %10.2f %9d\n",
               tmin[snapsLinSolv], tmax[snapsLinSolv], tavg[snapsLinSolv], counter[snapsLinSolv]);
    if (ioData->linearizedData.padeReconst == LinearizedData::TRUE) {
      com->fprintf(fp, "  Pade Reconstruction       : %10.2f %10.2f %10.2f %9d\n",
               tmin[padeReconstr], tmax[padeReconstr], tavg[padeReconstr], counter[padeReconstr]);
    }
    com->fprintf(fp, "  Correlation Matrix          : %10.2f %10.2f %10.2f %9d\n",
               tmin[correlMatrix], tmax[correlMatrix], tavg[correlMatrix], counter[correlMatrix]);
    com->fprintf(fp, "  SVD Solver                : %10.2f %10.2f %10.2f %9d\n",
               tmin[eigSolv], tmax[eigSolv], tavg[eigSolv], counter[eigSolv]);
    com->fprintf(fp, "  Gram-Schmidt                : %10.2f %10.2f %10.2f %9d\n",
               tmin[gramSchmidt], tmax[gramSchmidt], tavg[gramSchmidt], counter[gramSchmidt]);
    com->fprintf(fp, "\n");
  }
  else if (ioData->problem.alltype == ProblemData::_ROM_AEROELASTIC_)  {
    com->fprintf(fp, "ROM Solution                  : %10.2f %10.2f %10.2f         -\n",
              tmin[romSol], tmax[romSol], tavg[romSol]);
    com->fprintf(fp, "  ROM Construction            : %10.2f %10.2f %10.2f %9d\n",
               tmin[romConstr], tmax[romConstr], tavg[romConstr], counter[romConstr]);
    com->fprintf(fp, "  ROM Time Integration        : %10.2f %10.2f %10.2f %9d\n",
               tmin[romTimeInteg], tmax[romTimeInteg], tavg[romTimeInteg], counter[romTimeInteg]);
    com->fprintf(fp, "\n");
  }

  if (ioData->strucIntersect.intersectorName != 0) {
    com->fprintf(fp, "Eulerian FSI                  : %10.2f %10.2f %10.2f         -\n",
                 tmin[eulerFSI], tmax[eulerFSI], tavg[eulerFSI]);
    com->fprintf(fp, "  F-S Intersections           : %10.2f %10.2f %10.2f %9d\n", 
  	         tmin[intersect], tmax[intersect], tavg[intersect], 
  	         counter[intersect]);
    com->fprintf(fp,"\n");
  }

  com->fprintf(fp, "Communication/Synchronization : %10.2f %10.2f %10.2f         -\n", 
	       tmin[comm], tmax[comm], tavg[comm]);
  com->fprintf(fp, "  Local                       : %10.2f %10.2f %10.2f         -\n",
	       tmin[localCom], tmax[localCom], tavg[localCom]);
  com->fprintf(fp, "  Global                      : %10.2f %10.2f %10.2f         -\n", 
	       tmin[globalCom], tmax[globalCom], tavg[globalCom]);
  com->fprintf(fp, "  Inter                       : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[interCom], tmax[interCom], tavg[interCom], 
	       counter[interCom]);
  com->fprintf(fp, "\n");
  com->fprintf(fp, "I/O                           : %10.2f %10.2f %10.2f         -\n", 
	       tmin[io], tmax[io], tavg[io]);
  com->fprintf(fp, "  Binary Read                 : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[binread], tmax[binread], tavg[binread], 
	       counter[binread]);
  com->fprintf(fp, "  Binary Write                : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[binwrite], tmax[binwrite], tavg[binwrite], 
	       counter[binwrite]);

  com->fprintf(fp, "\n");
  com->fprintf(fp, "Total Simulation              : %10.2f %10.2f %10.2f         -\n", 
	       tmin[total], tmax[total], tavg[total]);
  com->fprintf(fp, "----------------------------------------------------------------------\n");

}

//------------------------------------------------------------------------------
