#include <TsDesc.h>
#include <Domain.h>
#include <cmath>
#include <RefVal.h>
#include <GeoSource.h>
#include <DistBcData.h>
#include <DistTimeState.h>
#include <DistGeoState.h>
#include <SpaceOperator.h>
#include <PostOperator.h>
#include <MeshMotionHandler.h>
#include <HeatTransferHandler.h>
#include <DistVector.h>
#include <MemoryPool.h>
#include <Timer.h>
#include <alloca.h>
#include <DistExactRiemannSolver.h>

extern int interruptCode;

//------------------------------------------------------------------------------

template<int dim>
TsDesc<dim>::TsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) : domain(dom),fluidIdDummy(dom->getNodeDistInfo())
{
  X = new DistSVec<double,3>(getVecInfo());
  A = new DistVec<double>(getVecInfo());
  Xs = new DistSVec<double,3>(getVecInfo());
 
  // Initialize the values
  *X = 0.0;
  *A = 0.0;
  *Xs = 0.0;

  fluidIdDummy = 0;
 
  V = new DistSVec<double,dim>(getVecInfo());
  R = new DistSVec<double,dim>(getVecInfo());
  Rinlet = new DistSVec<double,dim>(getVecInfo());
  Rreal = new DistSVec<double,dim>(getVecInfo());
  timer = domain->getTimer();
  com = domain->getCommunicator();

  problemType = ioData.problem.type;
  clippingType = ioData.ts.typeClipping;
  wallType = ioData.bc.wall.integration;
  wallRecType = ioData.bc.wall.reconstruction;
  timeStepCalculation = ioData.ts.timeStepCalculation;

  refVal = new RefVal(ioData.ref.rv);

  varFcn = new VarFcn(ioData);

  input = new TsInput(ioData);
  geoState = new DistGeoState(ioData, domain);
  // restart the geoState (positions of the mesh) At return X contains the last
  // position of the mesh.
  if(ioData.problem.framework==ProblemData::BODYFITTED || ioData.problem.framework==ProblemData::EMBEDDEDALE) 
    geoState->setup1(input->positions, X, A);
  else {
    char temp[1]; temp[0] = '\0';
    geoState->setup1(temp, X, A);
  }

  bcData = createBcData(ioData);

  spaceOp = new SpaceOperator<dim>(ioData, varFcn, bcData, geoState, domain, V);

  postOp = new PostOperator<dim>(ioData, varFcn, bcData, geoState, domain, V);

  data = new TsParameters(ioData);
  output = new TsOutput<dim>(ioData, refVal, domain, postOp);
  restart = new TsRestart(ioData, refVal);

  hth = createHeatTransferHandler(ioData, geoSource);

  riemann1 = new DistExactRiemannSolver<dim>(ioData, domain, varFcn);
// Included (MB)
  forceNorm = 0.0;
  failSafeFlag = false;
  if (ioData.sa.avgsIt) {
    forceNorms = new double[ioData.sa.avgsIt];
  }
  else {
    forceNorms = 0;
  }

  iForce = 0;
  iTotal = 0;

  modifiedGhidaglia = (ioData.schemes.bc.type==BoundarySchemeData::MODIFIED_GHIDAGLIA);

  if (ioData.sa.fixsol == 0)
    fixSol = 0;
  else if (ioData.sa.fixsol == 1)
    fixSol = 1;

	timeState = 0;
	mmh = 0; 
}

//------------------------------------------------------------------------------

template<int dim>
TsDesc<dim>::~TsDesc()
{
  if (X) delete X;
  if (Xs) delete Xs;
  if (A) delete A;
  if (V) delete V;
  if (R) delete R;
  if (Rinlet) delete Rinlet;
  if (Rreal) delete Rreal;
  if (data) delete data;
  if (input) delete input;
  if (output) delete output;
  if (restart) delete restart;
  if (refVal) delete refVal;
  if (varFcn) delete varFcn;
  if (timeState) delete timeState;
  if (bcData) delete bcData;
  if (geoState) delete geoState;
  if (spaceOp) delete spaceOp;
  if (postOp) delete postOp;
  if (mmh) delete mmh;
  if (hth) delete hth;
  if (forceNorms) delete forceNorms;
  if (riemann1) delete riemann1;
}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::printf(int verbose, const char *format, ...)
{

  if (com->cpuNum() == 0 && verbose <= com->getMaxVerbose()) {
    va_list args;
    va_start(args, format);
    vfprintf(stdout, format, args);
    ::fflush(stdout);
    va_end(args);
  }

}

//------------------------------------------------------------------------------

template<int dim>
VarFcn *TsDesc<dim>::createVarFcn(IoData &ioData)
{

  VarFcn *vf = 0;
  fprintf(stderr,"ERROR: obsolete function createVarFcn(...) called!\n");
  exit(-1);
  return vf;

}

//------------------------------------------------------------------------------

template<int dim>
DistBcData<dim> *TsDesc<dim>::createBcData(IoData &ioData)
{

  DistBcData<dim> *bc = 0;

  if (ioData.eqs.type == EquationsData::NAVIER_STOKES && 
      ioData.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
    if (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
	ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES)
      bc = new DistBcDataSA<dim>(ioData, varFcn, domain, *X);
    else if (ioData.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE)
      bc = new DistBcDataKE<dim>(ioData, varFcn, domain, *X);
  } 
  else
    bc = new DistBcDataEuler<dim>(ioData, varFcn, domain, *X);

  if (!bc) {
    com->fprintf(stderr, "*** Error: no valid choice for BCs\n");
    exit(1);
  }

  return bc;

}

//------------------------------------------------------------------------------

template<int dim>
MeshMotionHandler *TsDesc<dim>::
createMeshMotionHandler(IoData &ioData, GeoSource &geoSource, MemoryPool *mp)
{

  MeshMotionHandler *_mmh = 0;

  if (ioData.problem.type[ProblemData::AERO]) {
    if (ioData.problem.type[ProblemData::ACCELERATED])
      _mmh = new AccAeroMeshMotionHandler(ioData, varFcn, bcData->getInletPrimitiveState(),
					  geoSource.getMatchNodes(), domain, mp);
    else
      _mmh = new AeroMeshMotionHandler(ioData, varFcn, bcData->getInletPrimitiveState(),
				       geoSource.getMatchNodes(), domain, mp);
    //check that algorithm number is consistent with simulation in special case RK2-CD
    // if C0 and RK2 then RK2DGCL is needed!
    if(_mmh->getAlgNum() == 20 || _mmh->getAlgNum() == 21){
      if(ioData.ts.type == TsData::EXPLICIT &&
         (ioData.ts.expl.type == ExplicitData::RUNGE_KUTTA_2 ||
          ioData.ts.expl.type == ExplicitData::ONE_BLOCK_RK2 ||
          ioData.ts.expl.type == ExplicitData::ONE_BLOCK_RK2bis )){
        if(!(ioData.dgcl.normals    == DGCLData::EXPLICIT_RK2 &&
             ioData.dgcl.velocities == DGCLData::EXPLICIT_RK2_VEL)){
          com->fprintf(stderr, "***Error: Computation of the normals or velocities (%d,%d)\n", ioData.dgcl.normals, ioData.dgcl.velocities);
          com->fprintf(stderr, "***       is not consistent with Aeroelastic algorithm\n");
          exit(1);
        }
      }
    }
  }
  else if (ioData.problem.type[ProblemData::FORCED]) {
    if (ioData.problem.type[ProblemData::ACCELERATED])
      _mmh = new AccForcedMeshMotionHandler(ioData, varFcn, bcData->getInletPrimitiveState(), domain);
    else {
      if (ioData.forced.type == ForcedData::HEAVING)
        _mmh = new HeavingMeshMotionHandler(ioData, domain);
      else if (ioData.forced.type  == ForcedData::PITCHING)
        _mmh = new PitchingMeshMotionHandler(ioData, domain);
      else if (ioData.forced.type  == ForcedData::DEFORMING)
        _mmh = new DeformingMeshMotionHandler(ioData, domain);
    }
  }
  else if (ioData.problem.type[ProblemData::ACCELERATED])
    _mmh = new AccMeshMotionHandler(ioData, varFcn, bcData->getInletPrimitiveState(), domain);
  else if (ioData.problem.type[ProblemData::ROLL])
    _mmh = new RigidRollMeshMotionHandler(ioData, bcData->getInletAngles(), domain);
  else if (ioData.problem.type[ProblemData::RBM])
    _mmh = new RbmExtractor(ioData, domain);

  return _mmh;

}

//------------------------------------------------------------------------------

template<int dim>
HeatTransferHandler* TsDesc<dim>::createHeatTransferHandler(IoData& iod, GeoSource& gs)

{

  HeatTransferHandler* _hth = 0;

  if (iod.problem.type[ProblemData::THERMO])
    _hth = new HeatTransferHandler(iod, gs.getMatchNodes(), domain);

  return _hth;

}

//------------------------------------------------------------------------------
template<int dim>
double TsDesc<dim>::recomputeResidual(DistSVec<double,dim> &F, DistSVec<double,dim> &Finlet)
{

  return spaceOp->recomputeResidual(F,Finlet);

}

//------------------------------------------------------------------------------
template<int dim>
void TsDesc<dim>::setupTimeStepping(DistSVec<double,dim> *U, IoData &iod)
{

  geoState->setup2(timeState->getData());
  timeState->setup(input->solutions, *X, bcData->getInletBoundaryVector(), *U, iod);

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(mmh);
  DeformingMeshMotionHandler* _dmmh = dynamic_cast<DeformingMeshMotionHandler*>(mmh);
  HeavingMeshMotionHandler* _hmmh = dynamic_cast<HeavingMeshMotionHandler*>(mmh);
  PitchingMeshMotionHandler* _pmmh = dynamic_cast<PitchingMeshMotionHandler*>(mmh);

  if (_mmh)
    _mmh->setup(&restart->frequency, &data->maxTime, postOp, *X, *U);
  else if (_dmmh) 
    _dmmh->setup(*X);
  else if (_hmmh)
    _hmmh->setup(*X);
  else if (_pmmh)
    _pmmh->setup(*X);

  if (hth)
    hth->setup(&restart->frequency, &data->maxTime);

  *Xs = *X;

  initializeFarfieldCoeffs();
}

//------------------------------------------------------------------------------

template<int dim>
double TsDesc<dim>::computeTimeStep(int it, double *dtLeft, DistSVec<double,dim> &U, double angle)
{
  double t0 = timer->getTime();
  //com->fprintf(stderr,"data->residual = %lf, restart->residual = %lf.\n",data->residual, restart->residual);
  this->data->allowstop = this->timeState->allowcflstop;
  timeState->unphysical = data->unphysical;
  data->computeCflNumber(it - 1, data->residual / restart->residual, angle);
  int numSubCycles = 1;

  double dt = 0.0;
  if(failSafeFlag == false){
    if(timeStepCalculation == TsData::CFL || it==1)
      dt = timeState->computeTimeStep(data->cfl, data->dualtimecfl, dtLeft, &numSubCycles, *geoState, *X, *A, U);
    else  //time step size with error estimation
      dt = timeState->computeTimeStep(it, dtLeft, &numSubCycles);
  }
  else //if time step is repeated
    dt = this->timeState->computeTimeStepFailSafe(dtLeft, &numSubCycles);

  if(timeStepCalculation == TsData::ERRORESTIMATION && it == 1)
    this->timeState->setDtMin(dt * data->getCflMinOverCfl0());

  if (problemType[ProblemData::UNSTEADY])
    com->printf(5, "Global dt: %g (remaining subcycles = %d)\n", dt*refVal->time, numSubCycles);
  timer->addFluidSolutionTime(t0);
  timer->addTimeStepTime(t0);

  return dt;
}

//------------------------------------------------------------------------------

template<int dim>
double TsDesc<dim>::computePositionVector(bool *lastIt, int it, double t, DistSVec<double,dim> &U)
{
  double dt = 0.0;

  if (mmh && mmh->structureSubcycling()) {
    double dtleft = 0.0;
    double dtf = this->computeTimeStep(it+1, &dtleft, U); 
    mmh->storeFluidSuggestedTimestep(dtf);
  }
    
  if (mmh) {
    double t0 = timer->getTime();
    dt = mmh->updateStep1(lastIt, it, t, bcData->getVelocityVector(), *Xs, &data->maxTime);
    timer->addMeshSolutionTime(t0);
  }

  if (hth) {
    double dth = hth->updateStep1(lastIt, it, bcData->getTemperatureVector());
    if (!mmh)
      dt = dth;
  }

  if (mmh) {
    double t0 = timer->getTime();
    mmh->updateStep2(lastIt, it, t, bcData->getVelocityVector(), *Xs);
    timer->addMeshSolutionTime(t0);
  }

  if (hth) {
    hth->updateStep2(lastIt, it, bcData->getTemperatureVector());
  }

  return dt;

}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::interpolatePositionVector(double dt, double dtLeft)
{
  if (!mmh) return;

//  EmbeddedMeshMotionHandler* _mmh = dynamic_cast<EmbeddedMeshMotionHandler*>(mmh);
//  if (_mmh) return;

  geoState->interpolate(dt, dtLeft, *Xs, *X);

}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::computeMeshMetrics(int it)
{
//  EmbeddedMeshMotionHandler* _mmh = dynamic_cast<EmbeddedMeshMotionHandler*>(mmh);

  if (mmh) {
    if (it >= 0) com->fprintf(stderr, "GeoState Computing for it %d\n", it);
    double t0 = timer->getTime();
    geoState->compute(timeState->getData(), bcData->getVelocityVector(), *X, *A);
    timer->addMeshMetricsTime(t0);
    timer->addFluidSolutionTime(t0);
  }

  if (mmh || hth) 
    bcData->update(*X);

}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::updateStateVectors(DistSVec<double,dim> &U, int it)
{

  geoState->update(*X, *A);
  timeState->update(U);

  spaceOp->updateFixes(); 
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
bool TsDesc<dim>::checkForLastIteration(IoData &ioData, int it, double t, double dt, DistSVec<double,dim> &U)
{

  if ((ioData.eqs.type == EquationsData::NAVIER_STOKES) && (ioData.sa.fres)) {
    if (!problemType[ProblemData::UNSTEADY]) {
      bool forceconv = false;
      if (ioData.sa.avgsIt)
        forceconv = monitorAvgForceConvergence(ioData, it, U);
      else
        forceconv = monitorForceConvergence(ioData, it, U);
      bool fluidconv = monitorConvergence(it, U);
      if ((forceconv || fluidconv) || it >= data->maxIts) {
        if (forceconv)
          com->fprintf(stderr,"\n***** Residual of the aerodynamic force is satisfied \n\n");
        else if (fluidconv)
          com->fprintf(stderr,"\n***** Residual of the fluid solution is satisfied \n\n");
        else if (it >= data->maxIts)
          com->fprintf(stderr,"\n***** Maximum number of iteration is reached \n\n");
        return true;
      }
    }
  }
  else {
    if (!problemType[ProblemData::UNSTEADY] && monitorConvergence(it, U))
      return true;
  }

  if (!problemType[ProblemData::AERO] && !problemType[ProblemData::THERMO] && it >= data->maxIts) return true;

  if (problemType[ProblemData::UNSTEADY] )
    if(t >= data->maxTime - 0.01 * dt)
      return true;

  return false;

}

//------------------------------------------------------------------------------

template<int dim>
int TsDesc<dim>::checkSolution(DistSVec<double,dim> &U)
{

  int ierr = 0;

  if (dim == 6)
    ierr = domain->template 
      clipSolution<dim,1>(clippingType, wallType, varFcn, bcData->getInletConservativeState(), U);
  else if (dim == 7)
    ierr = domain->template 
      clipSolution<dim,2>(clippingType, wallType, varFcn, bcData->getInletConservativeState(), U);
  else
    ierr = domain->checkSolution(varFcn, U);

  return ierr;

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim> 
void TsDesc<dim>::fixSolution(DistSVec<double,dim> &U, DistSVec<double,dim> &dU)
{        
          
  if (fixSol == 1)
    domain->fixSolution(varFcn, U, dU);
    
}   

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::setupOutputToDisk(IoData &ioData, bool *lastIt, int it, double t, 
				    DistSVec<double,dim> &U)
{
  if (it == data->maxIts)
    *lastIt = true;
  else
    monitorInitialState(it, U);
  
  output->setMeshMotionHandler(ioData, mmh);
  output->openAsciiFiles();
  timer->setSetupTime();
  output->cleanProbesFile();

  if (it == 0) {
    // First time step: compute GradP before computing forces
    spaceOp->computeGradP(*X, *A, U);

    if (wallRecType==BcsWallData::CONSTANT)
      output->writeForcesToDisk(*lastIt, it, 0, 0, t, 0.0, restart->energy, *X, U);
    else //wallRecType == EXACT_RIEMANN
      output->writeForcesToDisk(*riemann1, *lastIt, it, 0, 0, t, 0.0, restart->energy, *X, U);

    output->writeLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, restart->energy, *X, U);
    output->writeHydroForcesToDisk(*lastIt, it, 0, 0, t, 0.0, restart->energy, *X, U);
    output->writeHydroLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, restart->energy, *X, U);
    output->writeResidualsToDisk(it, 0.0, 1.0, data->cfl);
    output->writeMaterialVolumesToDisk(it, 0.0, *A);
    output->writeCPUTimingToDisk(*lastIt, it, t, timer);
    output->writeBinaryVectorsToDisk(*lastIt, it, t, *X, *A, U, timeState);
    output->writeAvgVectorsToDisk(*lastIt, it, t, *X, *A, U, timeState);
    output->writeHeatFluxesToDisk(*lastIt, it, 0, 0, t, 0.0, restart->energy, *X, U);
    restart->writeKPtracesToDisk(ioData, *lastIt, it, t, *X, *A, U, timeState, domain, postOp);
    writeStateRomToDisk(it, 0.0);
    writeErrorToDisk(it, 0.0);
  }
}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::outputToDisk(IoData &ioData, bool* lastIt, int it, int itSc, int itNl, 
			       double t, double dt, DistSVec<double,dim> &U)
{

  com->globalSum(1, &interruptCode);
  if (interruptCode)
    *lastIt = true;

  double cpu = timer->getRunTime();
  double res = data->residual / restart->residual;

  output->writeLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, restart->energy, *X, U);
  output->writeHydroForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, restart->energy, *X, U);
  output->writeHydroLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, restart->energy, *X, U);
  output->writeResidualsToDisk(it, cpu, res, data->cfl);
  writeStateRomToDisk(it, cpu);
  output->writeMaterialVolumesToDisk(it, t, *A);
  output->writeCPUTimingToDisk(*lastIt, it, t, timer);
  writeErrorToDisk(it, cpu);
  output->writeBinaryVectorsToDisk(*lastIt, it, t, *X, *A, U, timeState);
  output->writeAvgVectorsToDisk(*lastIt, it, t, *X, *A, U, timeState);
  output->writeProbesToDisk(*lastIt, it, t, *X, *A, U, timeState,fluidIdDummy);
  restart->writeToDisk<dim,1>(com->cpuNum(), *lastIt, it, t, dt, *timeState, *geoState);
  output->writeHeatFluxesToDisk(*lastIt, it, itSc, itNl, t, cpu, restart->energy, *X, U);

  restart->writeKPtracesToDisk(ioData, *lastIt, it, t, *X, *A, U, timeState, domain, postOp);

  this->output->updatePrtout(t);
  if (*lastIt) {
    timer->setRunTime();
    if (com->getMaxVerbose() >= 2)
      timer->print(domain->getStrTimer());

    output->closeAsciiFiles();

  }

}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::outputForces(IoData &ioData, bool* lastIt, int it, int itSc, int itNl, 
			       double t, double dt, DistSVec<double,dim> &U)  {

  double cpu = timer->getRunTime();
  if (wallRecType==BcsWallData::CONSTANT)
    output->writeForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, restart->energy, *X, U);
  else //wallRecType==EXACT_RIEMANN
    output->writeForcesToDisk(*riemann1, *lastIt, it, itSc, itNl, t, cpu, restart->energy, *X, U);
}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::outputPositionVectorToDisk(DistSVec<double,dim> &U)
{

  *X = *Xs;

  if (mmh)  {
    int algNum = mmh->getAlgNum();
    if (algNum == 8)  return;
  }

  domain->writeVectorToFile(restart->positions[0], 0, 0.0, *Xs, &(refVal->tlength));

  if(mmh && mmh->getAlgNum() == 1)
    output->writeDisplacementVectorToDisk(1, 1.0, *X, U); 

  timer->setRunTime();
  if (com->getMaxVerbose() >= 2)
    timer->print(domain->getStrTimer());

  DistVec<double> As(getVecInfo());
  domain->computeControlVolumes(refVal->tlength, *Xs, As);

}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::resetOutputToStructure(DistSVec<double,dim> &U)
{

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(mmh);
  if (_mmh) 
    _mmh->resetOutputToStructure(postOp, *X, U);

}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::updateOutputToStructure(double dt, double dtLeft,
					  DistSVec<double,dim> &U)
{

  if (mmh) {
    double work[2];
    mmh->computeInterfaceWork(dt, postOp, geoState->getXn(), timeState->getUn(), *X, U, work);
    restart->energy[0] += work[0];
    restart->energy[1] += work[1];
  }

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(mmh);
  if (_mmh)
    _mmh->updateOutputToStructure(dt, dtLeft, postOp, *X, U);

  if (hth)
    hth->updateOutputToStructure(dt, dtLeft, postOp, *X, U);

}

//------------------------------------------------------------------------------

template<int dim>
double TsDesc<dim>::computeResidualNorm(DistSVec<double,dim>& U)
{
  if (wallRecType==BcsWallData::CONSTANT)
    spaceOp->computeResidual(*X, *A, U, *R, timeState);
  else //wallRecTyp == ExactRiemann
    spaceOp->computeResidual(riemann1, *X, *A, U, *R, timeState);
    
  spaceOp->applyBCsToResidual(U, *R);

  double res = 0.0;
  if (data->resType == -1){
    res = (*R)*(*R);

 
  }else{ 
    int iSub;
    const DistInfo& distInfo = R->info();
#ifndef MPI_OMP_REDUCTION
    double* allres = reinterpret_cast<double*>(alloca(sizeof(double) * distInfo.numGlobSub));
    for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) 
      allres[iSub] = 0.0;
#endif

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
    for (iSub=0; iSub<distInfo.numLocSub; ++iSub) {
      double (*r)[dim] = R->subData(iSub);
      bool* flag = R->getMasterFlag(iSub);
      double locres = 0.0;
      for (int i=0; i<R->subSize(iSub); ++i) {
	if (flag[i])
	  locres += r[i][data->resType]*r[i][data->resType];
      }
#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif
    }

#ifdef MPI_OMP_REDUCTION
    distInfo.com->globalSum(1, &res);
#else
    distInfo.com->globalSum(distInfo.numGlobSub, allres);
    res = 0.0;
    for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) 
      res += allres[iSub];
#endif
  }

  return sqrt(res);

}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::monitorInitialState(int it, DistSVec<double,dim> &U)
{

  //com->printf(2, "State vector norm = %.12e\n", sqrt(U*U));

  if (!problemType[ProblemData::UNSTEADY]) {
    double trhs = timer->getTimeSyncro();
    data->residual = computeResidualNorm(U);
    trhs = timer->getTimeSyncro() - trhs;
    if (it == 0)
      restart->residual = data->residual;
    if (data->resType == -1)
      com->printf(2, "Spatial residual norm = %.12e\n", data->residual);
    else
      com->printf(2, "Spatial residual norm[%d] = %.12e\n", data->resType, data->residual);
    com->printf(2, "Time for one residual evaluation: %f s\n", trhs);
  }

  com->printf(2, "\n");

}

//------------------------------------------------------------------------------

template<int dim>
bool TsDesc<dim>::monitorConvergence(int it, DistSVec<double,dim> &U)
{
  data->residual = computeResidualNorm(U);
  if ((problemType[ProblemData::AERO] || problemType[ProblemData::THERMO]) && (it == 1 || it == 2))
    restart->residual = data->residual;

  if (data->residual == 0.0 || data->residual < data->eps * restart->residual) 
    return true;
  else
    return false;

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
bool TsDesc<dim>::monitorForceConvergence(IoData &ioData, int it, DistSVec<double,dim> &U)
{

  double forceNorm_n;
  double resForce;

  int nSurfs = postOp->getNumSurf();
    
  Vec3D x0, F, M;

  Vec3D *Fi = new Vec3D[nSurfs];
  Vec3D *Mi = new Vec3D[nSurfs];
  Vec3D *Fv = new Vec3D[nSurfs];
  Vec3D *Mv = new Vec3D[nSurfs];

  x0[0] = ioData.output.transient.x0;
  x0[1] = ioData.output.transient.y0;
  x0[2] = ioData.output.transient.z0;

  postOp->computeForceAndMoment(x0, *this->X, U, 0, Fi, Mi, Fv, Mv);

  F = 0.0;
  M = 0.0;
  
  F = Fi[0] + Fv[0];
  M = Mi[0] + Mv[0];

  forceNorm_n = sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2]);

  resForce = fabs(forceNorm_n - forceNorm) / forceNorm_n;

  forceNorm = forceNorm_n;

  com->fprintf(stderr,"\n***** It = %d, Force residual = %e, Target = %e\n\n", it, resForce, ioData.sa.fres);

  if ((resForce == 0.0) || (resForce < ioData.sa.fres))
    return true;
  else
    return false;

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
bool TsDesc<dim>::monitorAvgForceConvergence(IoData &ioData, int it, DistSVec<double,dim> &U)
{

  double avgForceNorm;
  double forceNorm_n;
  double resForce;

  if (forceNorms == 0)
    forceNorms = new double[ioData.sa.avgsIt];

  int nSurfs = postOp->getNumSurf();

  Vec3D x0, F, M;

  Vec3D *Fi = new Vec3D[nSurfs];
  Vec3D *Mi = new Vec3D[nSurfs];
  Vec3D *Fv = new Vec3D[nSurfs];
  Vec3D *Mv = new Vec3D[nSurfs];

  x0[0] = ioData.output.transient.x0;
  x0[1] = ioData.output.transient.y0;
  x0[2] = ioData.output.transient.z0;

  postOp->computeForceAndMoment(x0, *this->X, U, 0, Fi, Mi, Fv, Mv);

  F = 0.0;
  M = 0.0;

  F = Fi[0] + Fv[0];
  M = Mi[0] + Mv[0];

  forceNorm_n = sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2]);

  if (iTotal < ioData.sa.avgsIt) {
    forceNorms[iForce] = forceNorm_n;
    ++iForce;
    avgForceNorm = 0.0;
    for (int i = 0; i < iForce; ++i)
      avgForceNorm += forceNorms[i];
    avgForceNorm *= 1.0/double(iForce);
  }
  else {
    if (iForce == ioData.sa.avgsIt)
      iForce = 0;
    forceNorms[iForce] = forceNorm_n;
    ++iForce;
    avgForceNorm = 0.0;
    for (int i = 0; i < ioData.sa.avgsIt; ++i)
      avgForceNorm += forceNorms[i];
    avgForceNorm *= 1.0/double(ioData.sa.avgsIt);
  }

  if (iTotal)
    resForce = fabs(forceNorm_n - avgForceNorm) / avgForceNorm;
  else
    resForce = 1.0;

  ++iTotal;

  com->fprintf(stderr,"\n***** It = %d, Force residual = %e, Target = %e\n\n", it, resForce, ioData.sa.fres);

  if ((resForce == 0.0) || (resForce < ioData.sa.fres))
    return true;
  else
    return false;

}

//-----------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::updateGhostFluid(DistSVec<double,dim> &U, Vec3D& totalForce, double dt)
{
/*
  if (eulerFSI)
    eulerFSI->updateGhostFluid(X, U, totalForce, dt);
*/
}

//-----------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::printNodalDebug(int globNodeId, int identifier, DistSVec<double,dim> *U, DistVec<int> *Id, DistVec<int> *Id0)
{ //Kevin:For debug only!
  int nSub = domain->getNumLocSub();
  SubDomain **sub = domain->getSubDomain();
  for(int iSub=0; iSub<nSub; iSub++) {
    int* locToGlob = sub[iSub]->getNodeMap();
    for(int i=0; i<(*U)(iSub).size(); i++)
      if(locToGlob[i]+1==globNodeId) { 
        fprintf(stderr,"*** %d Node %d: U = %e %e %e %e %e ", identifier, globNodeId, 
                (*U)(iSub)[i][0], (*U)(iSub)[i][1], (*U)(iSub)[i][2], (*U)(iSub)[i][3], (*U)(iSub)[i][4]);
        if(Id)
          fprintf(stderr,", Id = %d ", (*Id)(iSub)[i]);
        if(Id0)
          fprintf(stderr,", Id0 = %d ", (*Id0)(iSub)[i]);
        fprintf(stderr,"\n");
      }
  }
}

//----------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::writeBinaryVectorsToDiskRom(bool lastIt, int it, double t,
		DistSVec<double,dim> *F1 = NULL, DistSVec<double,dim> *F2 = NULL, VecSet< DistSVec<double,dim> > *F3 = NULL)

{

  output->writeBinaryVectorsToDiskRom(lastIt, it, t, F1, F2, F3);

}

//----------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::computeDistanceToWall(IoData &ioData)
{
  // Nothing to do here by default.
}

//----------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::updateFarfieldCoeffs(double dt)
{
  if(!modifiedGhidaglia) return;
  int nSub = domain->getNumLocSub();
  SubDomain **sub = domain->getSubDomain();
  int iSub;
#pragma omp parallel for
  for(iSub=0; iSub<nSub; iSub++)
    sub[iSub]->updateFarfieldCoeffs(dt);
}

//----------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::updateBoundaryExternalState()
{
  if(!modifiedGhidaglia) return;
  int nSub = domain->getNumLocSub();
  SubDomain **sub = domain->getSubDomain();
  int iSub;
#pragma omp parallel for
  for(iSub=0; iSub<nSub; iSub++)
    sub[iSub]->updateBoundaryExternalState();
}

//----------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::initializeFarfieldCoeffs()
{
  if(!modifiedGhidaglia) return;
  double *Vin = bcData->getInletPrimitiveState();
  double soundspeed = varFcn->computeSoundSpeed(Vin);
  double gamma;
  if (varFcn->getType() == VarFcnBase::STIFFENEDGAS ||
      varFcn->getType() == VarFcnBase::PERFECTGAS)
    gamma = varFcn->getGamma();
  else
    gamma = varFcn->getBetaWater();

  double HH_init = -2.0*soundspeed/(gamma - 1.0);
  //fprintf(stderr,"HH_init is set to %e.\n", HH_init);

  int nSub = domain->getNumLocSub();
  SubDomain **sub = domain->getSubDomain();
  int iSub;
#pragma omp parallel for
  for(iSub=0; iSub<nSub; iSub++)
    sub[iSub]->initializeFarfieldCoeffs(HH_init);
}

