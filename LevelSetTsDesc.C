#include <LevelSetTsDesc.h>

#include "IoData.h"
#include "GeoSource.h"
#include "Domain.h"
#include "LevelSet.h"
#include <DistExactRiemannSolver.h>
#include <FluidSelector.h>

#include <cmath>
                                                                                                        
#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
#endif


#ifdef TYPE_MAT
#define MatScalar TYPE_MAT
#else
#define MatScalar double
#endif

#ifdef TYPE_PREC
#define PrecScalar TYPE_PREC
#else
#define PrecScalar double
#endif

//------------------------------------------------------------------------------

template<int dim, int dimLS>
LevelSetTsDesc<dim,dimLS>::
LevelSetTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  TsDesc<dim>(ioData, geoSource, dom), Phi(this->getVecInfo()), V0(this->getVecInfo()),
  PhiV(this->getVecInfo()),
  fluidSelector(ioData.eqs.numPhase, ioData, dom),umax(this->getVecInfo()), programmedBurn(NULL),Utilde(this->getVecInfo())

{
  multiPhaseSpaceOp = new MultiPhaseSpaceOperator<dim,dimLS>(ioData, this->varFcn, this->bcData, this->geoState, 
                                                             this->domain, this->V);
  this->timeState = new DistTimeState<dim>(ioData, multiPhaseSpaceOp, this->varFcn, this->domain, this->V);

  LS = new LevelSet<dimLS>(ioData, this->domain);
  riemann = new DistExactRiemannSolver<dim>(ioData,this->domain,this->varFcn);

  frequencyLS = ioData.mf.frequency;
  interfaceType = ioData.mf.interfaceType;

  Prate = ioData.implosion.Prate;
  Pinit = ioData.implosion.Pinit;
  tmax = (ioData.bc.inlet.pressure - Pinit)/Prate;

  requireSpecialBDF = false;

  int numBurnableFluids = ProgrammedBurn::countBurnableFluids(ioData);
  //std::cout << "Num burnable fluids = " << numBurnableFluids << std::endl;
  if (numBurnableFluids > 0) {
    programmedBurn = new ProgrammedBurn(ioData,this->X);
    this->fluidSelector.attachProgrammedBurn(programmedBurn);
  }
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
LevelSetTsDesc<dim,dimLS>::~LevelSetTsDesc()
{

  if (LS) delete LS;
  if (riemann) delete riemann;
}

//------------------------------------------------------------------------------
template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::setupTimeStepping(DistSVec<double,dim> *U, IoData &ioData)
{

  this->geoState->setup2(this->timeState->getData());

  // initalize solution
  this->timeState->setup(this->input->solutions, *this->X, this->bcData->getInletBoundaryVector(), *U, ioData);
  LS->setup(this->input->levelsets, *this->X, *U, Phi, ioData,&fluidSelector,this->varFcn);
  fluidSelector.initializeFluidIds(Phi, LS->Phinm1, LS->Phinm2); //setup fluidId in fluidSelector

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(this->mmh);
  if (_mmh) 
    _mmh->setup(&this->restart->frequency, &this->data->maxTime, this->postOp, *this->X, *U, fluidSelector.fluidId);
 
  *this->Xs = *this->X;

  this->timer->setSetupTime();
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
double LevelSetTsDesc<dim,dimLS>::computeTimeStep(int it, double *dtLeft,
                                            DistSVec<double,dim> &U)
{
  double t0 = this->timer->getTime();

  this->data->computeCflNumber(it - 1, this->data->residual / this->restart->residual);

  umax = 0.0;
  int numSubCycles = 1;
  double dt = this->timeState->computeTimeStep(this->data->cfl, dtLeft,
                            &numSubCycles, *this->geoState, *this->A, U, *(fluidSelector.fluidId),&umax);

  if (this->problemType[ProblemData::UNSTEADY])
    this->com->printf(5, "Global dt: %g (remaining subcycles = %d)\n",
                      dt*this->refVal->time, numSubCycles);

  this->timer->addFluidSolutionTime(t0);
  this->timer->addTimeStepTime(t0);

  return dt;

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::updateStateVectors(DistSVec<double,dim> &U, int it)
{

  this->geoState->update(*this->X, *this->A);

  if(frequencyLS > 0 && it%frequencyLS == 0){
//    this->com->printf(5, "LevelSet norm before reinitialization = %e\n", Phi.norm());
    LS->conservativeToPrimitive(Phi,PhiV,U);
    LS->reinitializeLevelSet(*this->X, PhiV);
    LS->primitiveToConservative(PhiV,Phi,U);
//    this->com->printf(5, "LevelSet norm after reinitialization = %e\n", Phi.norm());
    LS->update(Phi);

    // If we are doing 3BDF, reinitialization destroys unm1
    // Create a new version using F(U) = dU/dt
    if (this->timeState->useNm1()) {
      DistSVec<double,dimLS>& Phinm1 = LS->getPhinm1();
      this->multiPhaseSpaceOp->computeResidualLS(*this->X, *this->A, Phi, *fluidSelector.fluidId, U, Phinm1);
      Phinm1 = -1.0*Phinm1;
      requireSpecialBDF = true;
    }      
  } else {
    requireSpecialBDF = false;
    LS->update(Phi);
  }

  if (programmedBurn) {

    programmedBurn->setFluidIds(currentTime, *fluidSelector.fluidId,U);
  }

  fluidSelector.update();
 
  this->timeState->update(U,Utilde,  *(fluidSelector.fluidIdn), fluidSelector.fluidIdnm1, riemann);


}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
int LevelSetTsDesc<dim,dimLS>::checkSolution(DistSVec<double,dim> &U)
{

  int ierr = this->domain->checkSolution(this->varFcn, *this->A, U, *fluidSelector.fluidId, *fluidSelector.fluidIdn);

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::setupOutputToDisk(IoData &ioData, bool *lastIt,
                          int it, double t, DistSVec<double,dim> &U)
{
  if (it == this->data->maxIts)
    *lastIt = true;
  else
    monitorInitialState(it, U); // Phi?

  this->output->setMeshMotionHandler(ioData, this->mmh);
  this->output->openAsciiFiles();
  this->output->cleanProbesFile();

  this->timer->setSetupTime();

  if (it == 0) {
    this->output->writeForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, fluidSelector.fluidId);
    this->output->writeLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, fluidSelector.fluidId);
    this->output->writeHydroForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, fluidSelector.fluidId);
    this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, fluidSelector.fluidId);
    this->output->writeResidualsToDisk(it, 0.0, 1.0, this->data->cfl);
    this->output->writeMaterialVolumesToDisk(it, 0.0, *this->A, fluidSelector.fluidId);
    this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState, *fluidSelector.fluidId,&Phi);
    this->output->writeHeatFluxesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, fluidSelector.fluidId);
  }

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::outputToDisk(IoData &ioData, bool* lastIt, int it,
                                               int itSc, int itNl,
                                               double t, double dt, DistSVec<double,dim> &U)
{
  this->com->globalSum(1, &interruptCode);
  if (interruptCode)
    *lastIt = true;

  double cpu = this->timer->getRunTime();
  double res = this->data->residual / this->restart->residual;
                                                                                                      
  this->output->writeLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, fluidSelector.fluidId);
  this->output->writeHydroForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, fluidSelector.fluidId);
  this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, fluidSelector.fluidId);
  this->output->writeResidualsToDisk(it, cpu, res, this->data->cfl);
  this->output->writeMaterialVolumesToDisk(it, t, *this->A, fluidSelector.fluidId);
  this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState,*fluidSelector.fluidId,&Phi);
  this->output->writeProbesToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState,*fluidSelector.fluidId,&Phi);
  this->restart->writeToDisk(this->com->cpuNum(), *lastIt, it, t, dt, *this->timeState, *this->geoState, LS);

  this->output->updatePrtout(t);
  this->restart->updatePrtout(t);
  if (*lastIt) {
    this->timer->setRunTime();
    if (this->com->getMaxVerbose() >= 2)
      this->timer->print(this->domain->getStrTimer());
    this->output->closeAsciiFiles();
  }

}
//------------------------------------------------------------------------------

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::outputForces(IoData &ioData, bool* lastIt, int it, int itSc, int itNl,
                               double t, double dt, DistSVec<double,dim> &U)  {

  double cpu = this->timer->getRunTime();
  this->output->writeForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, fluidSelector.fluidId);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::resetOutputToStructure(DistSVec<double,dim> &U)
{
  this->com->printf(5,"LevelSetTsDesc<dim,dimLS>::resetOutputToStructure\n");

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(this->mmh);
  if (_mmh) 
    _mmh->resetOutputToStructure(this->postOp, *this->X, U, fluidSelector.fluidId);

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::updateOutputToStructure(double dt, double dtLeft,
					  DistSVec<double,dim> &U)
{

  this->com->printf(5,"LevelSetTsDesc<dim,dimLS>::resetOutputToStructure\n");
  if (this->mmh) {
    double work[2];
    this->mmh->computeInterfaceWork(dt, this->postOp, this->geoState->getXn(), 
                                    this->timeState->getUn(), *this->X, U, 
                                    work, fluidSelector.fluidIdn, fluidSelector.fluidId);
    this->restart->energy[0] += work[0];
    this->restart->energy[1] += work[1];
  }

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(this->mmh);
  if (_mmh)
    _mmh->updateOutputToStructure(dt, dtLeft, this->postOp, *this->X, U, fluidSelector.fluidId);

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
bool LevelSetTsDesc<dim,dimLS>::IncreasePressure(double dt, double t, DistSVec<double,dim> &U)
{
  if(Pinit<0.0 || Prate<0.0) return true; // no setup for increasing pressure

  if(t>tmax && t-dt>tmax) {this->com->fprintf(stdout, "max pressure reached\n"); return true;} // max pressure was reached, so now we solve
  else{ // max pressure not reached, so we do not solve and we increase pressure and let structure react
    this->com->fprintf(stdout, "about to increase pressure to value of %e\n", Pinit+t*Prate);
    this->domain->IncreasePressure(Pinit+t*Prate, this->varFcn, U);
    return false;
  }


}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::avoidNewPhaseCreation(DistSVec<double,dimLS> &localPhi)
{

  this->domain->avoidNewPhaseCreation(localPhi, LS->Phin);

}

//------------------------------------------------------------------------------

template<int dim,int dimLS>
void LevelSetTsDesc<dim,dimLS>::fixSolution(DistSVec<double,dim>& U,DistSVec<double,dim>& dU) {

  if (this->fixSol == 1)
    this->domain->fixSolution(this->varFcn,U,dU,fluidSelector.fluidId);
}


template<int dim,int dimLS>
void LevelSetTsDesc<dim,dimLS>::setCurrentTime(double t,DistSVec<double,dim>& U) { 

  currentTime = t;

  if (programmedBurn)
    programmedBurn->setCurrentTime(t,multiPhaseSpaceOp->getVarFcn(), U,*(fluidSelector.fluidId),*(fluidSelector.fluidIdn));
}
