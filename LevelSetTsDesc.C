#include <LevelSetTsDesc.h>

#include "IoData.h"
#include "GeoSource.h"
#include "Domain.h"
#include "LevelSet.h"
#include <DistExactRiemannSolver.h>
#include <FluidSelector.h>

#include <math.h>
                                                                                                        
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
  PhiV(this->getVecInfo()), boundaryFlux(this->getVecInfo()),
  computedQty(this->getVecInfo()), interfaceFlux(this->getVecInfo()),
  fluidSelector(ioData.eqs.numPhase, ioData, dom) 
{

  multiPhaseSpaceOp = new MultiPhaseSpaceOperator<dim,dimLS>(ioData, this->varFcn, this->bcData, this->geoState, 
                                                             this->domain, this->V);
  this->timeState = new DistTimeState<dim>(ioData, multiPhaseSpaceOp, this->varFcn, this->domain, this->V);

  LS = new LevelSet<dimLS>(ioData, this->domain);
  riemann = new DistExactRiemannSolver<dim>(ioData,this->domain,this->varFcn);

  //multiphase conservation check
  boundaryFlux  = 0.0;
  interfaceFlux = 0.0;
  computedQty   = 0.0;
  tmpDistSVec   = 0;
  tmpDistSVec2  = 0;
  expected = new double * [dimLS+2];
  computed = new double * [dimLS+2];
  for(int j=0; j<dimLS+2; j++){
    expected[j] = new double[dim];
    computed[j] = new double[dim];
    for(int i=0; i<dim; i++){
      expected[j][i] = 0.0;
      computed[j][i] = 0.0;
    }
  }

  frequencyLS = ioData.mf.frequency;
  interfaceType = ioData.mf.interfaceType;

  Prate = ioData.mf.Prate;
  Pinit = ioData.mf.Pinit;
  tmax = (ioData.bc.inlet.pressure - Pinit)/Prate;

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
LevelSetTsDesc<dim,dimLS>::~LevelSetTsDesc()
{

  if (LS) delete LS;
  if (riemann) delete riemann;

  if(tmpDistSVec)  delete tmpDistSVec;
  if(tmpDistSVec2) delete tmpDistSVec2;
  for(int j=0; j<dimLS+2; j++){
    delete expected[j];
    delete computed[j];
  }
  delete expected;
  delete computed;

}

//------------------------------------------------------------------------------
template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::setupTimeStepping(DistSVec<double,dim> *U, IoData &ioData)
{

  this->geoState->setup2(this->timeState->getData());

  // initalize solution
  this->timeState->setup(this->input->solutions, *this->X, this->bcData->getInletBoundaryVector(), *U, ioData);
  LS->setup(this->input->levelsets, *this->X, *U, Phi, ioData);
  fluidSelector.initializeFluidIds(Phi, LS->Phinm1, LS->Phinm2); //setup fluidId in fluidSelector

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(this->mmh);
  if (_mmh) 
    _mmh->setup(&this->restart->frequency, &this->data->maxTime, this->postOp, *this->X, *U, &(fluidSelector.fluidId));
  
  *this->Xs = *this->X;

  //this->timer->setSetupTime();
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
double LevelSetTsDesc<dim,dimLS>::computeTimeStep(int it, double *dtLeft,
                                            DistSVec<double,dim> &U)
{
  double t0 = this->timer->getTime();

  this->data->computeCflNumber(it - 1, this->data->residual / this->restart->residual);

  int numSubCycles = 1;
  double dt = this->timeState->computeTimeStep(this->data->cfl, dtLeft,
                            &numSubCycles, *this->geoState, *this->A, U, fluidSelector.fluidId);

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

  LS->update(Phi);
  fluidSelector.update();
  
  this->timeState->update(U, fluidSelector.fluidIdn, fluidSelector.fluidIdnm1, riemann);

  if(frequencyLS > 0 && it%frequencyLS == 0){
    LS->conservativeToPrimitive(Phi,PhiV,U);
    LS->reinitializeLevelSet(*this->geoState,*this->X, *this->A, U, PhiV);
    LS->primitiveToConservative(PhiV,Phi,U);
  }

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
int LevelSetTsDesc<dim,dimLS>::checkSolution(DistSVec<double,dim> &U)
{

  int ierr = this->domain->checkSolution(this->varFcn, *this->A, U, fluidSelector.fluidId, fluidSelector.fluidIdn);

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::conservationErrors(DistSVec<double,dim> &U, int it)
{
  if(!tmpDistSVec)  tmpDistSVec  = new DistSVec<double,dim>(this->getVecInfo());
// computes the total mass, momentum and energy
// 1- in the whole domain regardless of which phase they belong to
// 2- in the positive phase (phi>=0.0 <=> sign= 1)
// 3- in the negative phase (phi< 0.0 <=> sign=-1)
// The computed total domain mass, momentum and energy can be compared
//     to expected values (since only what goes out and comes in need to be
//     to be accounted for).
// When considering one phase, only the mass can be compared to an expected
//     value since the flux is zero at the interface between the two fluids.
//     For the momentum and the energy, the physical flux depends on the 
//     pressure at this interface. More to be said here (numerical flux among
//     other things)

  // total mass, total mass of fluid1, total mass of fluid2
  // note that there are two types of conservation errors:
  // 1 - fluxes do not cancel at the interface
  // 2 - populating of cell that changes fluid
  // Both have an effect on the total mass
  // Only the type-2 error has an effect on the total mass of fluid-i

  computedQty = U;
  computedQty *= *this->A;
  computedQty.sum(computed[0]);

  for(int i=1; i<dimLS+2; i++){ // for each separate fluid
    int fluidIdTarget = i-1; //use fluidSelector to obtain them?
    this->domain->restrictionOnPhi(computedQty, fluidSelector.fluidId, *tmpDistSVec, fluidIdTarget); //tmpDistSVec reinitialized to 0.0 inside routine
    tmpDistSVec->sum(computed[i]);
  }

  // expected total mass, mass in fluid1, mass in fluid2
  // an expected mass is computed iteratively, that is
  // m_{n+1} = m_{n} + fluxes_at_boundaries
  if(it==0){
    for (int j=0; j<dimLS+2; j++)
      for (int i=0; i<dim; i++)
        expected[j][i] = computed[j][i];
    return;
  }

  const double dt = this->timeState->getTime();
  double bcfluxsum[dim];

  boundaryFlux.sum(bcfluxsum);
  for (int i=0; i<dim; i++)
    expected[0][i] += dt*bcfluxsum[i];

  for(int j=1; j<dimLS+2; j++){
    int fluidIdTarget = j-1; //use fluidSelector to obtain them?
    this->domain->restrictionOnPhi(boundaryFlux, fluidSelector.fluidId, *tmpDistSVec, fluidIdTarget); //tmpDistSVec reinitialized to 0.0 inside routine
    tmpDistSVec->sum(bcfluxsum);
    for (int i=0; i<dim; i++)
      expected[j][i] += dt*bcfluxsum[i];
  }

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::setupOutputToDisk(IoData &ioData, bool *lastIt,
                          int it, double t, DistSVec<double,dim> &U)
{
  conservationErrors(U,it);
 
  if (it == this->data->maxIts)
    *lastIt = true;
  else
    monitorInitialState(it, U); // Phi?

  this->output->setMeshMotionHandler(ioData, this->mmh);

  this->output->openAsciiFiles();

  if (it == 0) {
    this->output->writeForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &(fluidSelector.fluidId));
    this->output->writeLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &(fluidSelector.fluidId));
    this->output->writeHydroForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &(fluidSelector.fluidId));
    this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &(fluidSelector.fluidId));
    this->output->writeResidualsToDisk(it, 0.0, 1.0, this->data->cfl);
    this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState, Phi, (fluidSelector.fluidId));
    this->output->writeConservationErrors(ioData, it, t, dimLS+2, expected, computed);
    this->output->writeHeatFluxesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &(fluidSelector.fluidId));
  }

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::outputToDisk(IoData &ioData, bool* lastIt, int it,
                                               int itSc, int itNl,
                                               double t, double dt, DistSVec<double,dim> &U)
{
  conservationErrors(U,it);

  this->com->globalSum(1, &interruptCode);
  if (interruptCode)
    *lastIt = true;

  double cpu = this->timer->getRunTime();
  double res = this->data->residual / this->restart->residual;
                                                                                                      
  this->output->writeLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &(fluidSelector.fluidId));
  this->output->writeHydroForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &(fluidSelector.fluidId));
  this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &(fluidSelector.fluidId));
  this->output->writeResidualsToDisk(it, cpu, res, this->data->cfl);
  this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState, Phi, fluidSelector.fluidId);
  this->output->writeConservationErrors(ioData, it, t, dimLS+2, expected, computed);
  this->restart->writeToDisk(this->com->cpuNum(), *lastIt, it, t, dt, *this->timeState, *this->geoState, LS);

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
  this->output->writeForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &(fluidSelector.fluidId));
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::resetOutputToStructure(DistSVec<double,dim> &U)
{
  this->com->printf(5,"LevelSetTsDesc<dim,dimLS>::resetOutputToStructure\n");

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(this->mmh);
  if (_mmh) 
    _mmh->resetOutputToStructure(this->postOp, *this->X, U, &(fluidSelector.fluidId));

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
                                    work, &(fluidSelector.fluidIdn), &(fluidSelector.fluidId));
    this->restart->energy[0] += work[0];
    this->restart->energy[1] += work[1];
  }

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(this->mmh);
  if (_mmh)
    _mmh->updateOutputToStructure(dt, dtLeft, this->postOp, *this->X, U, &(fluidSelector.fluidId));

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
