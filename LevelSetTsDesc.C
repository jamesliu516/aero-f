#include <LevelSetTsDesc.h>
#include <DistExactRiemannSolver.h>


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

template<int dim>
LevelSetTsDesc<dim>::
LevelSetTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  TsDesc<dim>(ioData, geoSource, dom), Phi(this->getVecInfo()), Vg(this->getVecInfo()),
	PhiV(this->getVecInfo())
{

  this->timeState = new DistTimeState<dim>(ioData, this->spaceOp, this->varFcn, this->domain, this->V);

  LS = new LevelSet(ioData, this->domain);
  riemann = new DistExactRiemannSolver<dim>(ioData,this->domain);

  Vgf = 0;
  Vgfweight = 0;
  if(ioData.mf.typePhaseChange == MultiFluidData::EXTRAPOLATION){
    Vgf = new DistSVec<double,dim>(this->getVecInfo());
    Vgfweight = new DistVec<double>(this->getVecInfo());
    *Vgf =-1.0;
    *Vgfweight =0.0;
  }

  frequencyLS = ioData.mf.frequency;
  interfaceType = ioData.mf.interfaceType;
}

//------------------------------------------------------------------------------

template<int dim>
LevelSetTsDesc<dim>::~LevelSetTsDesc()
{

  if (LS) delete LS;
  if (riemann) delete riemann;
  if (Vgf) delete Vgf;
  if (Vgfweight) delete Vgfweight;

}

//------------------------------------------------------------------------------
template<int dim>
void LevelSetTsDesc<dim>::setupTimeStepping(DistSVec<double,dim> *U, IoData &ioData)
{

  this->geoState->setup2(this->timeState->getData());


  // initalize solution
  fprintf(stdout, "init0\n");
  this->timeState->setup(this->input->solutions, *this->X, this->bcData->getInletBoundaryVector(), *U, ioData);
  fprintf(stdout, "init1\n");
  LS->setup(this->input->levelsets, *this->X, *U, Phi, ioData);
  fprintf(stdout, "init2\n");

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(this->mmh);
  if (_mmh)
    _mmh->setup(&this->restart->frequency, &this->data->maxTime, this->postOp, *this->X, *U, &Phi);

  *this->Xs = *this->X;

  //this->timer->setSetupTime();
}

//------------------------------------------------------------------------------

template<int dim>
double LevelSetTsDesc<dim>::computeTimeStep(int it, double *dtLeft,
                                            DistSVec<double,dim> &U)
{
  double t0 = this->timer->getTime();

  this->data->computeCflNumber(it - 1, this->data->residual / this->restart->residual);

  int numSubCycles = 1;
  double dt = this->timeState->computeTimeStep(this->data->cfl, dtLeft,
                            &numSubCycles, *this->geoState, *this->A, U, Phi);

  if (this->problemType[ProblemData::UNSTEADY])
    this->com->printf(5, "Global dt: %g (remaining subcycles = %d)\n",
                      dt*this->refVal->time, numSubCycles);

  this->timer->addFluidSolutionTime(t0);
  this->timer->addTimeStepTime(t0);

  return dt;

}

//------------------------------------------------------------------------------

template<int dim>
void LevelSetTsDesc<dim>::updateStateVectors(DistSVec<double,dim> &U, int it)
{

  this->geoState->update(*this->X, *this->A);

  LS->update(Phi);

  this->timeState->update(U, LS->Phin, LS->Phinm1, LS->Phinm2,
                          Vgf, Vgfweight, riemann);

  if(frequencyLS > 0 && it%frequencyLS == 0){
    LS->conservativeToPrimitive(Phi,PhiV,U);
    LS->reinitializeLevelSet(*this->geoState,*this->X, *this->A, U, PhiV);
    LS->primitiveToConservative(PhiV,Phi,U);
  }

}

//------------------------------------------------------------------------------

template<int dim>
int LevelSetTsDesc<dim>::checkSolution(DistSVec<double,dim> &U)
{

  int ierr = this->domain->checkSolution(this->varFcn, *this->A, U, Phi);

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim>
void LevelSetTsDesc<dim>::setupOutputToDisk(IoData &ioData, bool *lastIt,
                                                    int it, double t, DistSVec<double,dim> &U)
{
 
  if (it == this->data->maxIts)
    *lastIt = true;
  else
    monitorInitialState(it, U); // Phi?

  this->output->setMeshMotionHandler(ioData, this->mmh);

  this->output->openAsciiFiles();

  if (it == 0) {
    this->output->writeForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &Phi);
    this->output->writeLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &Phi);
    this->output->writeHydroForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &Phi);
    this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &Phi);
    this->output->writeResidualsToDisk(it, 0.0, 1.0, this->data->cfl);
    this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, Phi);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void LevelSetTsDesc<dim>::outputToDisk(IoData &ioData, bool* lastIt, int it,
                                               int itSc, int itNl,
                                               double t, double dt, DistSVec<double,dim> &U)
{

  this->com->globalSum(1, &interruptCode);
  if (interruptCode)
    *lastIt = true;

  double cpu = this->timer->getRunTime();
  double res = this->data->residual / this->restart->residual;
                                                                                                      
  this->output->writeLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &Phi);
  this->output->writeHydroForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &Phi);
  this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &Phi);
  this->output->writeResidualsToDisk(it, cpu, res, this->data->cfl);
  this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, Phi);
  this->restart->writeToDisk(this->com->cpuNum(), *lastIt, it, t, dt, *this->timeState, *this->geoState, LS);

  if (*lastIt) {
    this->timer->setRunTime();
    if (this->com->getMaxVerbose() >= 2)
      this->timer->print(this->domain->getStrTimer());
    this->output->closeAsciiFiles();
  }

}
//------------------------------------------------------------------------------

template<int dim>
void LevelSetTsDesc<dim>::outputForces(IoData &ioData, bool* lastIt, int it, int itSc, int itNl,
                               double t, double dt, DistSVec<double,dim> &U)  {

  double cpu = this->timer->getRunTime();
  this->output->writeForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &Phi);
}

//------------------------------------------------------------------------------

template<int dim>
void LevelSetTsDesc<dim>::resetOutputToStructure(DistSVec<double,dim> &U)
{
  this->com->printf(5,"LevelSetTsDesc<dim>::resetOutputToStructure\n");

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(this->mmh);
  if (_mmh) 
    _mmh->resetOutputToStructure(this->postOp, *this->X, U, &Phi);

}

//------------------------------------------------------------------------------

template<int dim>
void LevelSetTsDesc<dim>::updateOutputToStructure(double dt, double dtLeft,
					  DistSVec<double,dim> &U)
{

  this->com->printf(5,"LevelSetTsDesc<dim>::resetOutputToStructure\n");
  if (this->mmh) {
    double work[2];
    this->mmh->computeInterfaceWork(dt, this->postOp, this->geoState->getXn(), 
                                    this->timeState->getUn(), *this->X, U, 
                                    work, &(LS->Phin), &Phi);
    this->restart->energy[0] += work[0];
    this->restart->energy[1] += work[1];
  }

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(this->mmh);
  if (_mmh)
    _mmh->updateOutputToStructure(dt, dtLeft, this->postOp, *this->X, U, &Phi);

}

//------------------------------------------------------------------------------
